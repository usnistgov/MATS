import numpy as np
from numba import jit, prange, float64, complex128
from math import tanh, log10, exp, sqrt, pi, cos, sin, atan2

# Import the industry standard CPF from your existing file
# Assuming numba_mHTP.py is in the same folder as CPF.py
from .CPF import cpf_accurate

@jit(float64(float64, float64, float64), nopython=True, fastmath=True, cache=True)
def beta_jit(GammaD, NuOptRe, alpha):
    """ JIT-compiled version of the Beta correction function """
    if alpha < 5.0:
        val = (0.0534 + 0.1585 * exp(-0.4510 * alpha)) * \
              tanh((1.9595 + alpha * (-0.1258 + alpha * (0.0056 + alpha * 0.0050))) * \
              log10(NuOptRe / GammaD) + \
              (-0.0546 + alpha * (0.0672 + alpha * (-0.0125 + alpha * 0.0003)))) + \
              0.9466 - 0.1585 * exp(-0.4510 * alpha)
        return val
    else:
        return 1.0

@jit(float64[:](float64, float64, float64, float64, float64, float64, 
                float64, float64, float64[:], float64, float64, float64, float64), 
     nopython=True, fastmath=True, parallel=False, cache=True)
def mHTprofile_vector_numba(nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, 
                            NuOptRe, NuOptIm, nu_grid, 
                            Ylm, Xlm, alpha, disp_val):
    """
    Numba-optimized mHTP. 
    Uses SCALAR loops to avoid allocating temporary arrays (X, Y, Z, etc.) 
    Resulting in massive speedups over the numpy-vectorized version.
    """
    
    # 1. Setup Constants
    nuD = 1.2011224087864498 * GammaD
    inv_nuD = 1.0 / nuD
    
    # Beta correction
    if alpha < 5.0 and NuOptRe != 0.0:
        nuR = NuOptRe * beta_jit(GammaD, NuOptRe, alpha)
    else:
        nuR = NuOptRe
        
    # Line Mixing: 1 - i*Y + X (Args were swapped in original profile.py matching?)
    # Original: LM = args[1] + 1.0 - args[0]*1j 
    # Mapped: args[1] -> Xlm, args[0] -> Ylm
    LM_real = 1.0 + Xlm
    LM_imag = -Ylm
    
    c2_real = Gamma2
    c2_imag = Delta2
    # c0 = Gamma0+Delta0*1j - 1.5*c2 + nuR + NuOptIm*1j
    c0_real = Gamma0 - 1.5 * Gamma2 + nuR
    c0_imag = Delta0 - 1.5 * Delta2 + NuOptIm
    
    n_points = len(nu_grid)
    result = np.zeros(n_points, dtype=np.float64)

    # 2. Case Selection based on c2 (Speed Dependence)
    abs_c2 = sqrt(c2_real**2 + c2_imag**2)
    
    if abs_c2 > 1.0e-9:
        # --- Speed Dependent Case ---
        denom_norm = c2_real**2 + c2_imag**2
        inv_c2_real = c2_real / denom_norm
        inv_c2_imag = -c2_imag / denom_norm
        
        # Y = 0.25 * (nuD/c2)^2
        # Calculate Y once (it is independent of nu)
        temp_real = nuD * inv_c2_real
        temp_imag = nuD * inv_c2_imag
        Y_real = 0.25 * (temp_real**2 - temp_imag**2)
        Y_imag = 0.25 * (2 * temp_real * temp_imag)
        abs_Y = sqrt(Y_real**2 + Y_imag**2)
        
        # Pre-calc sqrt(Y)
        if abs_Y > 0:
            r_Y = sqrt(abs_Y)
            theta_Y = atan2(Y_imag, Y_real) * 0.5
            csqY_real = r_Y * cos(theta_Y)
            csqY_imag = r_Y * sin(theta_Y)
        else:
            csqY_real = 0.0; csqY_imag = 0.0

        # --- Optimization: Estimate min(abs(X)) ---
        # The original code scans the whole array min(abs(X)). 
        # Analytically, X is smallest at line center (nu=nu0).
        # X_center = c0 / c2
        X_cen_real = c0_real * inv_c2_real - c0_imag * inv_c2_imag
        X_cen_imag = c0_real * inv_c2_imag + c0_imag * inv_c2_real
        vecabsX = sqrt(X_cen_real**2 + X_cen_imag**2)

        # Loop over grid
        for i in range(n_points):
            nu = nu_grid[i]
            
            # X = ((nu0 - nu)*1j + c0) / c2
            # num = c0 + i*(nu0-nu)
            num_real = c0_real
            num_imag = (nu0 - nu) + c0_imag
            
            X_real = num_real * inv_c2_real - num_imag * inv_c2_imag
            X_imag = num_real * inv_c2_imag + num_imag * inv_c2_real
            
            abs_X = sqrt(X_real**2 + X_imag**2)
            
            # Logic: if abs(Y) > vecabsX * 1e-15
            # Note: We use the pre-calculated vecabsX to match original logic
            if abs_Y > vecabsX * 1.0e-15:
                # z2 = sqrt(X+Y) + csqY
                sum_real = X_real + Y_real
                sum_imag = X_imag + Y_imag
                r_sum = sqrt(sqrt(sum_real**2 + sum_imag**2))
                theta_sum = atan2(sum_imag, sum_real) * 0.5
                
                z2_real = r_sum * cos(theta_sum) + csqY_real
                z2_imag = r_sum * sin(theta_sum) + csqY_imag
                
                if vecabsX > abs_Y * 3e-8:
                    z1_real = z2_real - 2 * csqY_real
                    z1_imag = z2_imag - 2 * csqY_imag
                else:
                    z1_real = num_real * inv_nuD
                    z1_imag = num_imag * inv_nuD
                
                # cpf takes (-z.imag, z.real)
                w1 = cpf_accurate(-z1_imag, z1_real)
                w2 = cpf_accurate(-z2_imag, z2_real)
                
                # A = 1.772... / nuD * (w1-w2)
                pre = 1.772453850905516 * inv_nuD
                A_real = pre * (w1.real - w2.real)
                A_imag = pre * (w1.imag - w2.imag)
                
            else:
                # Approximation for small Y
                r_X = sqrt(abs_X)
                theta_X = atan2(X_imag, X_real) * 0.5
                rX_real = r_X * cos(theta_X)
                rX_imag = r_X * sin(theta_X)
                
                if vecabsX < 4000.0:
                    wX = cpf_accurate(-rX_imag, rX_real)
                    
                    # Term: rX * wX
                    term_real = rX_real*wX.real - rX_imag*wX.imag
                    term_imag = rX_real*wX.imag + rX_imag*wX.real
                    
                    top_real = 2.0 - 3.5449077018110318 * term_real
                    top_imag = -3.5449077018110318 * term_imag
                    
                    # A = top / c2
                    A_real = top_real * inv_c2_real - top_imag * inv_c2_imag
                    A_imag = top_real * inv_c2_imag + top_imag * inv_c2_real
                else:
                    # Large X approx A = (1 - 1.5/X) / (X*c2)
                    # For safety in far wings or optimization, can set to 0 or Voigt limit
                    A_real = 0.0; A_imag = 0.0
            
            # --- Final Intensity ---
            # I = 0.318... * LM / (1/A - (nuR + i*NuOptIm))
            
            norm_A = A_real**2 + A_imag**2
            if norm_A > 0:
                invA_real = A_real / norm_A
                invA_imag = -A_imag / norm_A
                
                denom_real = invA_real - nuR
                denom_imag = invA_imag - NuOptIm
                
                norm_denom = denom_real**2 + denom_imag**2
                if norm_denom > 0:
                    res_real = (LM_real * denom_real + LM_imag * denom_imag) / norm_denom
                    res_imag = (LM_imag * denom_real - LM_real * denom_imag) / norm_denom
                    
                    factor = 0.3183098861837907
                    if disp_val > 0.5:
                        result[i] = -res_imag * factor
                    else:
                        result[i] = res_real * factor

    else:
        # --- Speed Independent Case (Standard Voigt/Rautian) ---
        for i in range(n_points):
            nu = nu_grid[i]
            z_real = c0_real * inv_nuD
            z_imag = ((nu0 - nu) + c0_imag) * inv_nuD
            
            w = cpf_accurate(-z_imag, z_real)
            
            pre = 1.772453850905516 * inv_nuD
            A_real = pre * w.real
            A_imag = pre * w.imag
            
            norm_A = A_real**2 + A_imag**2
            if norm_A > 0:
                invA_real = A_real / norm_A
                invA_imag = -A_imag / norm_A
                
                denom_real = invA_real - nuR
                denom_imag = invA_imag - NuOptIm
                
                norm_denom = denom_real**2 + denom_imag**2
                if norm_denom > 0:
                    res_real = (LM_real * denom_real + LM_imag * denom_imag) / norm_denom
                    res_imag = (LM_imag * denom_real - LM_real * denom_imag) / norm_denom
                    
                    factor = 0.3183098861837907
                    if disp_val > 0.5:
                        result[i] = -res_imag * factor
                    else:
                        result[i] = res_real * factor
                
    return result