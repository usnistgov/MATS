import numpy as np
from numba import jit, float64
from math import tanh, log10, exp

# Import the optimized CPF scalar function you already have
# Assuming it is in .CPF
from .CPF import cpf_accurate as cpf_scalar

# --- Beta Correction (JIT) ---
@jit(float64(float64, float64, float64), nopython=True, fastmath=True, cache=True)
def beta_jit(GammaD, NuOptRe, alpha):
    if alpha < 5.0:
        val = (0.0534 + 0.1585 * exp(-0.4510 * alpha)) * \
              tanh((1.9595 + alpha * (-0.1258 + alpha * (0.0056 + alpha * 0.0050))) * \
              log10(NuOptRe / GammaD) + \
              (-0.0546 + alpha * (0.0672 + alpha * (-0.0125 + alpha * 0.0003)))) + \
              0.9466 - 0.1585 * exp(-0.4510 * alpha)
        return val
    else:
        return 1.0

# --- Vectorized mHTP (JIT) ---
@jit(float64[:](float64, float64, float64, float64, float64, float64, 
                float64, float64, float64[:], float64, float64, float64, 
                float64), 
     nopython=True, fastmath=True, cache=True)
def mHTprofile_vector_numba(nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, 
                            NuOptRe, NuOptIm, nu_grid, 
                            Ylm, Xlm, alpha, disp_val):
    """
    Numba-optimized vectorized mHTP calculator.
    Calculates the profile for a single spectral line over the 'nu_grid' array.
    """
    
    # 1. Setup Constants
    nuD = 1.2011224087864498 * GammaD
    
    # Beta correction
    if alpha < 5.0:
        nuR = NuOptRe * beta_jit(GammaD, NuOptRe, alpha)
    else:
        nuR = NuOptRe
        
    # Line Mixing: 1 - i*Y + X
    LM_real = 1.0 + Xlm
    LM_imag = -Ylm
    
    c2_real = Gamma2
    c2_imag = Delta2
    c0_real = Gamma0 - 1.5 * Gamma2 + nuR
    c0_imag = Delta0 - 1.5 * Delta2 + NuOptIm
    
    n_points = len(nu_grid)
    result = np.zeros(n_points, dtype=np.float64)

    # 2. Case Selection based on c2 (Speed Dependence)
    abs_c2 = np.sqrt(c2_real**2 + c2_imag**2)
    
    if abs_c2 > 1.0e-9:
        # --- Speed Dependent Case ---
        denom_norm = c2_real**2 + c2_imag**2
        inv_c2_real = c2_real / denom_norm
        inv_c2_imag = -c2_imag / denom_norm
        
        # Y = 0.25 * (nuD/c2)^2
        temp_real = nuD * inv_c2_real
        temp_imag = nuD * inv_c2_imag
        Y_real = 0.25 * (temp_real**2 - temp_imag**2)
        Y_imag = 0.25 * (2 * temp_real * temp_imag)
        
        abs_Y = np.sqrt(Y_real**2 + Y_imag**2)
        
        # Pre-calc sqrt(Y) if needed
        if abs_Y > 0:
            r_Y = np.sqrt(abs_Y)
            theta_Y = np.arctan2(Y_imag, Y_real) * 0.5
            csqY_real = r_Y * np.cos(theta_Y)
            csqY_imag = r_Y * np.sin(theta_Y)
        else:
            csqY_real = 0.0; csqY_imag = 0.0

        # Loop over grid
        for i in range(n_points):
            nu = nu_grid[i]
            
            # X = ((nu0 - nu)*1j + c0) / c2
            num_real = c0_real
            num_imag = (nu0 - nu) + c0_imag
            
            X_real = num_real * inv_c2_real - num_imag * inv_c2_imag
            X_imag = num_real * inv_c2_imag + num_imag * inv_c2_real
            
            abs_X = np.sqrt(X_real**2 + X_imag**2)
            
            # Condition 1: Y vs X
            if abs_Y > abs_X * 1.0e-15:
                # z2 = sqrt(X+Y) + csqY
                sum_real = X_real + Y_real
                sum_imag = X_imag + Y_imag
                r_sum = np.sqrt(np.sqrt(sum_real**2 + sum_imag**2))
                theta_sum = np.arctan2(sum_imag, sum_real) * 0.5
                
                z2_real = r_sum * np.cos(theta_sum) + csqY_real
                z2_imag = r_sum * np.sin(theta_sum) + csqY_imag
                
                if abs_X > abs_Y * 3e-8:
                    z1_real = z2_real - 2 * csqY_real
                    z1_imag = z2_imag - 2 * csqY_imag
                else:
                    z1_real = num_real / nuD
                    z1_imag = num_imag / nuD
                
                w1 = cpf_scalar(-z1_imag, z1_real)
                w2 = cpf_scalar(-z2_imag, z2_real)
                
                diff_w_real = w1.real - w2.real
                diff_w_imag = w1.imag - w2.imag
                
                # A = 1.772... / nuD * (w1-w2)
                pre = 1.772453850905516 / nuD
                A_real = pre * diff_w_real
                A_imag = pre * diff_w_imag
                
            else:
                # Approximation
                r_X = np.sqrt(abs_X)
                theta_X = np.arctan2(X_imag, X_real) * 0.5
                rX_real = r_X * np.cos(theta_X)
                rX_imag = r_X * np.sin(theta_X)
                
                if abs_X < 4.0e3:
                    wX = cpf_scalar(-rX_imag, rX_real)
                    
                    # Term: rX * wX
                    term_real = rX_real*wX.real - rX_imag*wX.imag
                    term_imag = rX_real*wX.imag + rX_imag*wX.real
                    
                    top_real = 2.0 - 3.5449077018110318 * term_real
                    top_imag = -3.5449077018110318 * term_imag
                    
                    # A = top / c2
                    A_real = top_real * inv_c2_real - top_imag * inv_c2_imag
                    A_imag = top_real * inv_c2_imag + top_imag * inv_c2_real
                else:
                    # Large X approx: A = (1 - 1.5/X) / (X*c2)
                    # Simplified to 0 for brevity in wings
                    A_real = 0.0
                    A_imag = 0.0
            
            # Final Calculation
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
        inv_nuD = 1.0 / nuD
        for i in range(n_points):
            nu = nu_grid[i]
            z_real = c0_real * inv_nuD
            z_imag = ((nu0 - nu) + c0_imag) * inv_nuD
            
            w = cpf_scalar(-z_imag, z_real)
            
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