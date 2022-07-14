#Import Packages
# from .Utilities import *
#import re


class Generate_FitParam_File:
    """Class generates the parameter files used in fitting and sets up fit settings.

    Parameters
    ----------
    dataset : object
        Dataset object to be used in the fits
    param_linelist : dataframe
        parameter linelist dataframe name, where the dataframe has the appropriate column headers
    base_linelist : dataframe
        baseline parameter dataframe name generated from the dataset.generate_baseline_paramlist() function.
    CIA_linelist : dataframe, optional
        Future Function: CIA linelist dataframe name generated from the dataset.generate_CIA_paramlist() function.
    lineprofile : str
        lineprofile to use for the simulation. This sets values in the parameter linelist to 0 and forces them to not vary unless manually changed. Default values is VP, voigt profile. Options are VP, SDVP, NGP, SDNGP, HTP.
    linemixing : bool
        If False, then all linemixing parameters are set to 0 and not able to float.  Default is False.
    threshold_intensity : float
        This is the minimum line intensity that will be simulated. Default value is 1e-30.
    fit_intensity  : float
        This is the minimum line intensity for which parameters will be set to float. This can be overwritten manually. Default value is 1e-26.
    fit_window : float
        currently not used
    sim_window : float
        This is the region outside of the wavelength region of the dataset that will be simulated for analysis. Default value is 5 cm-1.
    base_linelist_savename : str
        filename that the baseline linelist will be saved as. Default is Baseline_LineList
    param_linelist_savename : str
        filename that the parameter linelist will be saved as. Default is Parameter_LineList.
    CIA_linelist_save_name : str
        Future Feature: filename that the CIA linelist will be saved as. Default is CIA_LineList.
    nu_constrain : bool
        True indicates that the line centers will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    sw_constrain : bool
        True indicates that the line intensities will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    gamma0_constrain : bool
        True indicates that the collisional width will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    delta0_constrain : bool
        True indicates that the shift will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    aw_constrain : bool
        True indicates that the speed dependent width will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    as_constrain : bool
        True indicates that the speed dependent shift will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    nuVC_constrain : bool
        True indicates that the dicke narrowing term will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    eta_constrain : bool
        True indicates that the correlation parameter will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    linemixing_constrain : bool
        True indicates that the first-order linemixing term will be a global variable across all spectra. False generates a new value for each spectrum in the datas
    """

    def __init__ (self, dataset, param_linelist, base_linelist, CIA_linelist = None,
                  lineprofile = 'VP', linemixing = False, threshold_intensity = 1e-30, fit_intensity = 1e-26, fit_window = 1.5, sim_window = 5,
                  param_linelist_savename = 'Parameter_LineList', base_linelist_savename = 'Baseline_LineList', CIA_linelist_savename = 'CIA_LineList',
                 nu_constrain = True, sw_constrain = True, gamma0_constrain = True, delta0_constrain = True, aw_constrain = True, as_constrain = True,
                 nuVC_constrain = True, eta_constrain =True, linemixing_constrain = True,
                 additional_columns = []):
        self.dataset = dataset
        self.param_linelist = param_linelist
        self.base_linelist = base_linelist
        if self.dataset.CIA_model == None:
            self.CIA_linelist = None
            self.CIA_linelist_savename = None
        else:
            self.CIA_linelist = None
            self.CIA_linelist_savename = None

        self.lineprofile = lineprofile
        self.linemixing = linemixing
        self.threshold_intensity = threshold_intensity
        self.fit_intensity = fit_intensity
        self.fit_window = fit_window
        self.sim_window = sim_window
        self.base_linelist_savename = base_linelist_savename
        self.param_linelist_savename = param_linelist_savename
        self.nu_constrain = nu_constrain
        self.sw_constrain = sw_constrain
        self.gamma0_constrain = gamma0_constrain
        self.delta0_constrain = delta0_constrain
        self.aw_constrain = aw_constrain
        self.as_constrain = as_constrain
        self.nuVC_constrain = nuVC_constrain
        self.eta_constrain = eta_constrain
        self.linemixing_constrain = linemixing_constrain
        self.additional_columns = additional_columns
    def get_dataset(self):
        return self.dataset
    def get_param_linelist(self):
        return self.param_linelist
    def get_base_linelist(self):
        return self.base_linelist
    def get_CIA_linelist(self):
        return self.CIA_linelist
    def generate_fit_param_linelist_from_linelist(self, vary_nu = {}, vary_sw = {},
                                   vary_gamma0 = {}, vary_n_gamma0 = {},
                                   vary_delta0 = {}, vary_n_delta0 = {},
                                   vary_aw = {}, vary_n_gamma2 = {},
                                   vary_as = {}, vary_n_delta2 = {},
                                   vary_nuVC = {}, vary_n_nuVC = {},
                                   vary_eta = {}, vary_linemixing = {}, vary_n_linemixing = {}):
        """Generates the parameter line list used in fitting and updates the fitting booleans to desired settings.



        Parameters
        ----------
        vary_nu : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope line centers should be floated.  Each key in the dictionary corresponds to a molecule, where the values are a dictionary corresponding to local_iso_ids as keys and boolean values acting as boolean flags. True means float the nu for that molecule and isotope. Example vary_nu = {7:{1:True, 2: False, 3: False}, 1:{1:False, 2: False, 3: False}}, would indicate that line centers of all main oxygen isotopes would be floated (if the line intensity is greater than the fit intensity), where the line centers for all other O2 isotopes and all water isotopes would not be varied. If these value are left blank, then all variable would have to be manually switched to float. The default is {}.
        vary_sw : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope line intensities should be floated.  Follows nu_vary example.  The default is {}.
        vary_gamma0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope collisional half-width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_gamma0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the half width  should be floated.  Follows nu_vary example.  . The default is {}.
        vary_delta0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope pressure shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_delta0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the pressure shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_aw : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope ratio of the speed dependent width to the half-width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_gamma2 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the speed dependent width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_as : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope the ratio of the speed dependent shift to the shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_delta2 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the speed dependent shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_nuVC : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope dicke narrowing should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_nuVC : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope  coefficient of the temperature dependence of the dicke narrowing should be floated.  Follows nu_vary example.  . The default is {}.
        vary_eta : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope correlation parameter for the VC and SD effects should be floated.  Follows nu_vary example.  . The default is {}.
        vary_linemixing : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope first-order line-mixing should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_linemixing : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope temperature dependence for the first-order line-mixing should be floated.  Follows nu_vary example.  . The default is {}.
        Returns
        -------
        param_linelist_df : dataframe
            returns dataframe based on parameter line list with addition of a vary and err column for every floatable parameter.  The vary columns are defined by the inputs and the fit_intensity value.  The err columns will be populated from fit results.  The dataframe is also saved as a .csv file.  line intensity will be normalized by the fit_intensity (set to the sw_scale_factor). The term 'sw' is now equal to the normalized value, such that in the simulation 'sw'*sw_scale_factor is used for the line intensity. Because line intensities are so small it is difficult to fit the value without normalization.

        """

        param_linelist_df = self.get_param_linelist().copy()
        #Intensity Thresholding
        param_linelist_df = param_linelist_df[param_linelist_df['sw'] > self.threshold_intensity] #intensity thresholding
        #Cutdown linelist to parameters within specified simulation window of minimum and maximum of dataset
        dataset_min, dataset_max = (self.dataset.get_spectra_extremes())
        extreme_dictionary = self.dataset.get_spectrum_extremes()
        param_linelist_df = param_linelist_df[param_linelist_df['nu'] < (dataset_max + self.sim_window)]
        param_linelist_df = param_linelist_df[param_linelist_df['nu'] > (dataset_min - self.sim_window)]
        #delete parameters not partaining to species (where relevant)
        diluent_list = []
        for spectrum in self.dataset.spectra:
            for diluent in spectrum.Diluent:
                if diluent not in diluent_list:
                    diluent_list.append(diluent)
        num_nominal_temps, list_nominal_temps = self.dataset.get_number_nominal_temperatures()
        column_list =  self.additional_columns.copy()
        column_list += ['molec_id', 'local_iso_id','elower', 'nu', 'sw']
        for diluent in diluent_list:
            column_list.append('gamma0_' + diluent)
            column_list.append('n_gamma0_' + diluent)
            column_list.append('delta0_' + diluent)
            column_list.append('n_delta0_' + diluent)
            column_list.append('SD_gamma_' + diluent)
            column_list.append('n_gamma2_' + diluent)
            column_list.append('SD_delta_' + diluent)
            column_list.append('n_delta2_' + diluent)
            column_list.append('nuVC_' + diluent)
            column_list.append('n_nuVC_' + diluent)
            column_list.append('eta_' + diluent)
            column_list.append('y_' + diluent)
            column_list.append('n_y_' + diluent)

        param_linelist_df = param_linelist_df[column_list]
        param_linelist_df = param_linelist_df.reset_index(drop = True)
        #Re-defines the Line intensity as sw*sw_scale_factor
        param_linelist_df['sw'] = param_linelist_df['sw'] / self.fit_intensity
        param_linelist_df['sw_scale_factor'] = [self.fit_intensity]*len(param_linelist_df)

        # Defines the Linecenter parameters in the event that the line center is held constant across all samples and is not
        ## Starting point is equal to the inital value
        order_nu = ['nu']
        param_linelist_df['nu_vary'] = len(param_linelist_df)*[False]
        param_linelist_df['nu_err'] = len(param_linelist_df)*[0]
        if self.nu_constrain:
            if vary_nu != {}:
                for molecule in vary_nu:
                    for isotope in vary_nu[molecule]:
                        param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) & (param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nu_vary'] = (vary_nu[molecule][isotope])
        else:
            for spec in self.dataset.get_list_spectrum_numbers():
                param_linelist_df['nu_' + str(spec)] = param_linelist_df['nu'].values
                order_nu.append('nu_' + str(spec))
                param_linelist_df['nu_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['nu_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                if vary_nu != {}:
                    for molecule in vary_nu:
                        for isotope in vary_nu[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) & (param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nu_' + str(spec) + '_vary'] = (vary_nu[molecule][isotope])

        # Defines the linestrength in the event that the sw is held constant across all samples and not
        ## Starting point is equal to the initial value
        order_sw = ['sw', 'sw_scale_factor']
        param_linelist_df['sw_vary'] = len(param_linelist_df)*[False]
        param_linelist_df['sw_err'] = len(param_linelist_df)*[0]
        if self.sw_constrain:
            if vary_sw != {}:
                 for molecule in vary_sw:
                    for isotope in vary_sw[molecule]:
                        param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'sw_vary'] = (vary_sw[molecule][isotope])
        else:
            for spec in self.dataset.get_list_spectrum_numbers():
                param_linelist_df['sw_' + str(spec)] = param_linelist_df['sw'].values
                order_sw.append('sw_' + str(spec))
                param_linelist_df['sw_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['sw_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                if vary_sw != {}:
                    for molecule in vary_sw:
                        for isotope in vary_sw[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'sw_' + str(spec) + '_vary'] = (vary_sw[molecule][isotope])

        #Loop through other parameters and then set things to 0 based on lineshape
        order_gamma0 = []
        order_delta0 = []
        order_SD_gamma = []
        order_SD_delta = []
        order_nuVC = []
        order_eta = []
        order_linemixing = []
        for diluent in diluent_list:
            #Gamma0 option for constrain and not constrained
            order_gamma0.append('gamma0_' + diluent)
            param_linelist_df['gamma0_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['gamma0_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.gamma0_constrain:
                if vary_gamma0 != {}:
                    for molecule in vary_gamma0:
                        for isotope in vary_gamma0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'gamma0_' +diluent + '_vary'] = (vary_gamma0[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_gamma0.append('gamma0_' +diluent + '_' +str(spec))
                    param_linelist_df['gamma0_' +diluent + '_' +str(spec)] = (param_linelist_df['gamma0_' + diluent].values)
                    param_linelist_df['gamma0_' + diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['gamma0_'+ diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if vary_gamma0 != {}:
                        for molecule in vary_gamma0:
                            for isotope in vary_gamma0[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'gamma0_' +diluent +'_'+str(spec) + '_vary'] = (vary_gamma0[molecule][isotope])
            order_gamma0.append('n_gamma0_' +diluent )

            #Delta0 option for constrain and not constrained
            order_delta0.append('delta0_' + diluent)
            param_linelist_df['delta0_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['delta0_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.delta0_constrain:
                if vary_delta0 != {}:
                    for molecule in vary_delta0:
                        for isotope in vary_delta0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'delta0_' +diluent + '_vary'] = (vary_delta0[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_delta0.append('delta0_' +diluent + '_' +str(spec))
                    param_linelist_df['delta0_' +diluent + '_' +str(spec)] = (param_linelist_df['delta0_' + diluent].values)
                    param_linelist_df['delta0_' +  diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[0]
                    param_linelist_df['delta0_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if vary_delta0 != {}:
                        for molecule in vary_delta0:
                            for isotope in vary_delta0[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'delta0_' +diluent +'_'+str(spec) + '_vary'] = (vary_delta0[molecule][isotope])
            order_delta0.append('n_delta0_' +diluent )

            #SD Gamma option for constrain and not constrained
            order_SD_gamma.append('SD_gamma_' + diluent )
            param_linelist_df['SD_gamma_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['SD_gamma_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.aw_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                    param_linelist_df.loc[:, 'SD_gamma_' + diluent] = 0
                else:
                    if vary_aw != {}:
                        for molecule in vary_aw:
                            for isotope in vary_aw[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_gamma_' +diluent + '_vary'] = (vary_aw[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_SD_gamma.append('SD_gamma_' +diluent + '_' +str(spec))
                    param_linelist_df['SD_gamma_' + diluent + '_'+ str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['SD_gamma_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                        param_linelist_df.loc[:, 'SD_gamma_' +diluent + '_' +str(spec)] = 0
                    else:
                        param_linelist_df['SD_gamma_' +diluent + '_' +str(spec)] = (param_linelist_df['SD_gamma_' + diluent].values)
                        if vary_aw != {}:
                            for molecule in vary_aw:
                                for isotope in vary_aw[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_gamma_' +diluent +'_'+str(spec) + '_vary'] = (vary_aw[molecule][isotope])
            order_SD_gamma.append('n_gamma2_' +diluent )

            #SD Delta option for constrain and not constrained
            order_SD_delta.append('SD_delta_' + diluent)
            param_linelist_df['SD_delta_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['SD_delta_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.as_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                    param_linelist_df.loc[:, 'SD_delta_' + diluent] = 0
                else:
                    if vary_as != {}:
                        for molecule in vary_as:
                            for isotope in vary_as[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_delta_' +diluent + '_vary'] = (vary_as[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_SD_delta.append('SD_delta_' +diluent + '_' +str(spec))
                    param_linelist_df['SD_delta_' + diluent + '_'+ str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['SD_delta_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                        param_linelist_df.loc[:, 'SD_delta_' +diluent + '_' +str(spec)] = 0
                    else:
                        param_linelist_df['SD_delta_' +diluent + '_' +str(spec)] = (param_linelist_df['SD_delta_' + diluent].values)

                        if vary_as != {}:
                            for molecule in vary_as:
                                for isotope in vary_as[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_delta_' +diluent +'_'+str(spec) + '_vary'] = (vary_as[molecule][isotope])
            order_SD_delta.append('n_delta2_' +diluent )

            #nuVC option for constrain and not constrained
            order_nuVC.append('nuVC_' + diluent)
            param_linelist_df['nuVC_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['nuVC_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.nuVC_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP'):
                    param_linelist_df.loc[:, 'nuVC_' + diluent ] = 0
                else:
                    if vary_nuVC != {}:
                        for molecule in vary_nuVC:
                            for isotope in vary_nuVC[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nuVC_' +diluent + '_vary'] = (vary_nuVC[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_nuVC.append('nuVC_' + diluent + '_'+str(spec))
                    param_linelist_df['nuVC_' + diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['nuVC_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP'):
                        param_linelist_df.loc[:, 'nuVC_' + diluent + '_'+str(spec)] = 0
                    else:
                        param_linelist_df['nuVC_' +diluent + '_' +str(spec)] = (param_linelist_df['nuVC_' + diluent].values)

                        if vary_nuVC != {}:
                            for molecule in vary_nuVC:
                                for isotope in vary_nuVC[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nuVC_' +diluent +'_'+str(spec) + '_vary'] = (vary_nuVC[molecule][isotope])
            order_nuVC.append('n_nuVC_' +diluent )

            #eta option for constrain and not constrained
            order_eta.append('eta_' + diluent)
            param_linelist_df['eta_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['eta_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.eta_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP') or (self.lineprofile == 'NGP') or (self.lineprofile == 'SDNGP'):
                    param_linelist_df.loc[:, 'eta_' + diluent] = 0
                else:
                    if vary_eta != {}:
                        for molecule in vary_eta:
                            for isotope in vary_eta[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'eta_' +diluent + '_vary'] = (vary_eta[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_eta.append('eta_' +diluent + '_' +str(spec))
                    param_linelist_df['eta_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['eta_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP') or (self.lineprofile == 'NGP') or (self.lineprofile == 'SDNGP'):
                        param_linelist_df.loc[:, 'eta_' +diluent + '_' +str(spec)] = 0
                    else:
                        param_linelist_df['eta_' +diluent + '_' +str(spec)] = (param_linelist_df['eta_' + diluent].values)
                        if vary_eta != {}:
                            for molecule in vary_eta:
                                for isotope in vary_eta[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'eta_' +diluent +'_'+str(spec) + '_vary'] = (vary_eta[molecule][isotope])

            # Linemixing
            order_linemixing.append('y_' + diluent)
            param_linelist_df['y_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['y_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.linemixing:
                if self.linemixing_constrain:
                    if vary_linemixing != {}:
                        for molecule in vary_linemixing:
                            for isotope in vary_linemixing[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'y_' +diluent + '_vary'] = (vary_linemixing[molecule][isotope])
                else:
                    for spec in self.dataset.get_list_spectrum_numbers():
                        order_linemixing.append('y_' +diluent + '_' +str(spec))
                        param_linelist_df['y_' +diluent + '_' +str(spec)] = (param_linelist_df['y_' + diluent].values)
                        param_linelist_df['y_' + diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[False]
                        param_linelist_df['y_'+ diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                        if vary_linemixing != {}:
                            for molecule in vary_linemixing:
                                for isotope in vary_linemixing[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'y_' +diluent +'_'+str(spec) + '_vary'] = (vary_linemixing[molecule][isotope])
            else:
                param_linelist_df['y_' + diluent] = 0
            order_linemixing.append('n_y_' +diluent )
            
            #Temperature Dependence
            if num_nominal_temps > 1:
                param_linelist_df['n_gamma0_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_gamma0_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_delta0_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_delta0_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_gamma2_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_gamma2_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_delta2_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_delta2_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_nuVC_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_nuVC_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_y_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_y_'+diluent+'_err'] = len(param_linelist_df)*[0]
                #n_Gamma0
                if vary_n_gamma0 != {}:
                    for molecule in vary_n_gamma0:
                        for isotope in vary_n_gamma0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_gamma0_' +diluent + '_vary'] = (vary_n_gamma0[molecule][isotope])
                #n_Delta0
                if vary_n_delta0 != {}:
                    for molecule in vary_n_delta0:
                        for isotope in vary_n_delta0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_delta0_' +diluent + '_vary'] = (vary_n_delta0[molecule][isotope])
                #n_Gamma2
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'NGP') :
                    if vary_n_gamma2 != {}:
                        for molecule in vary_n_gamma2:
                            for isotope in vary_n_gamma2[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_gamma2_' +diluent + '_vary'] = (vary_n_gamma2[molecule][isotope])
                #n_Delta2
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'NGP') :
                    if vary_n_delta2 != {}:
                        for molecule in vary_n_delta2:
                            for isotope in vary_n_delta2[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_delta2_' +diluent + '_vary'] = (vary_n_delta2[molecule][isotope])
                #n_nuVC
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'SDVP') :
                    if vary_n_nuVC != {}:
                        for molecule in vary_n_nuVC:
                            for isotope in vary_n_nuVC[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_nuVC_' +diluent + '_vary'] = (vary_n_nuVC[molecule][isotope])
                #n_y
                if not (self.linemixing) :
                    if vary_n_linemixing != {}:
                        for molecule in vary_n_linemixing:
                            for isotope in vary_n_linemixing[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_y_' +diluent + '_vary'] = (vary_n_linemixing[molecule][isotope])
        
        ordered_list = self.additional_columns.copy()
        ordered_list += ['molec_id', 'local_iso_id','elower']

        for item in order_nu:
            ordered_list.append(item)
            ordered_list.append(item + '_err')
            ordered_list.append(item + '_vary')
        for item in order_sw:
            ordered_list.append(item)
            if item != 'sw_scale_factor':
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
        for item in order_gamma0:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_delta0:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_nuVC:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_SD_gamma:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_SD_delta:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_eta:
            ordered_list.append(item)
            ordered_list.append(item + '_err')
            ordered_list.append(item + '_vary')
        for item in order_linemixing:          
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        param_linelist_df = param_linelist_df[ordered_list]
        param_linelist_df.to_csv(self.param_linelist_savename + '.csv') #

        #Warnings of SD_gamma, SD_delta, nuVC, eta are zero and floated
        for param in ordered_list:
            if ('SD_gamma' in param) & ('vary' not in param) & ('err' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating SD_gamma terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.')
            if ('SD_delta' in param) & ('vary' not in param) & ('err' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating SD_delta terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.')
            if ('nuVC' in param) & ('vary' not in param) & ('err' not in param) & ('n_nuVC' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating nuVC terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.  This is generally less of an issue than the case where SD_gamma is floated when the initial guess is zero.')

        return param_linelist_df

    def generate_fit_baseline_linelist(self, vary_baseline = True, vary_pressure = False, vary_temperature = False,vary_molefraction = {7:True, 1:False}, vary_xshift = False,
                                      vary_etalon_amp= False, vary_etalon_period= False, vary_etalon_phase= False,
                                      vary_ILS_res = False):
        """Generates the baseline line list used in fitting and updates the fitting booleans to desired settings.


        Parameters
        ----------
        vary_baseline : bool, optional
            If True, then sets all baseline parameters for all spectra to float. The default is True.
        vary_pressure : bool, optional
            If True, then the pressures for all spectra are floated.  This should be used with caution as the impact these parameters have on other floated parameters might lead to an unstable solution. The default is False.
        vary_temperature : bool, optional
            If True, then the temperatures for all spectra are floated. This should be used with caution as the impact these parameters have on other floated parameters might lead to an unstable solution.The default is False.
        vary_molefraction : dict, optional
            keys to dictionary correspond to molecule id where the value is boolean flag, which dictates whether to float the mole fraction. The default is {}. Example: {7: True, 1: False}
        vary_xshift : bool, optional
            If True, then sets x-shift parameters for all spectra to float. The default is False.
        vary_etalon_amp : bool, optional
            If True, then sets etalon amplitude parameters for all spectra to float. The default is False.
        vary_etalon_period : bool, optional
            If True, then sets etalon period parameters for all spectra to float. . The default is False.
        vary_etalon_phase : bool, optional
            If True, then sets etalon phase parameters for all spectra to float.. The default is False.
        vary_ILS_res : bool, optional
            If True, then sets ILS resolution parameters for all spectra to float.. The default is False.

        Returns
        -------
        base_linelist_df : dataframe
            returns dataframe based on baseline line list with addition of a vary and err column for every floatable parameter.  The vary columns are defined by the inputs.  The err columns will be populated from fit results.  The dataframe is also saved as a .csv file..



        """

        base_linelist_df = self.get_base_linelist().copy()
        parameters =  (list(base_linelist_df))
        baseline_param_order = ['Segment Number']

        #Generate Fit Baseline file
        for param in parameters:
            if ('Baseline Order' != param) and ('Segment Number' != param):
                base_linelist_df[param + '_err'] = 0
                base_linelist_df[param + '_vary']= False
                baseline_param_order += [param, param + '_err', param + '_vary']

            if 'Pressure' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_pressure)]
                if (vary_pressure):
                    print ('USE CAUTION WHEN FLOATING PRESSURES')
            if 'Temperature' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_temperature)]
                if (vary_temperature):
                    print ('USE CAUTION WHEN FLOATING TEMPERATURES')
            if 'x_shift' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_xshift)]
            if 'baseline' in param:
                order = ord(param.replace('baseline_', '')) - 97
                base_linelist_df.loc[base_linelist_df['Baseline Order']>= order, param + '_vary'] = vary_baseline
            if 'molefraction' in param:
                for molecule in vary_molefraction:
                    if (self.dataset.isotope_list[(molecule, 1)][4]) in param:
                        base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_molefraction[molecule])
            if 'amp' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_etalon_amp)
            if 'period' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_etalon_period)
            if 'phase' in param:
                base_linelist_df.loc[base_linelist_df[param.replace("phase", "period")]!=0, param + '_vary'] = (vary_etalon_phase)
            if '_res_' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_ILS_res)


        #base_linelist_df.drop(['Baseline Order'], axis=1, inplace = True)
        base_linelist_df = base_linelist_df[baseline_param_order]
        #base_linelist_df = base_linelist_df.reindex(sorted(base_linelist_df.columns), axis=1)
        base_linelist_df.to_csv(self.base_linelist_savename + '.csv', index = True)
        return base_linelist_df
