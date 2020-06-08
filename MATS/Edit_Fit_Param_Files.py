#Import Packages
from Utilities import *



class Edit_Fit_Param_Files:
    """Allows for the baseline and parameter files to be editted in jupyter notebook.  Can also edit everything in .csv
    
    Parameters
    ----------
    base_linelist_file : str
        name of the .csv file containing the baseline parameters generated in format specified in the generate fit parameters class.
    param_linelist_file : str
        name of the .csv file containing the linelist parameters generated in format specified in the generate fit parameters class.
    new_base_linelist_file : str
        name of the to save the baseline param list as after editing. IF the value is None (default) then it will edit the input.
    new_param_linelist_file  : str
        name of the to save the linelist param list as after editing. IF the value is None (default) then it will edit the input.
    """

    def __init__(self, base_linelist_file, param_linelist_file, new_base_linelist_file = None, new_param_linelist_file = None):
        self.base_linelist_file = base_linelist_file
        self.param_linelist_file = param_linelist_file
        if new_base_linelist_file == None:
            self.new_base_linelist_file = base_linelist_file
        else:
            self.new_base_linelist_file = new_base_linelist_file
        if new_param_linelist_file == None:
            self.new_param_linelist_file = param_linelist_file
        else:
            self.new_param_linelist_file = new_param_linelist_file
    def edit_generated_baselist(self):
        """Allows editting of baseline linelist in notebook        
        """

        base_linelist_df = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0)
        baseline_widget = qgrid.show_grid(base_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return baseline_widget
    def save_edited_baselist(self, baseline_widget):
        """Saves edits of baseline linelist in notebook        
        """

        base_linelist_df = baseline_widget.get_changed_df()
        base_linelist_df.to_csv(self.new_base_linelist_file + '.csv', index = False)
        return base_linelist_df
    def edit_generated_paramlist(self):
        """Allows editting of parameter linelist in notebook        
        """

        param_linelist_df = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        param_widget = qgrid.show_grid(param_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return param_widget
    def save_edited_paramlist(self, param_widget):
        """Saves edits of parameter linelist in notebook        
        """

        param_linelist_df = param_widget.get_changed_df()
        param_linelist_df.to_csv(self.new_param_linelist_file + '.csv') 
        return param_linelist_df