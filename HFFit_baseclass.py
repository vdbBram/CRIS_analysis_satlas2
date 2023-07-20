import os
import numpy as np
import pandas as pd
import satlas2
from typing import Union
from numpy.typing import ArrayLike

def calculate_variation_dict(I: float, Jlower: float, Jupper: float) -> dict:
    '''Create a dictionary with booleans indicating which ABC values should be varied during fitting.
    MADE FOR DIPOLE AND QUADRUPOLE, REMOVES OCTUPOLE HF PARAMETERS

    Parameters
    ----------
    I: float
        the spin
    Jlower: float
        the total angular momentum of the lower electronic level
    Jupper: float
        the total angular momentum of the lower electronic level
    Returns
    -------
    dict
    '''
    if I == 0:
        return {'Al': False, 'Au': False, 'Bu': False, 'Bl': False, 'Cu': False, 'Cl': False}
    elif (I == 0.5):
        return {'Bu': False, 'Bl': False, 'Cu': False, 'Cl': False}
    elif (Jlower == 0.5 and Jupper == 0.5):
        return {'Bu': False, 'Bl': False, 'Cu': False, 'Cl': False}
    elif (Jlower == 0.5):
        return {'Bl': False, 'Cu': False, 'Cl': False}
    elif (Jupper == 0.5):
        return {'Bu': False, 'Cu': False, 'Cl': False}
    else:
        return {'Cu': False, 'Cl': False}

# just to make labels that are 7/2 instead of 3.5
def convert_decimalstring_fractionstring(stringfloat: str) -> str:
    '''Converts a decimal string to the fraction in string format

    Parameters
    ----------
    stringfloat: str
        decimal in string format
    Returns
    -------
    str
    '''
    if type(stringfloat) == type(list()):
        stringfraction = list()
        for el in stringfloat:
            if el.as_integer_ratio()[1] == 1:
                stringfraction.append(str(el))

            else:
                temptop = str(el.as_integer_ratio()[0])
                tempbottom = str(el.as_integer_ratio()[1])
                stringfraction.append(r'$\frac{' + temptop + '}{' + tempbottom + '}$')
    else:
        if stringfloat.as_integer_ratio()[1] == 1:
            stringfraction = str(stringfloat)
        else:
            temptop = str(stringfloat.as_integer_ratio()[0])
            tempbottom = str(stringfloat.as_integer_ratio()[1])
            stringfraction = r'$\frac{' + temptop + '}{' + tempbottom + '}$'
    return stringfraction

def create_spin_label(inputstring: str) -> str:
    '''Creates a legend label in fraction format for a spin

    Parameters
    ----------
    inputstring: str
        spin as decimal in string format
    Returns
    -------
    str
    '''
    pre__,post__ = inputstring.split('__')[0],inputstring.split('__')[1]
    if len(post__.split('_')) > 1:
        val,dec_val = post__.split('_')[0],post__.split('_')[1]
        spin = float(val+'.'+dec_val)
    else:
        spin = int(post__)
    return 'I = ' + convert_decimalstring_fractionstring(spin)

# input_param should be dict of dict eg. for data with 2 isomers and background, first element in dict is first isomer and second element is second isomer with all the respective nuclear info and last element is the background
class HF_fitter:
    # all possible arguments for the input_param and plot_args
    _all_input_param = ['I', 'J', 'ABC', 'centroid', 'scale', 'background_values', 'background_bounds','racah','variation_HFparam', 'set_values', 'set_boundaries', 'set_priors', 
    'expressions', 'variation_param', 'share_params', 'fwhmg', 'fwhml']
    _all_plot_args = ['legend_fontsize', 'ax_fontsize', 'ticker_locator', 'ticker_fontsize', 'nrows', 'ncols', 'data_ecolor', 'data_fillstyle', 'data_markersize', 'data_fmt', 
    'range_x', 'range_y']
    _colors = ['orange', 'green', 'darkorchid']
    # save_path = {save_folder_param, save_folder_fig}
    def __init__(self, MASS: int, SCANS: list, x: ArrayLike, y: ArrayLike, xerr: ArrayLike, yerr: ArrayLike, bunches: Union[ArrayLike, int], input_param: dict, 
        save_paths: list, **plot_args) -> None:
        '''Initialises the class to fit HF spectra

        Parameters
        ----------
        MASS: int
            mass number of the isotope
        SCANS: list
            scan numbers of the data
        x: ArrayLike
            x values of data in relative frequency format in MHz
        y: ArrayLike
            total counts on that x value
        xerr: ArrayLike
            xerr values of data in MHz
        yerr: ArrayLike
            error on y 
        bunches: ArrayLike
            the amount of bunches on each x. If using total counts for fitting, supply np.ones(len(x))
        input_param: dict
            dictionary containing all fitting information
        save_paths: list 
            list of length 2 with the savepath for the HF parameters and the figure of the fit
        '''
        self.hf_models = []
        self._mass = MASS 
        self._scans = SCANS
        self._x = x
        self._y = y/bunches
        self._init_y = y 
        self._xerr = xerr
        self._yerr = yerr/bunches
        self._init_yerr = yerr
        self._bunches = bunches
        self._input_param = input_param
        self._plot_args = plot_args
        self.set_input_param_values()
        self.set_plot_args_values()
        self._savepaths_param = f'{save_paths[0]}{self._mass}_'
        self._savepath_fig = f'{save_paths[1]}{self._mass}'

    def make_save_paths(self) -> None:
        '''Checks if the savepath is valid and create the folder if neccesary'''
        save_path_param = {}
        save_path_fig = self._savepath_fig
        for model in self.hf_models:
            save_path_param[model.name] = self._savepaths_param + str(self._input_param[model.name]['I'])
            save_path_fig += '_'+str(self._input_param[model.name]['I'])
            if not os.path.isdir(save_path_param[model.name]):
                os.mkdir(save_path_param[model.name])
        self._savepaths_param = save_path_param
        self._savepath_fig = save_path_fig
        if not os.path.isdir(self._savepath_fig):
            os.mkdir(self._savepath_fig)
    
    def set_input_param_values(self) -> None:
        '''fills in the missing keys to the input_param dictionary to default values'''
        default_kwargs = {'background_bounds':False, 'background_values':False, 'variation_HFparam':False, 'expressions':False, 'variation_param':False, 
        'racah':False, 'set_values':False, 'set_boundaries':False, 'set_priors':False, 'share_params':False, 'fwhmg':50, 'fwhml':50}
        for key in self._input_param.keys():
            self._input_param[key] = {**default_kwargs, **self._input_param[key]}

    def set_plot_args_values(self) -> None:
        '''Fills in the missing keys to the plot_args dictionary to default values'''
        default_kwargs = {'legend_fontsize':20, 'ax_fontsize':25, 'ticker_locator':1000, 'ticker_fontsize':20, 'nrows':2, 'ncols':1, 
        'data_ecolor':'k', 'data_fillstyle':'none', 'data_markersize':15, 'data_fmt':'.'}
        self._plot_args = {**default_kwargs, **self._plot_args}

    def init_all_models(self, input_param: dict) -> None:
        '''initialises all satlas2.Model objects with the model name == the keys of input_param, and the parameters following the dictionary in input_param[model_name] 

        Parameters
        ----------
        input_param: dict
            Dictionary of dictionaries where the top dictionary represent the different models with as keys the model names, and the bottom dictionaries the different 
            parameters that are associated with each model. eg {model_name1: {'I':3, 'J':[0.5,1.5], etc.}, model_name2:{'background_values':[15,10], 'background_bounds':[0]}, etc.} 
        '''
        class CustomLlhFitter(satlas2.Fitter):
            def customLlh(self) -> ArrayLike:
                try:
                    bunches = self.bunches
                except:
                    bunches = self.getSourceAttr('bunches')
                    self.bunches = bunches
                try:
                    data_counts = self.data_counts
                except:
                    data_counts = self.temp_y * bunches
                    self.data_counts = data_counts
                model_rates = self.f()
                model_counts = model_rates * bunches
                returnvalue = data_counts * np.log(model_counts) - model_counts
                returnvalue[model_counts <= 0] = -np.inf
                priors = self.gaussianPriorResid()
                if len(priors) > 1:
                    priors = -0.5 * priors * priors
                    returnvalue = np.append(returnvalue, priors)
                return returnvalue
        self.fitter = CustomLlhFitter()
        self.models = []
        self.datasource = satlas2.Source(x = self._x, y = self._y, xerr = self._xerr, yerr = self._yerr, name = 'source', bunches = self._bunches)
        self.fitter.addSource(self.datasource)
        for i,key in enumerate(input_param.keys()):
            if input_param[key]['share_params'] != False:
                self.fitter.shareParams(input_param[key]['share_params'])
            self.models.append(self.init_model(input_param[key],key))
            self.datasource.addModel(self.models[i])

    def init_model(self, input_param: dict, model_name: str) -> satlas2.Model: 
        '''Initialises a satlas2.Model with the name == model_name and parameters given in the input_param dictionary

        Parameters
        ----------
        input_param: dict
            Dictionary with the different parameters that are associated with this model. eg {'I':3, 'J':[0.5,1.5], etc.}
        model_name: str
            name of the model that is initialised
        Returns
        -------
        satlas2.Model
        '''
        if input_param['background_bounds'] == False and input_param['background_values'] == False:
            model = satlas2.HFS(I = input_param['I'], 
                          J = input_param['J'],   
                          A = input_param['ABC'][0:2],
                          B = input_param['ABC'][2:4],
                          C = input_param['ABC'][4:6],
                          df = input_param['centroid'],
                          scale = input_param['scale'],
                          racah = input_param['racah'],
                          fwhmg = input_param['fwhmg'],
                          fwhml = input_param['fwhml'],
                          name = model_name)
            self.hf_models.append(model)
            model = self.add_hfs_constraints(model, input_param)
            return model
        elif input_param['background_bounds'] == False:
            model = satlas2.Polynomial(input_param['background_values'], name = model_name)
            model = self.add_constraints(model, input_param)
            return model
        else:
            model = satlas2.PiecewiseConstant(input_param['background_values'], input_param['background_bounds'], name = model_name)
            model = self.add_constraints(model, input_param)
            return model

    def add_hfs_constraints(self, model: satlas2.Model, input_param: dict) -> satlas2.Model:
        '''adds variation constraints to the HF parameters of the supplied model

        Parameters
        ----------
        model: satlas2.Model
            the model to apply the constraints to
        input_param: dict
            Dictionary with the different parameters that are associated with this model. eg {'I':3, 'J':[0.5,1.5], etc.}
        Returns
        -------
        satlas2.Model
        '''
        model = self.add_constraints(model, input_param)
        if input_param['variation_HFparam'] == False: # if you do not want to constrain any HF params manually, constrains Cl and Cu automatically
            for key in calculate_variation_dict(input_param['I'], input_param['J'][0], input_param['J'][1]):
                model.params[key].vary = False
        else:
            for key in input_param['variation_HFparam'].keys():
                model.params[key].vary = False
        return model

    def add_constraints(self, model: satlas2.Model, input_param: dict) -> satlas2.Model:
        '''adds constraints to the parameters of the supplied model

        Parameters
        ----------
        model: satlas2.Model
            the model to apply the constraints to
        input_param: dict
            Dictionary with the different parameters that are associated with this model. eg {'I':3, 'J':[0.5,1.5], etc.}
        Returns
        -------
        satlas2.Model
        '''
        if input_param['set_values'] != False:
            for param in input_param['set_values']:
                model.params[param].vary = False
                model.params[param].value = input_param['set_values'][param]
        if input_param['variation_param'] != False: # if you do not want to constrain any params manually
            for key in input_param['variation_param'].keys():
                model.params[key].vary = False
        if input_param['set_boundaries'] != False: # to set boundaries on your fit params
            for param in input_param['set_boundaries'].keys():
                model.params[param].min = input_param['set_boundaries'][param].get('min', -np.inf)
                model.params[param].max = input_param['set_boundaries'][param].get('max', np.inf)
        if input_param['set_priors'] != False:
            for key in input_param['set_priors'].keys():
                value, sigma = input_param['set_priors'][key][0], input_param['set_priors'][key][1]
                self.fitter.setParamPrior(self.fitter.sources[0].name, model.name, key, value, sigma) 
        if input_param['expressions'] != False: # sets contraints on the intensities if given
            for key in input_param['expressions'].keys():
                parameter_name_to_fix = f'{self.fitter.sources[0][1].name}___{model.name}___{key}'
                factor = input_param['expressions'][key][0]
                parameter_name_depend = self.fitter.sources[0][1].name+'___'+model.name+'___'+input_param['expressions'][key][1]
                self.fitter.setExpr(parameter_name_to_fix, f'{factor}*{parameter_name_depend}')
        return model

    def fit_model(self, show_correl: bool = False, min_correl: float = 0.1, llh: bool = False, **llh_kwargs) -> satlas2.Fitter:
        '''fits the models supplied to the fitter object of this class and prints the report of the fit

        Parameters
        ----------
        show_correl: bool
            True to show the correlation coefficients
        min_correl: float
            minimum correlation for it to be shown
        llh: bool
            True if using MLE fitting
        Returns
        -------
        satlas2.Fitter
        '''
        self.fitter.fit(llh = llh, **llh_kwargs)
        print(self.fitter.reportFit(show_correl = show_correl))
        return self.fitter

    def save_results(self) -> pd.DataFrame:
        '''Saves fitted results to .csv file format and returns the pd.DataFrame
        Returns
        -------
        pd.DataFrame
        '''
        self.result_DataFrame = self.fitter.createResultDataframe()
        self.make_save_paths()
        # add redchi
        if len(self._scans) == 2:
            try: # C:\\Users...\\save_folder\\ + mass_I\\scan-4708.csv
                self.result_DataFrame.to_csv(self._savepaths_param[self.hf_models[0].name] + '\\scan-' + str(self._scans[0]) + '-' + str(self._scans[1]) + '.csv', sep = ';')
            except:
                os.mkdir(self._savepaths_param[self.hf_models[0].name])
                self.result_DataFrame.to_csv(self._savepaths_param[self.hf_models[0].name] + '\\scan-' + str(self._scans[0]) + '-' + str(self._scans[1]) + '.csv', sep = ';')
        else:
            try: # C:\\Users...\\save_folder\\ + mass_I\\scan-4708.csv
                self.result_DataFrame.to_csv(self._savepaths_param[self.hf_models[0].name] + '\\scan-' + str(self._scans[0]) + '.csv', sep = ';')
            except:
                os.mkdir(self._savepaths_param[self.hf_models[0].name])
                self.result_DataFrame.to_csv(self._savepaths_param[self.hf_models[0].name] + '\\scan-' + str(self._scans[0]) + '.csv', sep = ';')
        return self.result_DataFrame
