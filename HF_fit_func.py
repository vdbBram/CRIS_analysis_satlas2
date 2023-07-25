from interpolate_refscans_centroid import *
from HFFit_baseclass import *
from typing import Tuple

def check_ref(initial_param: dict, MASS: int, I_ref: float, mass_ref: int) -> bool:
    '''checks if the scan is the reference isomer scan

        Parameters
        ----------
        initial_param: dict
            Dictionary of dictionaries where the top dictionary represent the different models with as keys the model names, and the bottom dictionaries the different 
            parameters that are associated with each model. eg {model_name1: {'I':3, 'J':[0.5,1.5], etc.}, model_name2:{'background_values':[15,10], 'background_bounds':[0]}, etc.} 
        MASS: int
            mass number of the scan
        I_ref: float
            spin of the reference isomer
        mass_ref: int
            mass number of the reference isomer
        Returns
        -------
        bool
        '''
    for model_name in initial_param.keys():
        if model_name[0].lower() == 's':
            if  initial_param[model_name]['I'] == I_ref and MASS == mass_ref:
                return True
    return False

def param_test_plot(nb_datasets: int, fig: plt.figure, ax: Union[plt.Axes, ArrayLike], data: ArrayLike, initial_param: dict, range_x: list, **plot_kwargs) -> Union[plt.Axes,ArrayLike]:
    '''plots initial parameter guess together with the data

        Parameters
        ----------
        nb_datasets: int
            number of scans for this nuclear state
        fig: plt.figure
            figure to plot the data and estimate on
        ax: plt.Axes or ArrayLike
            the ax or axes to plot the estimates and data on
        data: ArrayLike
            the experimental data
        initial_param: dict
            Dictionary of dictionaries where the top dictionary represent the different models with as keys the model names, and the bottom dictionaries the different 
            parameters that are associated with each model. eg {model_name1: {'I':3, 'J':[0.5,1.5], etc.}, model_name2:{'background_values':[15,10], 'background_bounds':[0]}, etc.} 
        range_x: list
            minimum and maximum x values to show in the plot
        Returns
        -------
        plt.Axes or ArrayLike
        '''
    colors = (plt.get_cmap()(np.linspace(0,1,len(list(initial_param.keys())))))
    bunches = data[4]
    testsource = satlas2.Source(x = data[0], y = data[2]/bunches, xerr = data[1], yerr = data[3]/bunches, name = 'Test')
    if nb_datasets == 1:
        ax.errorbar(x = data[0], y = data[2]/bunches, xerr = data[1], yerr = data[3]/bunches, **plot_kwargs)
        ax.set_xlim(range_x[0],range_x[1])
        ax.set_xlabel('Observed frequency - transition frequency [MHz]')
        for i,model_name in enumerate(initial_param.keys()):
            if list(initial_param.keys())[i][0].lower() != 's':
                continue
            mod_init_param = initial_param[list(initial_param.keys())[i]]
            testmodel = satlas2.HFS(I = mod_init_param['I'], J = mod_init_param['J'], A = mod_init_param['ABC'][0:2], B = mod_init_param['ABC'][2:4], C = mod_init_param['ABC'][4:], 
            df = mod_init_param['centroid'], scale = mod_init_param['scale'], racah = False, name = 'HFSmodel_'+str(mod_init_param['I']))
            testsource.addModel(testmodel)
            freq_range = np.arange(range_x[0],range_x[1],1)
            response_testmodel =  testmodel.f(freq_range)# + tmodel.f(freq_range)
            ax.plot(freq_range, response_testmodel, '-', color = colors[i], label = 'Initial guess model I = ' + convert_decimalstring_fractionstring(mod_init_param['I']))
        ax.legend(fontsize = 20)
        tmodel = satlas2.Polynomial(initial_param[model_name]['background_values'], name = 'bkg')
        testsource.addModel(tmodel)
        ax.plot(freq_range,tmodel.f(freq_range), label = 'bkg')
        return ax
    for dataset_nb in range(nb_datasets):
        ax[dataset_nb].errorbar(x = data[0], y = data[2]/bunches, xerr = data[1], yerr = data[3]/bunches, **plot_kwargs)
        ax[dataset_nb].set_xlim(range_x[dataset_nb][0],range_x[dataset_nb][1])
        ax[dataset_nb].set_xlabel('Observed frequency - transition frequency [MHz]')
        for i,model_name in enumerate(initial_param.keys()):
            if list(initial_param.keys())[i][0].lower() != 's':
                continue
            mod_init_param = initial_param[list(initial_param.keys())[i]]
            testmodel = satlas2.HFS(I = mod_init_param['I'], J = mod_init_param['J'], A = mod_init_param['ABC'][0:2], B = mod_init_param['ABC'][2:4], C = mod_init_param['ABC'][4:], 
            df = mod_init_param['centroid'], scale = mod_init_param['scale'], racah = False, name = 'HFSmodel_'+str(mod_init_param['I']))
            testsource.addModel(testmodel)
            freq_range = np.arange(range_x[dataset_nb][0],range_x[dataset_nb][1],1)
            response_testmodel =  testmodel.f(freq_range)
            ax[dataset_nb].plot(freq_range, response_testmodel, '-', color = colors[i], label = 'Initial guess model I = ' + convert_decimalstring_fractionstring(mod_init_param['I']))
            ax[dataset_nb].legend()
    return ax

def set_plot_labels(ax: Union[plt.Axes, ArrayLike], fitter: HF_fitter, nb_datasets: int) -> Union[plt.Axes,ArrayLike]:
    '''sets the labels, fontsizes, etc. for the plots of the fit

        Parameters
        ----------
        ax: plt.Axes or ArrayLike
            the ax or axes to set the labels, etc. for
        fitter: HF_fitter obj
            Fitter object from HFFit_baseclass
        nb_datasets: int
            number of scans for this nuclear state
        Returns
        -------
        plt.Axes or ArrayLike
        '''
    if nb_datasets == 1:
        ax[0].set_ylabel('Events per bunch', fontsize = fitter._plot_args['ax_fontsize'])
        ax[0].set_xlim(fitter._plot_args['range_x'][0],fitter._plot_args['range_x'][1])
        ax[0].set_ylim(fitter._plot_args['range_y'][0],fitter._plot_args['range_y'][1])
        ax[0].xaxis.set_major_locator(ticker.MultipleLocator(fitter._plot_args['ticker_locator']))
        ax[0].tick_params(axis='both', which='major', labelsize=fitter._plot_args['ticker_fontsize'])
        ax[1].set_ylabel('Residuals', fontsize = fitter._plot_args['ax_fontsize'])
        ax[1].set_xlabel('Relative frequency (MHz)', fontsize = fitter._plot_args['ax_fontsize'])
        ax[1].xaxis.set_major_locator(ticker.MultipleLocator(fitter._plot_args['ticker_locator']))
        ax[1].tick_params(axis='both', which='major', labelsize=fitter._plot_args['ticker_fontsize'])
    else:
        for dataset_nb in range(nb_datasets):
            ax[0][dataset_nb].set_ylabel('Events per bunch', fontsize = fitter._plot_args['ax_fontsize'])
            ax[0][dataset_nb].set_xlim(fitter._plot_args['range_x'][dataset_nb][0],fitter._plot_args['range_x'][dataset_nb][1])
            ax[0][dataset_nb].set_ylim(fitter._plot_args['range_y'][0],fitter._plot_args['range_y'][1])
            ax[0][dataset_nb].xaxis.set_major_locator(ticker.MultipleLocator(fitter._plot_args['ticker_locator']))
            ax[0][dataset_nb].tick_params(axis='both', which='major', labelsize=fitter._plot_args['ticker_fontsize'])
            ax[1][dataset_nb].set_ylabel('Residuals', fontsize = fitter._plot_args['ax_fontsize'])
            ax[1][dataset_nb].set_xlabel('Relative frequency (MHz)', fontsize = fitter._plot_args['ax_fontsize'])
            ax[1][dataset_nb].xaxis.set_major_locator(ticker.MultipleLocator(fitter._plot_args['ticker_locator']))
            ax[1][dataset_nb].tick_params(axis='both', which='major', labelsize=fitter._plot_args['ticker_fontsize'])
    return ax

def param_plot(nb_datasets: int, fig: plt.figure, ax: Union[plt.Axes,ArrayLike], data: ArrayLike, HF_fitter: HF_fitter) -> Union[plt.Axes,ArrayLike]:
    '''plots fitted parameter together with the data

        Parameters
        ----------
        nb_datasets: int
            number of scans for this nuclear state
        fig: plt.figure
            figure to plot the data and estimate on
        ax: plt.Axes or ArrayLike
            the ax or axes to plot the estimates and data on
        data: ArrayLike
            the experimental data
        fitter: HF_fitter obj
            Fitter object from HFFit_baseclass
        Returns
        -------
        plt.Axes or ArrayLike
        '''
    colors = (plt.get_cmap()(np.linspace(0,1,len(list(HF_fitter._input_param.keys()))+1)))
    bunches = data[4]
    res = (HF_fitter._y - HF_fitter.datasource.evaluate(HF_fitter._x)) / np.sqrt((HF_fitter._yerr*HF_fitter._yerr)) # residuals
    if nb_datasets == 1:
        ax[0].errorbar(x = data[0], y = data[2]/bunches, xerr = data[1], yerr = data[3]/bunches, fmt = HF_fitter._plot_args['data_fmt'], fillstyle = HF_fitter._plot_args['data_fillstyle'], 
            markersize = HF_fitter._plot_args['data_markersize'], ecolor = HF_fitter._plot_args['data_ecolor'])
        ax[1].plot(HF_fitter._x, res, '.')
        for i,model_name in enumerate(HF_fitter._input_param.keys()):
            if list(HF_fitter._input_param.keys())[i][0].lower() != 's':
                continue
            if len(HF_fitter.hf_models) == 1:
                freq_range = np.arange(HF_fitter._plot_args['range_x'][0],HF_fitter._plot_args['range_x'][1],1)
                ax[0].plot(freq_range, HF_fitter.datasource.evaluate(freq_range), '-', color = HF_fitter._colors[0], label = create_spin_label(HF_fitter.hf_models[0].name))
                break
            else:
                freq_range = np.arange(HF_fitter._plot_args['range_x'][0],HF_fitter._plot_args['range_x'][1],1)
                ax[0].plot(freq_range, HF_fitter.hf_models[i].f(freq_range), '-', color = HF_fitter._colors[i+1])
                ax[0].plot(freq_range, HF_fitter.datasource.evaluate(freq_range), '-', color = HF_fitter._colors[0], label = create_spin_label(HF_fitter.hf_models[i].name))
        ax = set_plot_labels(ax, HF_fitter, nb_datasets)
        ax[0].legend(fontsize = HF_fitter._plot_args['legend_fontsize'])
        return ax
    for dataset_nb in range(nb_datasets):
        ax[0][dataset_nb].errorbar(x = data[0], y = data[2]/bunches, xerr = data[1], yerr = data[3]/bunches, fmt = HF_fitter._plot_args['data_fmt'], fillstyle = HF_fitter._plot_args['data_fillstyle'], 
            markersize = HF_fitter._plot_args['data_markersize'], ecolor = HF_fitter._plot_args['data_ecolor'])
        ax[1][dataset_nb].plot(HF_fitter._x, res, '.')
        for i,model_name in enumerate(HF_fitter._input_param.keys()):
            if list(HF_fitter._input_param.keys())[i][0].lower() != 's':
                continue
            if len(HF_fitter.hf_models) == 1:
                freq_range = np.arange(HF_fitter._plot_args['range_x'][dataset_nb][0],HF_fitter._plot_args['range_x'][dataset_nb][1],1)
                ax[0][dataset_nb].plot(freq_range, HF_fitter.datasource.evaluate(freq_range), '-', color = HF_fitter._colors[0], label = create_spin_label(HF_fitter.hf_models[0].name))
                break
            else:
                freq_range = np.arange(HF_fitter._plot_args['range_x'][dataset_nb][0],HF_fitter._plot_args['range_x'][dataset_nb][1],1)
                ax[0][dataset_nb].plot(freq_range, HF_fitter.hf_models[i].f(freq_range), '-', color = colors[i], label = create_spin_label(HF_fitter.hf_models[i].name))
                ax[0][dataset_nb].plot(freq_range, HF_fitter.datasource.evaluate(freq_range), '-', color = HF_fitter._colors[0])
        ax = set_plot_labels(ax, HF_fitter, nb_datasets)
        ax[0][dataset_nb].legend(fontsize = HF_fitter._plot_args['legend_fontsize'])
    return ax

def fitter(test: bool, save: bool, llh: bool, interp_method: str, data: list, initial_param: dict, plot_param: dict, MASS: int, mass_ref: int, I_ref: float, SCANS: list, 
    DATA_FOLDER: str, SAVEPATH_FIG: str, SAVEPATH_PARAM: str, **llh_kwargs) -> Tuple[plt.figure,plt.Axes,satlas2.Fitter]:
    '''plots fitted parameter together with the data

        Parameters
        ----------
        test: bool
            True if test initial estimates of parameter, False to fit
        save: bool
            True if save fit results and fit, False if not
        llh: bool
            True if MLE fitting or random walk, False if least square
        interp_method: str
            interpolation method for reference centroid, choose spline or gp
        data: Arraylike
            experimental data in list. If multiple scans, then each element is the data of 1 scan. If 1 scan then data in a list with 1 element
        initial_param: dict
            Dictionary of dictionaries where the top dictionary represent the different models with as keys the model names, and the bottom dictionaries the different 
            parameters that are associated with each model. eg {model_name1: {'I':3, 'J':[0.5,1.5], etc.}, model_name2:{'background_values':[15,10], 'background_bounds':[0]}, etc.} 
        plot_param: dict
            dictionary with all plot parameters and their values, eg. {'x_range':[-1000,2000], 'y_range': [0,50], etc.}
        MASS: int
            mass number for this scan
        mass_ref: int
            mass number of the reference isomer
        I_ref: float
            spin of the reference isomer
        SCANS: list
            list of scan(s) to fit
        DATA_FOLDER: str
            path to dopplershifted data
        SAVEPATH_FIG: str
            path to where the figures should be saved
        SAVEPATH_PARAM: str
            path to where the fitted HF parameters should be saved
        Returns
        -------
        tuple of plt.figure, list of plt.Axes, satlas2.Fitter
        '''
    nb_datasets = len(data)
    bool_plot = False
    if interp_method.lower() == 'gp':
        bool_plot = True
    if not check_ref(initial_param, MASS, I_ref, mass_ref):
        for nb in range(nb_datasets):
            data[nb][:,0] = data[nb][:,0] - interpolate_refscans_centroid(DATA_FOLDER + 'Final_analysis_output\\Fitted_HFS_param\\' + str(mass_ref) + '_' + str(I_ref) + '\\', method = interp_method, plot = bool_plot)(int(SCANS[nb]))
    data = np.concatenate(data).T
    if test:
        fig,ax = plt.subplots(nrows = 1, ncols = nb_datasets, figsize = (14,9))
        ax = param_test_plot(nb_datasets = nb_datasets, fig = fig, ax = ax, data = data, initial_param =  initial_param, 
            range_x = plot_param.get('range_x',[np.min(data[0]), np.max(data[0])]), fmt = 'r-', label = 'Dopplershifted data', ecolor = 'k', capsize = 2) 
        return fig,ax,0
    else:
        F = HF_fitter(MASS = MASS, SCANS = SCANS, x = data[0], y = data[2], xerr = data[1], yerr = data[3], bunches = data[4], input_param = initial_param, 
            nrows = 2, ncols = nb_datasets, **plot_param)
        fig,ax = plt.subplots(nrows = 2, ncols = nb_datasets, figsize = (14,9), sharex = 'col', sharey = plot_param.get('sharey',False))
        F.init_all_models(F._input_param)
        fitter = F.fit_model(fig = fig, axs = ax, show_correl = False, llh = llh, **llh_kwargs)
        ax = param_plot(nb_datasets, fig, ax, data, F)
        if save:
            result_df = F.save_results()
            fig.savefig(f'{F._savepath_fig}\\scan-{SCANS[0]}.png', dpi = 300, bbox_inches = 'tight')
            return fig,ax,result_df
        return fig,ax,F.fitter
