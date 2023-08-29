import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from uncertainties import unumpy
from uncertainties import ufloat as uf

def interpolate_refscans_centroid(path: str, method: str, plot: bool) -> interp1d:
	'''interpolates the reference scans following the supplied method

        Parameters
        ----------
        path: str
            path to the fitted HF parameters
        method: str
            method that should be used
        plot: bool
            True if you want to see the plot
        Returns
        -------
        interp1d
        '''
	if method.lower() == 'spline':
		return spline_refscans_centroid(path = path, plot = plot)
	elif method.lower() == 'gp':
		return GP_refscans_centroid(path = path, plot = plot, **kwargs)
	else:
		raise RuntimeError('implement that yourself, u lazy fk')

def spline_refscans_centroid(path: str, mass_ref: int, I_ref: float, plot: bool = False) -> interp1d:
	'''interpolates the reference scans following a linear spline

        Parameters
        ----------
        path: str
            path to the fitted HF parameters
        mass_ref: int
        	mass number of the reference isotope
        I_ref: float
        	Spin of the reference isotope
        plot: bool
            True if you want to see the plot
        Returns
        -------
        interp1d
        '''
	ref_centroid = []
	times = []
	time_path = path + 'Data\\' + str(mass_ref) + '\\'
	path = path + '\\Final_analysis_output\\Fitted_HFS_param\\' + str(mass_ref) + '_' + str(I_ref) + '\\' # make sure it can handle float spins
	for scan in os.listdir(path):
		times.append(np.median(pd.read_csv(time_path + 'scan_' + scan[5:-4] +'\\tagger_ds.csv', delimiter = ';', names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp']))
		df_data = pd.read_csv(path + scan, delimiter = ';')
		temp_spin = str(I_ref).split('.')
		if len(temp_spin) == 1:
			ref_centroid_val = (df_data[(df_data['Model'] == 'spin__'+str(I_ref)) & (df_data['Parameter'] == 'centroid')]['Value'][0])
			ref_centroid_std = (df_data[(df_data['Model'] == 'spin__'+str(I_ref)) & (df_data['Parameter'] == 'centroid')]['Stderr'][0])
		else:
			temp_spin = temp_spin[0] + '_' + temp_spin[1]
			ref_centroid_val = (df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Value'][0])
			ref_centroid_std = (df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Stderr'][0])
		ref_centroid.append(uf(ref_centroid_val,ref_centroid_std))
	times.append(2*times[-1])
	ref_centroid.append(uf(ref_centroid_val,ref_centroid_std))
	ref_centroid_spline = interp1d(times, unumpy.nominal_values(ref_centroid))
	if plot:
		test_x = np.arange(min(times), max(times), 100)
		test_y = ref_centroid_spline(test_x)
		plt.plot(times, unumpy.nominal_values(ref_centroid),'.', label = 'Reference centroid values')
		plt.plot(test_x, test_y, label = 'Spline')
		plt.xlabel('Scan')
		plt.ylabel('ref_centroid [MHz]')
		plt.legend(loc = 'upper right')
		plt.show()
	return ref_centroid_spline

def GP_refscans_centroid(path, plot, **kwargs):
	raise NotImplementedError('Bram should implement this at some point, the lazy fk xddd')