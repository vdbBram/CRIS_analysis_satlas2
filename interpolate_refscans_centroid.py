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

def spline_refscans_centroid(path: str, plot: bool = False) -> interp1d:
	'''interpolates the reference scans following a linear spline

        Parameters
        ----------
        path: str
            path to the fitted HF parameters
        plot: bool
            True if you want to see the plot
        Returns
        -------
        interp1d
        '''
	ref_centroid = list()
	scans = list()
	for scan in os.listdir(path):
		scans.append(int(scan[5:-4]))
		df_data = pd.read_csv(path + scan, delimiter = ';')
		if 's0' in df_data.columns[1]:
			add_key = 's0_'
		elif 's1' in df_data.columns[1]:
			add_key = 's1_'
		else:
			add_key = ''
		ref_centroid.append(uf(df_data[add_key + 'Centroid'][0], df_data[add_key + 'Centroid'][1]))
	ref_centroid_spline = interp1d(scans, unumpy.nominal_values(ref_centroid))
	if plot:
		test_x = np.arange(min(scans), max(scans), 1)
		test_y = ref_centroid_spline(test_x)
		plt.plot(scans, unumpy.nominal_values(ref_centroid),'.', label = 'Reference centroid values')
		plt.plot(test_x, test_y, label = 'Spline')
		plt.xlabel('Scan')
		plt.ylabel('ref_centroid [MHz]')
		plt.legend()
		plt.show()
	return ref_centroid_spline

def GP_refscans_centroid(path, plot, **kwargs):
	raise NotImplementedError('Bram should implement this at some point, the lazy fk xddd')