import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from uncertainties import unumpy
from uncertainties import ufloat as uf

def GP_refscans_centroid(path: str, name: str) -> pd.DataFrame:
    '''Returns the GPR of the reference centroids

    Parameters
    ----------
    path: string
        path to the GP .csv file
    name: string
        name of the GP .csv file

    Returns
    -------
    DataFrame
    '''
    df = pd.read_csv(path + name +'.csv', names = ['median_timestamp', 'MAPcentroid', 'MAPcentroid_err'])
    df['median_timestamp_copy'] = df['median_timestamp']
    df = df.set_index(['median_timestamp'])
    df = df.sort_index()
    return df

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
    # scans = []
    time_path = path + 'Data\\' + str(mass_ref) + '\\'
    path = path + 'GP_final_analysis_output\\Fitted_HFS_param\\' + str(mass_ref) + '_' + str(I_ref) + '\\' # make sure it can handle float spins
    for scan in os.listdir(path):
        times.append(np.median(pd.read_csv(time_path + 'scan_' + scan[5:-4] +'\\tagger_ds.csv', delimiter = ';', names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp']))
        df_data = pd.read_csv(path + scan, delimiter = ';')
        temp_spin = str(I_ref).split('.')
        if len(temp_spin) == 1:
            ref_centroid_val = (df_data[(df_data['Model'] == 'spin__'+str(I_ref)) & (df_data['Parameter'] == 'centroid')]['Value'][0])
            ref_centroid_std = (df_data[(df_data['Model'] == 'spin__'+str(I_ref)) & (df_data['Parameter'] == 'centroid')]['Stderr'][0])
        else:
            temp_spin = str(int(temp_spin[0])+int(2*int(temp_spin[1])/10)) + '_' + str(2)
            ref_centroid_val = (df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Value'][0])
            ref_centroid_std = (df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Stderr'][0])
        ref_centroid.append(uf(ref_centroid_val,ref_centroid_std))
        scans.append(int(scan[5:-4]))
    # print(scans)
    # print(unumpy.nominal_values(ref_centroid))
    # print(unumpy.std_devs(ref_centroid))
    # print((times-np.min(times))/3600)
    times.append(2*times[-1])
    # scans.append(int(scan[5:-4])+200)
    ref_centroid.append(uf(ref_centroid_val,ref_centroid_std))
    # print(np.average(unumpy.nominal_values(ref_centroid)), np.std(unumpy.nominal_values(ref_centroid)))
    ref_centroid_spline = interp1d(times, unumpy.nominal_values(ref_centroid))
    if plot:
        test_x = np.arange(min(times), max(times), 100)
        test_y = ref_centroid_spline(test_x)
        plt.errorbar((times-np.min(times))/3600, unumpy.nominal_values(ref_centroid), unumpy.std_devs(ref_centroid), fmt = '.', ecolor = 'k', capsize = 2, label = 'Reference centroid values')
        plt.plot((test_x-np.min(times))/3600, test_y, label = 'Spline')
        # plt.plot(scans,unumpy.nominal_values(ref_centroid))
        plt.xlabel('timestamp [h]')
        plt.ylabel('ref_centroid [MHz]')
        plt.legend(loc = 'upper right')
        plt.show()
    return ref_centroid_spline