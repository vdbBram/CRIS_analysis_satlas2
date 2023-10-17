import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import pymc3 as pm
from pymc3 import gp
import theano 
import theano.tensor as tt
import arviz as az

PATH = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_6-22\\Online Ag\\'
invalid_mass_folders = ['0', '1092', 'To_Serve', 'scanning.txt'] 
ref_mass = 109
ref_spin = 0.5
scale_redchi = False
name = 'GP_MAP_new'

def interpolate_refscans_centroid(path: str, mass: int, spin: float, scale_redchi: bool) -> pd.DataFrame:
    '''Reads and returns all centroids and median timestamp for the given mass and spin to a Pandas DataFrame 

        Parameters
        ----------
        path: str
            path to the folder with data, analysed data, etc
        mass: int
            mass number of the nucleus
        spin: float
            spin of the nuclear state
        scale_redchi: bool
            whether the error on the centroid should still be scaled with reduced chisquared or not

        Returns
        -------
        pd.DataFrame
        '''
    raw_data_path = f'{path}Data\\{mass}\\'
    centroid_path = f'{path}GP_final_analysis_output\\Fitted_HFS_param\\{mass}_{spin}\\'
    concat = []
    for scan in os.listdir(centroid_path):
        df_data = pd.read_csv(centroid_path + scan, delimiter = ';')
        median_time_scan = np.median(pd.read_csv(raw_data_path + '\\scan_' + str(scan[5:-4]) + '\\' + 'tagger_ds.csv' , sep = ';', header = None, names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp'])
        if not scale_redchi:
            temp_spin = str(spin).split('.')
            temp_spin = str(int(temp_spin[0])+int(2*int(temp_spin[1])/10)) + '_' + str(2)
            temp_data = {'timestamp': [median_time_scan],
                        'centroid': [(df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Value'][0])],
                        'centroid_err': (df_data[(df_data['Model'] == 'spin__'+str(temp_spin)) & (df_data['Parameter'] == 'centroid')]['Stderr'][0])}
        else:
            temp_data = {'timestamp': [median_time_scan],
                        'centroid': [df_data['Centroid'][0]],
                        'centroid_err': [df_data['Centroid'][1] * np.sqrt(df_data['redchi'][0])],
                        'red_chi': [df_data['redchi'][0]]}
        new_df = pd.DataFrame(temp_data)
        concat.append(new_df)        
    return_df = pd.concat(concat)
    return return_df

def get_median_timestamp(path: str) -> pd.DataFrame:
    '''Reads and returns all median timestamp for all scans to a Pandas DataFrame 

        Parameters
        ----------
        path: str
            path to the folder with data, analysed data, etc

        Returns
        -------
        pd.DataFrame
        '''
    data_path = f'{path}Data\\'
    concat_scan = []
    concat_mass = []
    for mass in os.listdir(data_path):
        if mass not in invalid_mass_folders:
            mass_path = f'{data_path}{mass}\\'
            for scan in os.listdir(mass_path):
                scan_path = f'{mass_path}{scan}\\'
                try:
                    median_time_scan = np.median(pd.read_csv(scan_path + 'tagger_ds.csv' , sep = ';', header = None, names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp'])
                except:
                    continue
                temp_data = {'timestamp': [median_time_scan]}
                scan_df = pd.DataFrame(temp_data)
                concat_scan.append(scan_df)   
            mass_df = pd.concat(concat_scan)
            concat_mass.append(mass_df)
    return_df = pd.concat(concat_mass)
    return return_df

if __name__ == '__main__':

    df = interpolate_refscans_centroid(PATH, ref_mass, ref_spin, scale_redchi)
    timestamp_min = np.min(df['timestamp'])
    df['timestamp'] = (df['timestamp'] - timestamp_min)/3600
    x = df['timestamp'].to_numpy()
    y = df['centroid'].to_numpy()
    avg_y = np.average(df['centroid'])
    y = y - avg_y
    yerr = df['centroid_err'].to_numpy()

    with pm.Model() as model:
        slow_trend_A = pm.HalfCauchy('A_slow', beta = 1000)
        length_scale_slow = pm.Gamma('length_slow', alpha = 160, beta = 1)
        cov_slow = slow_trend_A**2 *gp.cov.ExpQuad(1, length_scale_slow)
        gp_slow = pm.gp.Marginal(cov_func = cov_slow)

    #     good fast with period 20
    #     period_fast = pm.Gamma('period_fast', alpha = 20, beta = 1)
    #     length_fast = pm.Gamma('length_fast', alpha = 10, beta = 1)
    #     fast_trend_A = pm.Gamma('A_fast', alpha = 100, beta = 1)
    #     length_decay = pm.Gamma('length_decay', alpha = 10, beta = 0.5)
    #     cov_fast = fast_trend_A**2 *gp.cov.Periodic(1, period_fast, length_fast) * gp.cov.Matern52(1, length_decay)
    #     gp_fast = pm.gp.Marginal(cov_func = cov_fast)

        period_fast = pm.Gamma('period_fast', alpha = 12, beta = 0.5)
        length_fast = pm.Gamma('length_fast', alpha = 60, beta = 1)
        fast_trend_A = pm.Gamma('A_fast', alpha = 20, beta = 1)
        length_decay = pm.Gamma('length_decay', alpha = 60, beta = 1)
        cov_fast = fast_trend_A**2 *gp.cov.Periodic(1, period_fast, length_fast) * gp.cov.Matern52(1, length_decay)
        gp_fast = pm.gp.Marginal(cov_func = cov_fast)

        sigma_noise = pm.HalfNormal('sigma_noise', sigma = 10)
        cov_noise = gp.cov.WhiteNoise(sigma_noise)
        gp_noise = pm.gp.Marginal(cov_func = cov_noise)

        total_gp = gp_slow + gp_fast + gp_noise
        likelihood = total_gp.marginal_likelihood('centroid', X = x.reshape((len(x),1)), y = y, noise = yerr)

        trace = pm.sample(draws = 5000, tune = 500, chains = 4, cores = 2, return_inferencedata = True)
        MAP = pm.find_MAP(include_transformed = True)
    # with pm.Model() as model:
    #     slow_trend_A = pm.Gamma('A_slow', alpha = 80, beta = 1)
    #     length_scale_slow = pm.Gamma('length_slow', alpha = 160, beta = 1)
    #     cov_slow = slow_trend_A**2 *gp.cov.ExpQuad(1, length_scale_slow)
    #     gp_slow = pm.gp.Marginal(cov_func = cov_slow)

    #     period_fast = pm.Gamma('period_fast', alpha = 20, beta = 1)
    #     length_fast = pm.Gamma('length_fast', alpha = 10, beta = 1)
    #     fast_trend_A = pm.Gamma('A_fast', alpha = 15, beta = 1)
    #     length_decay = pm.Gamma('length_decay', alpha = 10, beta = 1)
    #     cov_fast = fast_trend_A**2 *gp.cov.Periodic(1, period_fast, length_fast) * gp.cov.Matern52(1, length_decay)
    #     gp_fast = pm.gp.Marginal(cov_func = cov_fast)

    #     sigma_noise = pm.HalfNormal('sigma_noise', sigma = 5)
    #     cov_noise = gp.cov.WhiteNoise(sigma_noise)
    #     gp_noise = pm.gp.Marginal(cov_func = cov_noise)

    #     total_gp = gp_slow + gp_fast + gp_noise

    #     likelihood = total_gp.marginal_likelihood('centroid', X = x.reshape((len(x),1)), y = y, noise = yerr)
    #     trace = pm.sample(draws = 5000, tune = 500, chains=4, cores=2, return_inferencedata=True)
    #     MAP = pm.find_MAP(include_transformed=True)

    fig = az.plot_trace(trace)

    print(az.summary(trace, hdi_prob = 0.95, round_to = 3))

    all_median_timestamp = get_median_timestamp(PATH)['timestamp'].sort_values()
    all_median_timestamp = (all_median_timestamp.to_numpy() - timestamp_min)/3600
    mu, var = total_gp.predict(all_median_timestamp[:, None], point = MAP, diag = True)
    mu = mu + avg_y
    sigma = np.sqrt(var)
    all_median_timestamp_og = (all_median_timestamp*3600) + timestamp_min

    import csv
    # with open(os.getcwd() + '\\' + name + '.csv', 'w', newline = '\n') as outfile:
    #     writer = csv.writer(outfile)
    #     for index in range(len(mu)):
    #         writer.writerow([all_median_timestamp_og[index],mu[index],sigma[index]])

    fig, ax = plt.subplots(1,1, figsize = (16,8))
    ax.errorbar(x, y+avg_y, yerr, fmt = 'r.', ecolor = 'k', capsize = 2, fillstyle = 'none', label = 'Experimental centroid values')
    ax.plot(all_median_timestamp, mu, c = "red", label = 'MAP')
    ax.fill_between(all_median_timestamp.ravel(), mu - sigma, mu + sigma, color = "red", alpha = 0.2, label = r'1-$\sigma$ interval')
    ax.fill_between(all_median_timestamp.ravel(), mu - 2*sigma, mu + 2*sigma, color = "red", alpha = 0.1, label = r'2-$\sigma$ interval')
    ax.set_xlabel("timestamp")
    ax.set_ylabel("Reference centroid [MHz]")
    plt.legend()
    plt.show()