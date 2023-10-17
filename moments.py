import os
import sys
import numpy as np
import pandas as pd
from uncertainties import ufloat as uf
from numpy.typing import ArrayLike
from typing import Union
# sys.path.append('C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_6-22\\Online Ag\\Analysis_satlas2')
from interpolate_refscans_centroid import *

class NuclearState:

    _CRIS_DATA_PATH = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_..-2021\\Data\\'
    _CRIS_ANALYSED_DATA_PATH = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_..-2021\\Final_analysis_output\\Fitted_HFS_param\\'
    _GP_path = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_6-22\\Online Ag\\Analysis_satlas2\\'
    _params_ = ['Al', 'Au', 'Bu', 'centroid', 'FWHMG', 'FWHML']
    _column_names = ['Source','Model','Parameter','Value','Stderr','Minimum','Maximum','Expression','Vary']
    _scan_nbs = []
    _scan_median_time = []

    def __init__(self, MASS: int, I: float) -> None:
        '''initialises class to calculate moments of a certain state of spin I of an isotope with massnumber MASS

        Parameters
        ----------
        MASS: int
            massnumber of the isotope
        I: float
            spin of the nuclear state
        '''
        self._MASS = MASS
        self._I = I
        self._I_ref = 0.5 #1 # for dipole moment ref
        self._mass_ref_u = 109 #110 # for dipole moment ref
        self._mass_ref_centroid = 109 # uf(108.9047558, 0.0000014) # in u for charge radius reference
        self._I_ref_centroid = 0.5
        self._separate_HFparam_dict = self.calculate_separate_HFparam_dict()
        self._HFparam_dict = self.calculate_HFparam_dict()
        self._A_ref =  self.calculate_HF_ref(self._mass_ref_u, self._I_ref)# in MHz
        self._Mu_ref = uf(0.1306906,0.0000002) #uf(2.7111, 0.0010) # in mu_n
        self._B_ref = uf(426, 18) # in MHz
        self._Q_ref = uf(1.4,0.20) # in barn from IBerkes 1984 from 110m but with updated EFG from stefanos (and correction for sternheimer shielding effects)
        self._M = uf(1956, 0) * 1000 #360
        self._F = uf(-4300, 0) #300
        self._K = 0.976
        self._Nuc_Prop_dict = self.calculate_Nuc_Prop_dict()

    def calculate_HF_ref(self, ref_mass: int, ref_I: float) -> uf:
        '''Calculates the Alower parameter of the reference state.

        Parameters
        ----------
        ref_mass: int
            massnumber of reference isotope
        ref_I: float
            spin of reference state

        Returns
        -------
        ufloat
        '''
        self._separate_HFparam_dict_ref = dict()
        for scanfile in os.listdir(self._CRIS_ANALYSED_DATA_PATH + str(ref_mass) + '_' + str(ref_I) + '\\'):
            single_state_HFparam_dict_ref = self.readData1State(scanfile, ref_mass, ref_I)
            for key in single_state_HFparam_dict_ref.keys():
                try:
                    self._separate_HFparam_dict_ref[key].append(single_state_HFparam_dict_ref[key])
                except:
                    self._separate_HFparam_dict_ref[key] = [single_state_HFparam_dict_ref[key]]
        return self.weightedAvg(self._separate_HFparam_dict_ref['Al'], self._separate_HFparam_dict_ref['Al_err'])
    
    def readData1State(self, scanfile: str, mass: int, I: float, median_times: ArrayLike = 0) -> dict:
        '''reads the fit parameters for a single scan file and returns a dictionary with keys the parameter names and parameter names + _err

        Parameters
        ----------
        scanfile: str
            name of the file where the HF parameters are saved
        mass: int
            massnumber of the isotope
        I: float
            spin of the state

        Returns
        -------
        dict
        '''
        df_data = pd.read_csv(self._CRIS_ANALYSED_DATA_PATH + '\\' + str(mass) + '_' + str(I) + '\\' + scanfile, delimiter = ';', names = self._column_names, header=0)
        models = df_data['Model'].unique()
        return_dict = dict()
        if len(str(I).split('.')) == 1:
            model_name = 'spin__' + str(I).split('.')[0]
        else:
            model_name = 'spin__' + str(int(2*I)) + '_2'
        for param in self._params_:
            return_dict[param] = df_data[(df_data['Model'] == model_name) & (df_data['Parameter'] == param)]['Value'].iloc[0]
            return_dict[param + '_err'] = df_data[(df_data['Model'] == model_name) & (df_data['Parameter'] == param)]['Stderr'].iloc[0]
            if param == 'centroid' and median_times != 0 and not (mass == self._mass_ref_centroid and I == self._I_ref_centroid):
                errGP = self.calculate_errGPCentroid(median_times)
                return_dict['centroid_err'] = np.sqrt((return_dict['centroid_err']*return_dict['centroid_err']) + (errGP*errGP))
        return return_dict

    def calculate_errGPCentroid(self, median_times: ArrayLike) -> ArrayLike:
        '''Returns the error from the Gaussian Process Regression on the reference centroids for a particular set of median times

        Parameters
        ----------
        median_times: ArrayLike
            what time the error on the GP should be returned

        Returns
        -------
        ArrayLike
        '''
        if type(median_times) != type([1,2,3]):
            GP_err = np.median(GP_refscans_centroid(self._GP_path, 'GP_MAP')['MAPcentroid_err'].loc[np.isclose(GP_refscans_centroid(self._GP_path, 'GP_MAP')['median_timestamp_copy'], median_times)])
            return GP_err
        else:
            GP1_err = np.median(GP_refscans_centroid(self._GP_path, 'GP_MAP')['MAPcentroid_err'].loc[np.isclose(GP_refscans_centroid(self._GP_path, 'GP_MAP')['median_timestamp_copy'], median_times[0])])
            GP2_err = np.median(GP_refscans_centroid(self._GP_path, 'GP_MAP')['MAPcentroid_err'].loc[np.isclose(GP_refscans_centroid(self._GP_path, 'GP_MAP')['median_timestamp_copy'], median_times[1])])
            GP_err = np.sqrt((GP1_err*GP1_err) + (GP2_err*GP2_err))
            return GP_err
    
    def calculate_separate_HFparam_dict(self) -> dict:
        '''calculates the HF parameters and returns a dictionary with keys the parameter names and parameter names + _err, the values are lists with all the parameter values/errors 

        Returns
        -------
        dict
        '''
        separate_HFparam_dict = dict()
        for scanfile in os.listdir(self._CRIS_ANALYSED_DATA_PATH + str(self._MASS) + '_' + str(self._I) + '\\'): #ADD LISTDIR
            self._scan_nbs.append(scanfile[5:-4])
            try:
                median_scan_time = np.median(pd.read_csv(self._CRIS_DATA_PATH + str(self._MASS) + '\\scan_' + str(scanfile[5:-4]) + '\\' + 'tagger_ds.csv' , sep = ';', header = None, names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp'])
                self._scan_median_time.append(median_scan_time)
                single_state_HFparam_dict = self.readData1State(scanfile, self._MASS, self._I, median_scan_time)
            except:
                median_scan_time_1 = np.median(pd.read_csv(self._CRIS_DATA_PATH + str(self._MASS) + '\\scan_' + str(scanfile[5:9]) + '\\' + 'tagger_ds.csv' , sep = ';', header = None, names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp'])
                median_scan_time_2 = np.median(pd.read_csv(self._CRIS_DATA_PATH + str(self._MASS) + '\\scan_' + str(scanfile[10:14]) + '\\' + 'tagger_ds.csv' , sep = ';', header = None, names = ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'])['timestamp'])
                self._scan_median_time.append([median_scan_time_1,median_scan_time_2])
                single_state_HFparam_dict = self.readData1State(scanfile, self._MASS, self._I, [median_scan_time_1,median_scan_time_2])
            for key in single_state_HFparam_dict.keys():
                try:
                    separate_HFparam_dict[key].append(single_state_HFparam_dict[key])
                except:
                    separate_HFparam_dict[key] = [single_state_HFparam_dict[key]]
        return separate_HFparam_dict

    def calculate_HFparam_dict(self) -> dict:
        '''Calculates the weighted average of all values measured for the HF parameters and returns a dictionary with this average and weighted error

        Returns
        -------
        dict
        '''
        HFparam_dict = dict()
        for key in self._params_:
            HFparam_dict[key] = self.weightedAvg(self._separate_HFparam_dict[key], self._separate_HFparam_dict[key + '_err'], self._separate_HFparam_dict.get('redchi',1))
        return HFparam_dict

    def calculate_Nuc_Prop_dict(self) -> dict:
        '''Calculates the nuclear properties and returns a dictionary with these values and weighted errors

        Returns
        -------
        dict
        '''
        Nuc_Prop_dict = dict()
        Nuc_Prop_dict['mu'] = self.calculate_Mu()
        Nuc_Prop_dict['q'] = self.calculate_Q()
        Nuc_Prop_dict['charge_radius'] = self.calculate_Charge_Radius()
        return Nuc_Prop_dict

    def weightedAvg(self, data: ArrayLike, data_err: ArrayLike, redchi: ArrayLike = 1, scale: bool = True) -> uf:
        '''Calculates the weighted average and weighted error scaled with reduced chi

        Parameters
        ----------
        data: ArrayLike
            the values for the average
        data_err: ArrayLike
            the errors for the average
        redchi: ArrayLike, default: 1
            reduced chi value
        scale: boolean, default: True
            whether to scale the average of this data with the reduced chisquared of this average and data

        Returns
        -------
        ufloat
        '''
        #### test this scaling
        variance = (np.array(data_err)**2) * (np.array(redchi))
        if 0 in variance:
            return 0
        weight = 1 / variance
        stand_dev = np.sqrt(1 / sum(weight))
        avg = sum(data * weight) / sum(weight)
        if scale and len(data) > 1:
            res_sq = ((avg-np.array(data))/np.array(data_err)) * ((avg-np.array(data))/np.array(data_err))
            red_res = np.sum(res_sq) / (len(data)-1)
            stand_dev = stand_dev * np.sqrt(red_res)
        return uf(avg, stand_dev)

    def calculate_Mu(self) -> Tuple[uf,float]:
        '''Calculates the nuclear magnetic dipole moment wrt the reference isotope

        Returns
        -------
        tuple of ufloat and float
        '''
        if self._I == 0:
            return uf(0, 0)
        val_stat_err = (self._I/self._I_ref) * (self._HFparam_dict['Al']/self._A_ref) * self._Mu_ref.nominal_value
        sys_err = val_stat_err.std_dev + ((self._I/self._I_ref) * (self._HFparam_dict['Al'].nominal_value/self._A_ref.nominal_value) * self._Mu_ref).std_dev
        return val_stat_err,sys_err

    def calculate_Q(self) -> Tuple[uf,float]:
        '''Calculates the nuclear electric quadrupole moment wrt the reference isotope

        Returns
        -------
        tuple of ufloat and float
        '''
        if self._I < 1:
            return uf(0, 0),0
        val_stat_err = (self._HFparam_dict['Bu']/self._B_ref.nominal_value) * self._Q_ref.nominal_value
        sys_err = val_stat_err.std_dev + ((self._HFparam_dict['Bu'].nominal_value/self._B_ref) * self._Q_ref).std_dev
        return val_stat_err,sys_err
    
    def calculate_Charge_Radius(self) -> Tuple[uf,float]:
        '''Calculates the change in mean-squared charge radius wrt the reference isotope

        Returns
        -------
        tuple of ufloat and float
        '''
        if self._MASS == self._mass_ref_centroid and self._I == self._I_ref_centroid:
            return uf(0,0),0
        else:
            try:
                val_stat_err = ((self._HFparam_dict.get('centroid',self._HFparam_dict['centroid']) / (self._F*self._K).nominal_value) - self._M.nominal_value * ((self._MASS-self._mass_ref_centroid) / (self._MASS*self._mass_ref_centroid)) / (self._F*self._K).nominal_value)
                sys_err = val_stat_err.std_dev + ((self._HFparam_dict.get('centroid',self._HFparam_dict['centroid']).nominal_value / (self._F*self._K)) - self._M * ((self._MASS-self._mass_ref_centroid) / (self._MASS*self._mass_ref_centroid)) / (self._F*self._K)).std_dev
                return val_stat_err,sys_err
            except:
                val_stat_err = ((self._HFparam_dict.get('Centroid',self._HFparam_dict['Centroid']) / (self._F*self._K).nominal_value) - self._M.nominal_value * ((self._MASS-self._mass_ref_centroid) / (self._MASS*self._mass_ref_centroid)) / (self._F*self._K).nominal_value)
                sys_err = val_stat_err.std_dev + ((self._HFparam_dict.get('Centroid',self._HFparam_dict['Centroid']).nominal_value / (self._F*self._K)) - self._M * ((self._MASS-self._mass_ref_centroid) / (self._MASS*self._mass_ref_centroid)) / (self._F*self._K)).std_dev
                return val_stat_err,sys_err
