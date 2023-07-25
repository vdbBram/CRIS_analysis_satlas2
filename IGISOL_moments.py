import os
import numpy as np
import pandas as pd
from uncertainties import ufloat as uf
from moments import NuclearState

class IGISOL_NuclearState(NuclearState):

    _ANALYSED_DATA_PATH = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\2021-2022\\JYV\\Silver 2020-2021\\Analysis\\Final_analysis_output\\Fitted_HFS_param\\'
    _CRIS_ANALYSED_DATA_PATH = 'C:\\Users\\u0148746\\OneDrive - KU Leuven\\Documents\\PhD\\CRIS\\Ag_6-22\\Online Ag\\Final_analysis_output\\Fitted_HFS_param\\'
    _params_ = ['Al', 'Au', 'Bu', 'Centroid', 'TotalFWHM']
    def __init__(self, MASS: int, I: float) -> None: 
        '''initialises class to calculate moments of a certain state of spin I of an isotope with massnumber MASS for my kankled datastructure of my masterthesis

        Parameters
        ----------
        MASS: int
            massnumber of the isotope
        I: float
            spin of the nuclear state
        '''
        self._I_ref_centroid = 0.5 
        self._mass_ref_centroid = 109
        self._ref_centroid = self.calculate_avg_centroid_ref(self._I_ref_centroid, self._mass_ref_centroid)
        super().__init__(MASS, I)

    def readData1State(self, scanfile: str, mass: int, I: float) -> tuple:
        '''reads the fit parameters for a single scan file wit given mass and spin I and returns a tuple of a dictionary with keys the parameter names and parameter names + _err and the reduced chi squared

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
        tuple
        '''
        df_data = pd.read_csv(self._CRIS_ANALYSED_DATA_PATH + '\\' + str(mass) + '_' + str(I) + '\\' + scanfile, delimiter = ';')
        red_chi = df_data['redchi'][0]
        return_dict = dict()
        if 's0' in df_data.columns[1]:
            add_key = 's0_'
        elif 's1' in df_data.columns[1]:
            add_key = 's1_'
        else:
            add_key = ''
        for param in self._params_:
            return_dict[param] = df_data[add_key + param][0]
            return_dict[param + '_err'] = df_data[add_key + param][1]
        return return_dict, red_chi

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
        for scanfile in os.listdir(self._CRIS_ANALYSED_DATA_PATH + str(ref_mass) + '_' + str(ref_I) + '\\'): #ADD LISTDIR
            single_state_HFparam_dict_ref, red_chi = self.readData1State(scanfile, ref_mass, ref_I)
            try:
                self._separate_HFparam_dict_ref['redchi'].append(red_chi)
            except:
                self._separate_HFparam_dict_ref['redchi'] = [red_chi]
            for key in single_state_HFparam_dict_ref.keys():
                try:
                    self._separate_HFparam_dict_ref[key].append(single_state_HFparam_dict_ref[key])
                except:
                    self._separate_HFparam_dict_ref[key] = [single_state_HFparam_dict_ref[key]]
        return self.weightedAvg(self._separate_HFparam_dict_ref['Al'], self._separate_HFparam_dict_ref['Al_err'], self._separate_HFparam_dict_ref['redchi'])

    def calculate_avg_centroid_ref(self, ref_I: float, ref_mass: int) -> uf:
        '''Calculates the average centroid of the reference state.

        Parameters
        ----------
        ref_I: float
            spin of reference state
        ref_mass: int
            massnumber of reference isotope

        Returns
        -------
        ufloat
        '''
        centroid_dict = dict()
        for scanfile in os.listdir(self._ANALYSED_DATA_PATH + 'mass-' + str(ref_mass) + 'spin-' + str(ref_I) + '\\'):
            if len(scanfile.split('_')) == 1:
                df_data = pd.read_csv(self._ANALYSED_DATA_PATH + 'mass-' + str(ref_mass) + 'spin-' + str(ref_I) + '\\' + scanfile, delimiter = ';')
                if 's0' in df_data.columns[1]:
                    add_key = 's0_'
                elif 's1' in df_data.columns[1]:
                    add_key = 's1_'
                else:
                    add_key = ''
                try:
                    centroid_dict['value'].append(df_data[add_key + 'Centroid'][0])
                    centroid_dict['err'].append(df_data[add_key + 'Centroid'][1])
                    centroid_dict['redchi'].append(df_data[add_key + 'redchi'][0])
                except:
                    centroid_dict['value'] = [df_data[add_key + 'Centroid'][0]]
                    centroid_dict['err'] = [df_data[add_key + 'Centroid'][1]]
                    centroid_dict['redchi'] = [df_data[add_key + 'redchi'][0]]
        ref_centroid = self.weightedAvg(centroid_dict['value'], centroid_dict['err'], centroid_dict['redchi'])
        return ref_centroid     
    
    def readData1State_IGISOL(self, scanfile: str, mass: int, I: float) -> tuple: 
        '''Reads the data for 1 state in my IGISOL kankled datastructure and returns a tuple with a dictionary with the parameters and values + errors and the reduced chi squared

        Parameters
        ----------
        scanfile: str
            filename of the scan 
        mass: int
            massnumber of the isotope
        I: float
            spin of the state

        Returns
        -------
        tuple
        '''
        df_data = pd.read_csv(self._ANALYSED_DATA_PATH + 'mass-' + str(mass) + 'spin-' + str(I) + '\\' + scanfile, delimiter = ';')
        red_chi = df_data['redchi'][0]
        return_dict = dict()
        if 's0' in df_data.columns[1]:
            add_key = 's0_'
        elif 's1' in df_data.columns[1]:
            add_key = 's1_'
        else:
            add_key = ''
        for param in self._params_:
            if param == 'Centroid':
                tcen = uf(df_data[add_key + param][0],df_data[add_key + param][1]) - self._ref_centroid
                return_dict[param] = tcen.nominal_value
                return_dict[param + '_err'] = tcen.std_dev
            else:
                return_dict[param] = df_data[add_key + param][0]
                return_dict[param + '_err'] = df_data[add_key + param][1]
        if mass == 109 and I == 0.5:
            return_dict['Centroid'] = 0
            return_dict['Centroid' + '_err'] = 0
        return return_dict, red_chi
    
    def calculate_separate_HFparam_dict(self) -> dict:
        '''calculates the HF parameters and returns a dictionary with keys the parameter names and parameter names + _err, the values are lists with all the parameter values/errors 
        and the reduced chi squared under key 'redchi'

        Returns
        -------
        dict
        '''
        separate_HFparam_dict = dict()
        try:
            for scanfile in os.listdir(self._ANALYSED_DATA_PATH + 'mass-' + str(self._MASS) + 'spin-' + str(self._I) + '\\'):
                if scanfile.split('_')[-1][0] != 'd':
                    single_state_HFparam_dict, red_chi = self.readData1State_IGISOL(scanfile, self._MASS, self._I)
                    try:
                        separate_HFparam_dict['redchi'].append(red_chi)
                    except:
                        separate_HFparam_dict['redchi'] = [red_chi]
                    for key in single_state_HFparam_dict.keys():
                        try:
                            separate_HFparam_dict[key].append(single_state_HFparam_dict[key])
                        except:
                            separate_HFparam_dict[key] = [single_state_HFparam_dict[key]]
            return separate_HFparam_dict
        except:
            print('Have you still not fixed the structure of your data analysis of IGISOL?... lazy fk')
            return self.special_calculate_separate_HFparam_dict()

    def special_calculate_separate_HFparam_dict(self) -> dict:
        '''calculates the HF parameters and returns a dictionary with keys the parameter names and parameter names + _err, the values are lists with all the parameter values/errors 

        Returns
        -------
        dict
        '''
        separate_HFparam_dict = dict()
        for folder in os.listdir(self._ANALYSED_DATA_PATH):
            if 'mass-' + str(self._MASS) in folder and '-' + str(self._I) in folder:
                path = self._ANALYSED_DATA_PATH + folder + '\\'
                if float(folder.split('-')[-1]) == float(self._I):
                    other_I = float(folder.split('-')[-2])
                else:
                    other_I = float(folder.split('-')[-1])
                break
        for scanfile in os.listdir(path):   
            if scanfile.split('_')[-1][0] != 'd':
                try:
                    single_state_HFparam_dict, red_chi = self.special_readData1State_IGISOL(scanfile, self._MASS, self._I, other_I)
                except:
                    single_state_HFparam_dict, red_chi = self.special_readData1State_IGISOL(scanfile, self._MASS, self._I, int(other_I))
                try:
                    separate_HFparam_dict['redchi'].append(red_chi)
                except:
                    separate_HFparam_dict['redchi'] = [red_chi]
                for key in single_state_HFparam_dict.keys():
                    try:
                        separate_HFparam_dict[key].append(single_state_HFparam_dict[key])
                    except:
                        separate_HFparam_dict[key] = [single_state_HFparam_dict[key]]
        return separate_HFparam_dict

    def special_readData1State_IGISOL(self, scanfile: str, mass: int, return_I: float, I: float) -> tuple:
        '''reads the fit parameters for the special cases where i had put 2 states in 1 .csv file for saving...
        It returns a tuple of a dictionary with keys the parameter names and parameter names + _err, the values are lists with all the parameter values/errors and, reduced chi squared 

        Returns
        -------
        dict
        '''
        if return_I < I:
            df_data = pd.read_csv(self._ANALYSED_DATA_PATH + 'mass-' + str(mass) + 'spin-' + str(return_I) + '-' + str(I) + '\\' + scanfile, delimiter = ';')
        else:
            df_data = pd.read_csv(self._ANALYSED_DATA_PATH + 'mass-' + str(mass) + 'spin-' + str(I) + '-' + str(return_I) + '\\' + scanfile, delimiter = ';')
        red_chi = df_data['redchi'][0]
        return_dict = dict()
        if return_I < I:
            add_key = 's0_'
        else:
            add_key = 's1_'
        for param in self._params_:
            if param == 'Centroid':
                tcen = uf(df_data[add_key + param][0],df_data[add_key + param][1]) - self._ref_centroid
                return_dict[param] = tcen.nominal_value
                return_dict[param + '_err'] = tcen.std_dev
            else:
                return_dict[param] = df_data[add_key + param][0]
                return_dict[param + '_err'] = df_data[add_key + param][1]
        return return_dict, red_chi