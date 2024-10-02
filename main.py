from Dopplershift_analysis import *
from HF_fit_func import *
from HFFit_constants import *

'''To do the dopplershifting'''
##################################################
# scan = 5383
# manual = False #[21.5,26]
# D = Dopplershift_data(mass = 58, scan = scan, voltage_scanning = True, wn_channel = 'wavenumber_1', wn_bounds = [0,99999999999])
# D.bin_size_voltage = 2 # optional, it will default to 4 V but this is how you can change it
# D.bin_size_MHz = 60 # optional, it will default to 30 MHz but this is how you can change it
# data = D.extract_raw_data(devices_to_read = D._devices, path = D._PATH)
# data = D.AdvCutNoise(data = data, threshold = 0.15)
# # data = D.SimpleCutNoise(data = data)
# try: 
#     data = D.gate_tof(data = data)
# except:
#     data = D.gate_tof(data = data, manual = manual)
# data = D.filter_scatter(data = data, filename = 'iscool2', method = 'spline_ISCOOL', ISCOOL_voltage_multiplier = 10001.645) # apply filter to ISCOOL, for savgol you need to include window_length and polyorder. other options are spline_ISCOOL, avg
# data = D.calibrate_CRIS_voltagemeter(data = data, calibration_factor = 1.005030) 
# data = D.apply_wavenumber_correction(data = data) # you can change the ref wavenumber channel by using eg. D._ref_channel = 'wavenumber_3', defaults to 'wavenumber_2' 
# data = D.gate_wavenumber(data = data, wavenumber = D._wn_channel) # gate the wn values in case laser lost lock or wavemeter trips, you might wait for 7min if you dont because of np.arange, etc.
# binned_data = D.bin_data(data = data, freq_multiplier = 2) # this is for voltage scans, use bin_wm_data() with same arguments for laser scanning
# fig,ax = plt.subplots(figsize = (14,9))
# fig = D.plot(data = binned_data, fig = fig, ax = ax, save = True, save_format = 'png', fmt = 'r-', label = 'Dopplershifted data') # add save = True if you want to save
# D.save_data(binned_data) # saves in D._SAVEPATH, change it with D._SAVEPATH = new_save_path if you want, or change it in the class variables, defaults to the one in the class
# plt.show()

'''To do 1 state in 1 figure fitting'''
# ##################################################
# MASS = 63
# SCANS = ['5382']
# initial_param = {'spin__0_5':{'I':spins_gs[MASS], 'J':[3,3], 'ABC':ABC_gs[MASS], 
# 'centroid':centroid_gs[MASS], 'scale':0.003, 'variation_HFparam': {'Au':False, 'Bl': False, 'Bu': False, 'Cu': False, 'Cl': False}, 'set_values': {'FWHML':90, 'FWHMG': 124},
# # 'expressions':{'Amp3_2to3_2':[1,'Amp3_2to5_2'], 
# # 				'Amp5_2to3_2':[1,'Amp5_2to7_2'],'Amp5_2to5_2':[1,'Amp5_2to7_2'], 
# # 				'Amp7_2to5_2':[1,'Amp7_2to9_2'], 'Amp7_2to7_2':[1,'Amp7_2to9_2'],
# # 				'Amp9_2to7_2':[1,'Amp9_2to9_2']}
# # 'expressions':{'Amp5_2to5_2':[1,'Amp5_2to7_2'], 
# # 				'Amp7_2to7_2':[1,'Amp7_2to5_2']}
# }, 'bkg':{'background_values':[0.002], 'variation_param':{'p0':False}}#, 'set_values':{'p1':0.00001,'p0':0.025}}#, 'set_boundaries':{'p0':{'min':0.0001}}}
# # }, 'bkg':{'background_values':[0.00001,0.025], 'set_values':{'p1':0.00001,'p0':0.025}} # for 55: 5256
# }
# data = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS[0] + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy()
# # for omitting data in case it REALLY is needed
# # index_not_omit=[index for index, val in enumerate(data.T[0]) if ((val > -1500) and (val < 2200))]# or ((val > -100) and (val < 2100))]
# # data=data[index_not_omit]
# plot_param = {'ticker_locator':500,'range_x':[-1500,2500], 'range_y':[0,np.max(1.1*data.T[2]/data.T[4])]}
# # fig,ax,result_df = fit_1state_1fig(test=True, save=False, llh = False, data=[data], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# # 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# fig,ax,fitter = fitter(test=False, save=True, llh = False, data=[data], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# # freq_range = np.arange(-1000,2000,1)
# # testmodel = satlas2.HFS(I = 0.5, J = [3,3], A = [273,0], B = [0,0], C = [0,0], 
# #         df = 200, scale = 0.005, racah = False, name = 'HFSmodel for 61Cr', fwhml = 200)
# # plt.plot(freq_range,testmodel.f(freq_range),'g-', label = 'HFModel for 61Cr')
# # plt.legend()
# plt.show()
##################################################



##################################################

'''To do 1 state in 2 figure fitting'''
##################################################
# MASS = 119
# SCANS = ['4664','4665']
# initial_param = {'spin__3_5':{'I':spins_m1[MASS], 'J':[0.5,1.5], 'ABC':ABC_m1[MASS], 'centroid':centroid_m1[MASS], 'scale':10},
# 'bkg':{'background_values':[0.5,0.5],'background_bounds':[0], 'set_values':{'value0':0.0000001}}
# }
# data_1 = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS[0] + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy()
# data_2 = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS[1] + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy()
# data = np.concatenate((data_1,data_2)).T
# # for omitting data in case it REALLY is needed
# # index_not_omit=[index for index, val in enumerate(data_1[:,0]) if ((val > -11100) and (val > -12500))]# or ((val > -100) and (val < 2100))]
# # data_1=data_1[index_not_omit]
# plot_param = {'ticker_locator':500, 'sharey':'row','range_x':[[-19500,-16000],[18500,21000]],'range_y':[0,np.max(1.1*data[2]/data[4])]
# }
# # fig,ax,result_df = fit_1state_2fig(test=False, save=False, llh = False, data=[data_1,data_2], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# # 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# # plt.show()
# fig,ax,result_df = fitter(test=False, save=False, llh = False, data=[data_1,data_2], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# plt.show()
##################################################

'''To do 2 state in 2 figure fitting'''
##################################################
# MASS = 120
# SCANS = ['4612','4611']
# initial_param = {'spin__4':{'I':spins_m1[MASS], 'J':[0.5,1.5], 'ABC':ABC_m1[MASS], 'centroid':centroid_m1[MASS], 'scale':10, 'share_params': ['FWHMG', 'FWHML']},
# 'spin__7':{'I':spins_m2[MASS], 'J':[0.5,1.5], 'ABC':ABC_m2[MASS], 'centroid':centroid_m2[MASS], 'scale':20,'expressions':{'Amp13_2to15_2':[4/10,'Amp13_2to13_2']}},
# 'bkg':{'background_values':[0.5,0.5],'background_bounds':[0]}
# # ,'intensity_constraints':[['s1_Amp13_2__15_2',3,10,'s1_Amp13_2__13_2']]#,['s0_Amp7_2__9_2',1,10,'s0_Amp7_2__5_2']]#,['s1_Amp15_2__15_2',3,4,'s1_Amp15_2__17_2'],['s1_Amp15_2__13_2',6,10,'s1_Amp15_2__17_2'],
# # ['s0_Amp7_2__7_2',2,4,'s0_Amp7_2__5_2'],['s0_Amp7_2__9_2',7,10,'s0_Amp7_2__5_2'],['s1_Amp13_2__13_2',2,4,'s1_Amp13_2__11_2'],['s1_Amp13_2__15_2',7,10,'s1_Amp13_2__11_2']]
# # ,'intensity_constraints':[['s0_Amp7_2__5_2',1,1,'s0_Amp7_2__7_2'],['s0_Amp7_2__5_2',1,1,'s0_Amp7_2__7_2'],['s0_Amp9_2__9_2',1,2,'s0_Amp9_2__11_2'],['s1_Amp13_2__15_2',1,2,'s1_Amp13_2__13_2']]
# }
# data_1 = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS[0] + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy()
# data_2 = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS[1] + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy()
# data = np.concatenate((data_1,data_2)).T
# # for omitting data in case it REALLY is needed
# # index_not_omit=[index for index, val in enumerate(data_1[0]) if ((val < -9000) or (val < -8000))]# or ((val > -100) and (val < 2100))]
# # data_1=data_1[:,index_not_omit] 
# plot_param = {'ticker_locator':500, 'sharey':'row','range_x':[[-15000,-12000],[12000,14500]],'range_y':[0,np.max(1.1*data[2]/data[4])]}
# # fig,ax,result_df = fit_2state_2fig(test=False, save=False, llh = False, data=[data_1,data_2], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# # 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# fig,ax,result_df = fitter(test=False, save=False, llh = False, data=[data_1,data_2], initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# plt.show()
##################################################


# special case
######################
# MASS = 116
# SCANS = '4603'
# initial_param = {'I':[4,7], 'J':[0.5,1.5], 'ABC':[ABC_m1[MASS],ABC_m2[MASS]], 'centroid':[centroid_m1[MASS],centroid_m2[MASS]], 'scale':[10,20], 'background_values':[0.5,0.5], 
# 'background_bounds': [0]
# ,'intensity_constraints':[['s0_Amp9_2__9_2',3,4,'s0_Amp9_2__11_2'],['s0_Amp9_2__7_2',6,10,'s0_Amp9_2__11_2'],['s1_Amp15_2__15_2',3,4,'s1_Amp15_2__17_2'],['s1_Amp15_2__13_2',6,10,'s1_Amp15_2__17_2'],
# ['s0_Amp7_2__7_2',2,4,'s0_Amp7_2__5_2'],['s0_Amp7_2__9_2',7,10,'s0_Amp7_2__5_2'],['s1_Amp13_2__13_2',2,4,'s1_Amp13_2__11_2'],['s1_Amp13_2__15_2',7,10,'s1_Amp13_2__11_2']]
# # ,'intensity_constraints':[['s0_Amp7_2__5_2',1,1,'s0_Amp7_2__7_2'],['s0_Amp7_2__5_2',1,1,'s0_Amp7_2__7_2'],['s0_Amp9_2__9_2',1,2,'s0_Amp9_2__11_2'],['s1_Amp13_2__15_2',1,2,'s1_Amp13_2__13_2']]
# }
# data = pd.read_csv(DATA_FOLDER + 'Analysis_output\\' + str(MASS) + '\\' + SCANS + '.csv', delimiter = ';',header = 0, names = data_columns).to_numpy().T
# # for omitting data in case it REALLY is needed
# index_not_omit=[index for index, val in enumerate(data[0]) if ((val < -8000) or (val > 10000))]# or ((val > -100) and (val < 2100))]
# data=data[:,index_not_omit]
# plot_param = {'range_x1':[-12500,-9500],'range_y1':[0,4],'range_x2':[12000,14000], 'range_y1':[-1,30], 'range_y2':[0,np.max(1.1*data[2]/data[4])],'ticker_locator':500, 'sharey':'row'}
# fig,ax = fit_2state_2fig(test=False, save=True, data=data, initial_param=initial_param, plot_param=plot_param, MASS=MASS, mass_ref=mass_ref, I_ref=I_ref, SCANS=SCANS, 
# 	DATA_FOLDER=DATA_FOLDER, SAVEPATH_FIG=SAVEPATH_FIG, SAVEPATH_PARAM=SAVEPATH_PARAM)
# plt.show()
