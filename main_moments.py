from plotting import *

P = Plotting()
fig,ax = plt.subplots(figsize = (12,9))

# low spin odd-odd Ag
# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {102:['2+lit'],104:['2+lit'],106:['1+','1+lit'],108:['1+'],110:['1+'],112:['2-'],114:['1+IG'],
# 	116:['1-','1-IG'],118:['0-','0-IG'],120:['0-','0-IG'],122:['0-']}, 
# 	param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = True, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4)

# # charge radii low spin Ag
# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {106:['1+'],107:['0.5-','0.5-IG'],108:['1+', '6+'],109:['0.5-IG', '0.5-'],110:['1+','6+'],111:['0.5-'],112:['2-'],113:['0.5-'],114:['1+IG'],
# 	115:['0.5-IG', '0.5-'],116:['1-','1-IG'],117:['0.5-IG', '0.5-'],118:['0-','0-IG'],119:['0.5-IG', '0.5-'],120:['0-','0-IG'],121:['0.5-IG', '0.5-'],122:['0-','3-','6+','7-','9-'],123:['0.5-']}, 
# 	param = 'Charge_radius', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4,
# 	fmt_func = P.get_shape_label, color_func = P.get_color_I, fillstyle_func = P.get_fillstyle_parity)
# secax = ax.secondary_xaxis('top', functions = (lambda x: x + 47, lambda x: x - 47))
# secax.set_xlabel('Mass number', fontsize = 25)
# secax.xaxis.set_major_locator(ticker.MultipleLocator(4))
# secax.tick_params(axis='both', which='major', labelsize=20)
# vspan = ax.plot([97-47,97-47],[-0.4,1], color = 'orange', label = 'Magic numbers', linewidth = 3)
# ax.plot([129-47,129-47],[-0.4,1], color = 'orange', linewidth = 3)
# # ax.set_ylim(-0.4,0.9)
# # ax.legend(fontsize = 16, ncol = 2, handles = [P._legend_handles[8],P._legend_handles[-1],P._legend_handles[1],P._legend_handles[2],P._legend_handles[0],P._legend_handles[5],P._legend_handles[6],P._legend_handles[7],P._legend_handles[4],P._legend_handles[3]]  + [vspan[0]], title = 'Spins', title_fontsize = 20)
# ax.legend(fontsize = 16, ncol = 2, handles = P._legend_handles, title = 'Spins', title_fontsize = 20)

# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {102:'2+lit', 104:'2+lit', 106:['1+lit', '1+'],108:['1+lit','1+'],110:['1+lit','1+'],112:['2-lit','2-'],114:['1+IG'],116:['1-','1-IG'],118:['0-','0-IG'],
# 	120:['0-','0-IG'],122:['0-']}, param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4)
# vspan = ax.plot([97-47,97-47],[-10,10], color = 'orange', label = 'Magic numbers', linewidth = 3)
# ax.plot([129-47,129-47],[-10,10], color = 'orange', linewidth = 3)
# x_in_exp = np.arange(106,111,2) - 49
# y_in_exp = np.array([uf(4.8924,0.0040),uf(4.936,0.018),uf(4.336,0.012)])/2
# in_exp = ax.errorbar(x = x_in_exp, y = unumpy.nominal_values(y_in_exp), yerr = unumpy.std_devs(y_in_exp), label = r'2$^{+}$ In CRIS',fmt = 'ks', markersize=10, markeredgecolor='k')
# ax.set_xlim([49,83])
# ax.set_ylim(-0.4,3.5)
# ax.legend(loc = 'upper right', fontsize = 16, ncol = 2, handles = [P._legend_handles[8],P._legend_handles[9],P._legend_handles[6],P._legend_handles[7],P._legend_handles[2],P._legend_handles[5],P._legend_handles[4],P._legend_handles[3],P._legend_handles[0], in_exp, vspan[0]])

# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {116:['4+','4+IG'],118:['4+','4+IG'],120:['4+','4+IG']}, param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, 
# 	n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4)


# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {102:'5+lit', 104:'5+lit', 106:['6+lit'],108:['6+lit','6+'],110:['6+lit','6+'],116:['4+', '7-','7-IG'],118:['4+','7-','7-IG','6-IG','8-IG'],
# 	120:['4+','7-','7-IG'],122:['7-']}, param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4)
# vspan = ax.plot([97-47,97-47],[-10,10], color = 'orange', label = 'Magic numbers', linewidth = 3)
# ax.plot([129-47,129-47],[-10,10], color = 'orange', linewidth = 3)
# # x_in_exp_3 = np.array([126,128]) -49
# x_in_exp_5 = np.array([114,116,118,120,122,130]) - 49
# x_in_exp_6 = np.arange(102,105,2) - 49
# x_in_exp_7 = np.array([106,108,110,112]) - 49 # 112 very short lived 0.69 microsec
# x_in_exp_8 = np.array([116,118,120,122,124,126,128,130]) - 49
# x_in_exp_10 = np.array([130]) - 49
# # y_in_exp_3 = np.array([uf(4.000,0.016),uf(4.135,0.033)])/3
# y_in_exp_5 = np.array([uf(4.6653,0.0018),uf(4.2270,0.0034),uf(4.229,0.002),uf(4.2926,0.0024),uf(4.3156,0.0028),uf(5.0138,0.0044)])/5
# y_in_exp_6 = (np.array([uf(4.5648,0.0083),uf(4.5862,0.0088)])/6)
# y_in_exp_7 = (np.array([uf(4.916,0.007),uf(4.561,0.003),uf(4.713,0.008),uf(4.73,0.04)])/7)
# y_in_exp_8 = np.array([uf(3.215,0.011),uf(3.321,0.011),uf(3.692,0.004),uf(3.781,0.006),uf(3.888,0.009),uf(4.061,0.004),uf(4.6549,0.0033),uf(5.0138,0.0044)])/8
# y_in_exp_10 = np.array(uf(5.0138,0.0044))/10
# # in_exp_3 = ax.errorbar(x = x_in_exp_3, y = unumpy.nominal_values(y_in_exp_3), yerr = unumpy.std_devs(y_in_exp_3), label = r'3$^{+}$ In Adam thesis',fmt = 'rv-', markersize=10, markeredgecolor='k')
# in_exp_5 = ax.errorbar(x = x_in_exp_5, y = unumpy.nominal_values(y_in_exp_5), yerr = unumpy.std_devs(y_in_exp_5), label = r'5$^{+}$ In Literature',fmt = 'gv-', markersize=10, markeredgecolor='k')
# in_exp_6 = ax.errorbar(x = x_in_exp_6, y = unumpy.nominal_values(y_in_exp_6), yerr = unumpy.std_devs(y_in_exp_6), label = r'6$^{+}$ In CRIS',fmt = 'kv-', markersize=10, markeredgecolor='k')
# in_exp_7 = ax.errorbar(x = x_in_exp_7, y = unumpy.nominal_values(y_in_exp_7), yerr = unumpy.std_devs(y_in_exp_7), label = r'7$^{+}$ In Literature',fmt = 'mv-', markersize=10, markeredgecolor='k')
# in_exp_8 = ax.errorbar(x = x_in_exp_8, y = unumpy.nominal_values(y_in_exp_8), yerr = unumpy.std_devs(y_in_exp_8), label = r'8$^{-}$ In CRIS',fmt = 'bv-', markersize=10, markeredgecolor='k', fillstyle = 'none')
# in_exp_10 = ax.errorbar(x = x_in_exp_10, y = unumpy.nominal_values(y_in_exp_10), yerr = unumpy.std_devs(y_in_exp_10), label = r'10$^{-}$ In Adam thesis',fmt = 'rv-', markersize=10, markeredgecolor='r', fillstyle = 'none')
# ax.set_xlim([49,83])
# ax.set_ylim(0.3,1.1)
# ax.legend(loc = 'upper right', fontsize = 12, ncol = 2, handles = [P._legend_handles[0], in_exp_5, P._legend_handles[2], in_exp_6, P._legend_handles[1], P._legend_handles[3], P._legend_handles[4],in_exp_7, in_exp_8, in_exp_10,vspan[0]], title = 'Spins', title_fontsize = 20)
# # ax.legend(loc = 'upper right', fontsize = 12, ncol = 2, handles = [P._legend_handles[0], in_exp_5, P._legend_handles[2], P._legend_handles[1], P._legend_handles[3], P._legend_handles[4],in_exp_7, in_exp_8,vspan[0]], title = 'Spins', title_fontsize = 20)

# not useful
# ax = P.plot_HF_factors(ax = ax , dict_mass_spin = {103:['3.5+lit'],105:['0.5-lit','3.5+lit'],107:['0.5-lit', '0.5-','3.5+lit', '3.5+'],109:['0.5-lit','3.5+lit', '3.5+'],
# 	111:['0.5-lit','0.5-','3.5+'],113:['0.5-lit','0.5-', '3.5+IG', '3.5+'],115:['0.5-IG','0.5-','3.5+IG', '3.5+'],117:['0.5-IG','0.5-','3.5+IG', '3.5+'],
# 	119:['0.5-IG','0.5-','3.5+IG', '3.5+'],121:['0.5-IG','0.5-','3.5+IG', '3.5+'],123:['0.5-','3.5+']}, param =  'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, 
# 	n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4)


ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {105:['0.5-lit'],107:['0.5-lit', '0.5-'],109:['0.5-lit','0.5-'],111:['0.5-lit','0.5-'],113:['0.5-lit','0.5-'],115:['0.5-IG','0.5-'],
	117:['0.5-IG','0.5-'],119:['0.5-IG','0.5-'],121:['0.5-IG','0.5-'],123:['0.5-']}, param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, 
	legend_ncols = 3, ax_major_ticker_distance = 4, fmt_func = P.get_shape_label, color_func = P.get_color_label, fillstyle_func = P.get_fillstyle_label, label_addon = 'Ag ')
vspan = ax.plot([97-47,97-47],[-10,10], color = 'orange', label = 'Magic number', linewidth = 3)
ax.plot([129-47,129-47],[-10,10], color = 'orange', linewidth = 3)
ax.plot([50-47,140-47],np.array([-0.2642833333333333,-0.2642833333333333])/0.5, '-', color = 'green', linewidth = 3)
ax.set_xlim([49,83])
# ax.set_ylim(-1.2,0.4)
x_DFT = np.arange(50,83,2)
y_DFT1 = np.array([-0.24288,-0.31569,-0.30511,-0.38870,-0.43423, -0.42638, -0.47881, -0.48407, -0.45119, -0.47987, -0.24751, -0.31334, -0.38657, -0.44348, -0.43235, -0.32783, -0.24214])/0.5
y_DFT2 = np.array([-0.26433,-0.26281,-0.24067,-0.30280,-0.33281, -0.31900, -0.34705, -0.35133, -0.32482, -0.36329, -0.26433, -0.28648, -0.34501, -0.39039, -0.38055, -0.28410, -0.26433])/0.5
DFT1 = ax.plot(x_DFT, y_DFT1, label = r'DFT HF with time-odd fields')
DFT2 = ax.plot(x_DFT, y_DFT2, '--', label = r'DFT HF without time-odd fields')
ab_initio_in_x = np.arange(50,83,2)
ab_initio_in_y_1 = np.array([-0.177,-0.227,-0.262,-0.287,-0.283,-0.344,-0.343,-0.343,-0.337,-0.441,-0.514,-0.537,-0.551,-0.558,-0.558 ,-0.538 ,-0.268])/0.5
ab_initio_in_y_2 = np.array([-0.178,-0.22,-0.264,-0.289,-0.301,-0.307,-0.301,-0.289,-0.355,-0.454,-0.515,-0.535,-0.546,-0.551,-0.546,-0.703,-0.34])/0.5
ab_initio_in_1 = ax.plot(ab_initio_in_x, ab_initio_in_y_1, '-.', label = 'VS-IMSRG 1.8/2.0(EM)')
# ab_initio_in_2 = ax.plot(ab_initio_in_x, ab_initio_in_y_2, '-',label = r'VS-IMSRG N$^{2}$LO$_{GO}$')
x_in = np.arange(103,112,2) - 49
y_in = np.array([uf(-0.1138,0.0024), uf(-0.1327,0.0027), uf(-0.1474,0.0046), uf(-0.1672,0.0028), uf(-0.185,0.033)])/0.5
x_in_lit = np.arange(113,132,2) - 49
y_in_lit = np.array([uf(-0.21,0.01), uf(-0.2405,0.0038), uf(-0.276,0.027), uf(-0.342,0.012), uf(-0.3600,0.0041), uf(-0.4047,0.0054), uf(-0.450,0.017), uf(-0.4355,0.0024), uf(-0.3871,0.0006), uf(-0.051,0.003)])/0.5
exp_in = ax.errorbar(x = x_in, y = unumpy.nominal_values(y_in), yerr = unumpy.std_devs(y_in), label = r'In CRIS$^{[1,2]}$',fmt = 'k*', markersize=15, markeredgecolor='k')
exp_in_lit = ax.errorbar(x = x_in_lit, y = unumpy.nominal_values(y_in_lit), yerr = unumpy.std_devs(y_in_lit), fmt = 'k*', markersize=15, markeredgecolor='k')
# exp_in_lit = ax.errorbar(x = x_in_lit, y = unumpy.nominal_values(y_in_lit), yerr = unumpy.std_devs(y_in_lit), label = r'In Literature',fmt = 'm*', markersize=15, markeredgecolor='k')
x_rh = np.array([103])-45
y_rh = np.array([uf(-0.08840,0.0002)])/0.5
exp_rh = ax.errorbar(x_rh,unumpy.nominal_values(y_rh), yerr = unumpy.std_devs(y_rh), label = 'Rhodium Z = 45', fmt = 'm*', markersize = 15)
x_nb = np.array([91])-41
y_nb = np.array([uf(-0.101,0.002)])/0.5
exp_nb = ax.errorbar(x_nb,unumpy.nominal_values(y_nb), yerr = unumpy.std_devs(y_nb), label = 'Niobium Z = 41', fmt = 'gs', fillstyle = 'none', markersize = 12)
x_y = np.array([89,91,93,95,97])-39
y_y = np.array([uf(-0.1374154,0.0000003),uf(-0.1641,0.0008),uf(-0.12,0.03),uf(-0.16,0.03),uf(-0.12,0.01)])/0.5
exp_y = ax.errorbar(x_y,unumpy.nominal_values(y_y), yerr = unumpy.std_devs(y_y), label = 'Yttrium Z = 39', fmt = 's', color = 'magenta', fillstyle = 'none', markersize = 12)

# ax.legend(fontsize = 14, ncol = 3, handles = P._legend_handles + [exp_in, ab_initio_in_1[0], DFT1[0], DFT2[0], vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)
# ax.legend(fontsize = 14, ncol = 3, handles = P._legend_handles + [exp_in, exp_rh, ab_initio_in_1[0], DFT1[0], DFT2[0], vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)
ax.legend(fontsize = 14, ncol = 3, handles = P._legend_handles + [exp_in, exp_rh, exp_nb, ab_initio_in_1[0], DFT1[0], DFT2[0], vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)
ax.legend(fontsize = 14, ncol = 3, handles = P._legend_handles + [exp_in, exp_rh, exp_nb, exp_y, ab_initio_in_1[0], DFT1[0], DFT2[0], vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)
# ax.legend(fontsize = 16, ncol = 2, handles = P._legend_handles + [exp_in, exp_in_lit, vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)
# ax.legend(fontsize = 16, ncol = 2, handles = P._legend_handles + [vspan[0]], loc = 'upper left', title = r'Spin $\frac{1}{2}^-$', title_fontsize = 18)

# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {103:['3.5+lit'],105:['3.5+lit'],107:['3.5+lit', '3.5+'],109:['3.5+lit', '3.5+'],111:['3.5+'],113:['3.5+IG', '3.5+'],
# 	115:['3.5+IG', '3.5+'],117:['3.5+IG', '3.5+'],119:['3.5+IG', '3.5+'],121:['3.5+IG', '3.5+'],123:['3.5+']}, param = 'g-factor', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, 
# 	legend_ncols = 3, ax_major_ticker_distance = 4, fmt_func = P.get_shape_I, color_func = P.get_color_label, label_addon = 'Ag ')
# vspan = ax.plot([82,82],[-10,10], color = 'orange', label = 'Magic numbers', linewidth = 3)
# ax.plot([50,50],[-10,10], color = 'orange', linewidth = 3)
# ax.plot([30,140],np.array([6.79285,6.79285])/4.5, '-', color = 'green', linewidth = 3)
# ax.set_xlim([49,83])
# ax.set_ylim(0.6,1.6)
# x_DFT = np.arange(99,132,2) - 49
# y_DFT = np.array([6.32224, 5.63223, 5.71883, 5.60755, 5.67050, 5.56458, 5.50175, 5.68074, 5.38345, 5.55459, 5.56822, 5.36118, 5.34454, 5.37163, 5.42497, 5.56245, 6.30967])/4.5
# DFT = ax.plot(x_DFT, y_DFT, label = r'DFT HF core coupling all I$^{\pi}$')
# x_in = np.arange(101,112,2) -49
# y_in = np.array([uf(5.8674,0.0082), uf(5.7599,0.0043), uf(5.6668,0.0050), uf(5.6082,0.0081), uf(5.5283,0.0055), uf(5.5153,0.0052)])/4.5
# x_in_lit = np.arange(113,132,2) - 49
# y_in_lit = np.array([uf(5.5264,0.0019), uf(5.5408,0.0002), uf(5.5286,0.0043), uf(5.499,0.062), uf(5.575,0.062), uf(5.442,0.061), uf(5.496,0.024), uf(5.5321,0.0014), uf(5.5961,0.0023), uf(6.312,0.014)])/4.5
# # exp_in_lit = ax.errorbar(x = x_in_lit, y = unumpy.nominal_values(y_in_lit), yerr = unumpy.std_devs(y_in_lit), label = r'$\frac{9}{2}^+$In CRIS Literature',fmt = 'ms', markersize=10, markeredgecolor='k')
# exp_in_lit = ax.errorbar(x = x_in_lit, y = unumpy.nominal_values(y_in_lit), yerr = unumpy.std_devs(y_in_lit), label = r'In CRIS$^{[1,2]}$',fmt = 'ms', markersize=10, markeredgecolor='k')
# exp_in = ax.errorbar(x = x_in, y = unumpy.nominal_values(y_in), yerr = unumpy.std_devs(y_in), fmt = 'ms', markersize=10, markeredgecolor='k')
# # exp_in = ax.errorbar(x = x_in, y = unumpy.nominal_values(y_in), yerr = unumpy.std_devs(y_in), label = r'$\frac{9}{2}^+$In CRIS',fmt = 'ks', markersize=10, markeredgecolor='k')
# x_raf = np.arange(97,103,2) - 47
# y_raf = np.array([uf(6.13,0.12), uf(5.81,0.03), uf(5.57,0.04)])/4.5
# exp_raf = ax.errorbar(x = x_raf, y = unumpy.nominal_values(y_raf), yerr = unumpy.std_devs(y_raf), label = r'$\frac{9}{2}^+$ Literature', fmt = 'b.', markersize=15)
# x_rh = np.array([99,101,103,105]) - 45
# y_rh92 = np.array([uf(5.62,0.06),uf(5.43,0.06)])/4.5
# y_rh72 = np.array([uf(4.50,0.05),uf(4.41,0.05)])/3.5
# exp_rh_lit = plt.errorbar(x_rh, unumpy.nominal_values(list(y_rh92)+list(y_rh72)), yerr = unumpy.std_devs(list(y_rh92)+list(y_rh72)), label = 'Rhodium Z = 45', fmt = 's', color = 'red', fillstyle = 'none', markersize = 10)
# x_tc = np.array([93,95,99]) - 43
# y_tc = np.array([uf(6.32,0.06),uf(5.94,0.06),uf(5.6847,0.0004)])/4.5
# exp_tc_lit = plt.errorbar(x_tc, unumpy.nominal_values(y_tc), yerr = unumpy.std_devs(y_tc), label = 'Technetium Z = 43', fmt = '*', color = 'limegreen', markersize = 13)
# x_nb = np.array([89,91,93,95,97,99,101,103])-41
# y_nb92 = np.array([uf(6.216,0.005), uf(6.521,0.002), uf(6.1705,0.0003), uf(6.141,0.005), uf(6.153,0.005), uf(5.97,0.03)])/4.5 # with NMR-on
# y_nb52 = np.array([uf(3.190,0.002), uf(3.137,0.004)])/2.5
# # exp_nb_lit = plt.errorbar(x_nb, unumpy.nominal_values(list(y_nb92)+list(y_nb52)), yerr = unumpy.std_devs(list(y_nb92)+list(y_nb52)), label = 'Niobium Z = 41', fmt = 'gs', fillstyle = 'none', markersize = 11)
# x_y = np.array([89,91,93,97,99,101])-39
# y_y92 = np.array([uf(6.37,0.04),uf(5.96,0.04),uf(6.04,0.03),uf(5.88,0.02)])/4.5
# y_y52 = np.array([uf(3.18,0.02),uf(3.22,0.02)])/2.5
# # exp_y_lit = plt.errorbar(x_y, unumpy.nominal_values(list(y_y92)+list(y_y52)), yerr = unumpy.std_devs(list(y_y92)+list(y_y52)), label = 'Yttrium Z = 39', fmt = 'k*', fillstyle = 'none', markersize = 12)
# x_rb = np.array([97]) - 37
# y_rb32 = np.array([uf(1.841,0.002)])/1.5
# # exp_rb_lit = plt.errorbar(x_rb, unumpy.nominal_values(y_rb32), yerr = unumpy.std_devs(y_rb32), label = 'Rubidium Z = 37', fmt = 'm*', markersize = 12)

# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, DFT[0], vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, exp_rh_lit, DFT[0], vspan[0]], loc = 'lower left')
# ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, exp_rh_lit, exp_tc_lit, DFT[0], vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, exp_rh_lit, exp_tc_lit, exp_nb_lit, DFT[0], vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, exp_rh_lit, exp_tc_lit, exp_nb_lit, exp_y_lit, DFT[0], vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_in_lit, exp_rh_lit, exp_tc_lit, exp_nb_lit, exp_y_lit, exp_rb_lit, DFT[0], vspan[0]], loc = 'lower left')
# # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [exp_raf, exp_in, exp_in_lit, exp_nb_lit, exp_tc_lit, exp_rh_lit, DFT[0], vspan[0]], loc = 'lower left')
# # # ax.legend(fontsize = 20, ncol = 2, handles = P._legend_handles + [ vspan[0]], loc = 'lower left')

# x50 = np.array([41,43,47])
# x52 = np.array([41,43,47,49])
# x54 = np.array([41,45,47,49])
# x56 = np.array([41,43,45,47,49])
# x58 = np.array([41,45,47,49])
# x60 = np.array([41,45,47,49])
# x62 = np.array([41,47,49])
# yn50 = np.array([uf(6.521,0.002)/4.5,uf(6.32,0.06)/4.5,uf(6.13,0.12)/4.5])
# yn52 = np.array([uf(6.1705,0.0003)/4.5,uf(5.94,0.06)/4.5,uf(5.81,0.03)/4.5,uf(5.8674,0.0082)/4.5])
# yn54 = np.array([uf(6.141,0.005)/4.5,uf(5.62,0.06)/4.5,uf(5.57,0.04)/4.5,uf(5.7599,0.0043)/4.5])
# yn56 = np.array([uf(6.153,0.005)/4.5,uf(5.6847,0.0004)/4.5,uf(5.43,0.06)/4.5,P.get_HF_moments_xy(103, '3.5+', 'g-factor', False, 47, 'Literature')[1],uf(5.6668,0.0050)/4.5])#uf(4.8,0.4)/4.5])
# yn58 = np.array([uf(5.97,0.03)/4.5,uf(4.50,0.05)/3.5,P.get_HF_moments_xy(105, '3.5+', 'g-factor', False, 47, 'Literature')[1],uf(5.6082,0.0081)/4.5])
# yn60 = np.array([uf(3.190,0.002)/2.5,uf(4.41,0.05)/3.5,P.get_HF_moments_xy(107, '3.5+', 'g-factor', False, 47, 'CRIS')[1],uf(5.5283,0.0055)/4.5])
# yn62 = np.array([uf(3.137,0.004)/2.5,P.get_HF_moments_xy(109, '3.5+', 'g-factor', False, 47, 'CRIS')[1],uf(5.5153,0.0052)/4.5])

# ax.plot(x50,unumpy.nominal_values(yn50), 'rx-', markersize = 10, label = 'N = 50')
# ax.plot(x52,unumpy.nominal_values(yn52), 'gx-', markersize = 10, label = 'N = 52')
# ax.plot(x54,unumpy.nominal_values(yn54), 'mx-', markersize = 10, label = 'N = 54')
# ax.plot(x56,unumpy.nominal_values(yn56), 'kx-', markersize = 10, label = 'N = 56')
# ax.plot(x58,unumpy.nominal_values(yn58), 'yx-', markersize = 10, label = 'N = 58')
# ax.plot(x60,unumpy.nominal_values(yn60), 'bx-', markersize = 10, label = 'N = 60')
# ax.plot(x62,unumpy.nominal_values(yn62), 'C1x-', markersize = 10, label = 'N = 62')

# ax.set_ylim(1.2,1.5)
# ax.legend(fontsize = 16, ncol = 2, title = 'Neutron number', title_fontsize = 20)
# ax.set_ylabel('g-factor', fontsize = 25)
# ax.set_xlabel('Proton number', fontsize = 25)
# ax.tick_params(axis='both', which='major', labelsize=20)

# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {107:['0.5-'],109:['0.5-'],111:['0.5-'],113:['0.5-'],115:['0.5-IG','0.5-'],117:['0.5-IG','0.5-'],
# 	119:['0.5-IG','0.5-'],121:['0.5-IG','0.5-'],123:['0.5-']}, param = 'Charge_radius', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, 
# 	legend_ncols = 3, ax_major_ticker_distance = 4, fmt_func = P.get_shape_I, color_func = P.get_color_label)
# ax.axvspan(81.9,82.1)
# fig.tight_layout()


# # not really useful
# fig,ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12,9), sharex = True)
# ax[1] = P.plot_HF_factors(ax = ax[1], dict_mass_spin = {105:['0.5-lit'],107:['0.5-lit', '0.5-'],109:['0.5-lit'],111:['0.5-lit','0.5-'],113:['0.5-lit','0.5-'],115:['0.5-IG','0.5-'],117:['0.5-IG','0.5-'],
# 	119:['0.5-IG','0.5-'],121:['0.5-IG','0.5-'],123:['0.5-']}, param = 'Mu', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, 
# 	legend_ncols = 3, ax_major_ticker_distance = 4, fmt_func = P.get_shape_I, color_func = P.get_color_label)
# ax[0] = P.plot_HF_factors(ax = ax[0], dict_mass_spin = {103:['3.5+lit'],105:['3.5+lit'],107:['3.5+lit', '3.5+'],109:['3.5+lit', '3.5+'],111:['3.5+'],113:['3.5+IG', '3.5+'],
# 	115:['3.5+IG', '3.5+'],117:['3.5+IG', '3.5+'],119:['3.5+IG', '3.5+'],121:['3.5+IG', '3.5+'],123:['3.5+']}, param = 'Mu', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, 
# 	legend_ncols = 3, ax_major_ticker_distance = 4, fmt_func = P.get_shape_I, color_func = P.get_color_label)
# ax[0].axvspan(81.9,82.1)
# ax[1].axvspan(81.9,82.1)
# t = 5.283327777777778 #* 0.84
# ax[0].axhspan(t-0.005,t+0.005, color = 'g')
# ax[0].spines.bottom.set_visible(False)
# ax[1].spines.top.set_visible(False)
# ax[0].xaxis.tick_top()
# ax[0].tick_params(labeltop=False)  
# ax[1].xaxis.tick_bottom()
# kwargs = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=12,
#               linestyle="none", color='k', mec='k', mew=1, clip_on=False)
# ax[0].plot([0, 1], [0, 0], transform=ax[0].transAxes, **kwargs)
# ax[1].plot([0, 1], [1, 1], transform=ax[1].transAxes, **kwargs)
# fig.tight_layout()

# not so useful
# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {106:['1+'],107:['0.5-','3.5+','0.5-IG'],108:['1+','6+'],109:['3.5+','0.5-IG'],110:['1+','6+'],111:['0.5-','3.5+'],112:['2-'],113:['0.5-','3.5+','3.5+IG'],114:['1+IG'],
# 	115:['0.5-IG', '0.5-','3.5+IG', '3.5+'],117:['0.5-IG', '0.5-','3.5+IG', '3.5+'],119:['0.5-IG', '0.5-','3.5+IG', '3.5+'],121:['0.5-IG', '0.5-','3.5+IG', '3.5+'],123:['0.5-','3.5+']},
# 	param = 'Charge_radius', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4, 
# 	fmt_func = P.get_shape_label, color_func = P.get_color_I)


# ax = P.plot_HF_factors(ax, dict_mass_spin = {108:['6+'],110:['6+'],116:['7-'],118:['7-'],120:['7-'],122:['7-']}, param = 'Q', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4, 
# 	fmt_func = P.get_shape_label, color_func = P.get_color_I)


# ax = P.plot_HF_factors(ax = ax, dict_mass_spin = {107:['3.5+'],109:['3.5+'],111:['3.5+'],113:['3.5+','3.5+IG'],115:['3.5+','3.5+IG'],117:['3.5+','3.5+IG'],
# 	119:['3.5+','3.5+IG'],121:['3.5+','3.5+IG'],123:['3.5+']},
# 	param = 'Q', getter_xy = P.get_HF_moments_xy, scale_I = False, n_proton = 47, legend_ncols = 3, ax_major_ticker_distance = 4, 
# 	fmt_func = P.get_shape_I, color_func = P.get_color_label)
# secax = ax.secondary_xaxis('top', functions = (lambda x: x + 47, lambda x: x - 47))
# secax.set_xlabel('Mass number', fontsize = 25)
# secax.xaxis.set_major_locator(ticker.MultipleLocator(4))
# secax.tick_params(axis='both', which='major', labelsize=20)


plt.show()