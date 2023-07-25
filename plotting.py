import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import numpy as np
import pandas as pd
from moments import NuclearState as NS
from IGISOL_moments import IGISOL_NuclearState as IG_NS
from uncertainties import unumpy
from uncertainties import ufloat as uf
from typing import Callable, Union
from numpy.typing import ArrayLike

class Plotting:

	_dict_spin_style = {'0+':['blue','.'], '0-':['blue','*'], '0.5+':['black','.'], '0.5-':['black','*'], '1+':['red','.'], '1-':['red','*'], 
						'1.5+':['blue','.'], '1.5-':['blue','*'], '2+':['orange','.'], '2-':['orange','*'], '2.5+':['green','.'], '2.5-':['green','*'], 
						'3+':['orange','.'], '3-':['orange','*'], '3.5+':['black','.'], '3.5-':['orange','*'], '4+':['purple','.'], '4-':['purple','*'], 
						'4.5+':['purple','.'], '4.5-':['purple','x'], '5+':['green','.'], '5-':['cyan','*'], '5.5+':['cyan','.'], '5.5-':['cyan','*'],
						'6+':['black','.'], '6-':['black','*'], '6.5+':['black','.'], '6.5-':['black','*'], '7+':['magenta','.'], '7-':['magenta','*'],
						'7.5+':['magenta','.'], '7.5-':['magenta','*']}
	_markersize = 15
	_axtitle_fontsize = 25
	_legend_fontsize = 16
	_legend_title_fontsize = 16
	_legend_entries_done = list()
	_legend_handles = list()
	_param_to_extract = ['Al', 'Au', 'Bl', 'Bu', 'Centroid']
	_moment_to_extract = ['Mu', 'Q', 'Charge_radius', 'g-factor']
	_lit_moment = {'101_4.5':{'Mu':uf(5.627,0.011)}, '102_2':{'Mu':uf(4.1,0.3)}, '102_5':{'Mu':uf(4.6,0.7)}, '103_3.5':{'Mu':uf(4.432,0.002)}, 
				   '104_2':{'Mu':uf(3.691,0.003)}, '104_5':{'Mu':uf(3.916,0.008)}, '105_0.5':{'Mu':uf(-0.1014,0.0010)}, '105_3.5':{'Mu':uf(4.414,0.013)}, 
				   '106_1':{'Mu':uf(2.8,0.2)}, '106_6':{'Mu':uf(3.705,0.004)}, '107_0.5':{'Mu':uf(-0.11367965,0.00000015)}, '107_3.5':{'Mu':uf(4.398,0.005)},
				   '108_1':{'Mu':uf(2.6884,0.0007)}, '108_6':{'Mu':uf(3.580,0.020)}, '109_0.5':{'Mu':uf(-0.1306906,0.0000002)}, '109_3.5':{'Mu':uf(4.400,0.006)},
				   '110_1':{'Mu':uf(2.7271,0.0008)}, '110_6':{'Mu':uf(3.588,0.003)}, '111_0.5':{'Mu':uf(-0.146,0.002)}, '112_2':{'Mu':uf(-0.0547,0.0005)},
				   '113_0.5':{'Mu':uf(-0.159,0.002)}}

	def check_param_valid(self, getter_xy: Callable[[str,float,str,bool,int,str], tuple], param: str) -> bool:
		'''Check wether parameter names are valid

        Parameters
        ----------
        getter_xy: str
            name of getter function eg. whether you choose HF parameters or moments
        param: str
            name of parameter

        Returns
        -------
        bool
        '''
		if getter_xy == self.get_HF_moments_xy:
			if param in self._moment_to_extract:
				return True
			return False
		elif getter_xy == self.get_HF_parameters_xy:
			if param in self._param_to_extract:
				return True
			return False
		print('wtf did you put in getter_xy xddddd')
		return False

	def plot_HF_factors(self, ax: plt.Axes, dict_mass_spin: dict, param: str, getter_xy: Callable[[str,float,str,bool,int,str], tuple], scale_I: bool = False, n_proton: int = 0, **plot_kwargs) -> plt.Axes:
		'''reads the fit parameters for a single scan file and returns a dictionary with keys the parameter names and parameter names + _err

        Parameters
        ----------
        ax: plt.Axes
            axes object to plot the data in
        dict_mass_spin: int
            dictionary with keys mass and values a list of string with the spin, parity and facility/literature
        param: str
            parameter to plot
        getter_xy: Callable[[str,float,str,bool,int,str], tuple]
        	fuction which gets the values for the param argument
        scale_I: bool, default: False
        	whether to scale the parameter with I or not
        n_proton: int, default: 0
        	the amount of protons to subtract from the x-axis. eg to go from mass number to neutron number

        Returns
        -------
        plt.Axes
        '''
		if not self.check_param_valid(getter_xy ,param):
			raise RuntimeError('use param = Al, Au, Bl, Bu or Centroid and getter_xy = get_HF_parameters_xy, or param = Mu, Q, Charge_radius or g-factor and getter_xy = get_HF_moments_xy')
		if param == 'Charge_radius':
			scale_I = False
		for mass in dict_mass_spin.keys():
			if type(dict_mass_spin[mass]) == type([1,2,3]): # if it is a list
				for I in dict_mass_spin[mass]:
					if 'lit' in I:
						I = I[:-3]
						x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'Literature')
						ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'Literature', **plot_kwargs)
					elif 'IG' in I:
						I = I[:-2] 
						x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'IGISOL')
						ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'IGISOL', **plot_kwargs)
					else:
						x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'CRIS')
						ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'CRIS', **plot_kwargs)
			else:
				I = dict_mass_spin[mass]
				if 'lit' in I:
					I = I[:-3]
					x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'Literature')
					ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'Literature', **plot_kwargs)
				elif 'IG' in I:
					I = I[:-2]
					x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'IGISOL')
					ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'IGISOL', **plot_kwargs)
				else:
					x,y = getter_xy(mass = mass, I = I, param = param, scale_I = scale_I, n_proton = n_proton, string_nuclear_state_class = 'CRIS')
					ax = self.ax_parameters(ax, x_value = x, y_value = y, I = I, label = 'CRIS', **plot_kwargs)
		return self.plot_layout(ax = ax, param = param, n_proton = n_proton, legend_handles = self._legend_handles, **plot_kwargs)

	def get_HF_parameters_xy(self, mass: int, I: float, param: str, scale_I: bool, n_proton: int, string_nuclear_state_class: str) -> tuple:
		'''getter function for HF parameters, returning a tuple with the x eg mass or neutron number and y the values + errors of the selected parameter

        Parameters
        ----------
        mass: int
            mass number of the isotope
        I: float
            spin of the state
        param: str
            parameter to plot
        scale_I: bool
        	whether to scale the parameter with I or not
        n_proton: int
        	the amount of protons to subtract from the x-axis. eg to go from mass number to neutron number
        string_nuclear_state_class: str
			the facility at which the measurement was taken

        Returns
        -------
        tuple
        '''
		if string_nuclear_state_class == 'IGISOL':
			try:
				y = IG_NS(mass, float(I[:-1])).HFparam_dict[param] 
			except:
				y = IG_NS(mass, int(I[:-1])).HFparam_dict[param]
			x = self.adapt_x_nproton(mass, n_proton)
			if scale_I:
				y = self.scale_y_to_spin(y, I[:-1])
			return x,y
		elif string_nuclear_state_class == 'CRIS':
			try:
				y = NS(mass, float(I[:-1])).HFparam_dict[param] 
			except:
				y = NS(mass, int(I[:-1])).HFparam_dict[param]
			x = self.adapt_x_nproton(mass, n_proton)
			if scale_I:
				y = self.scale_y_to_spin(y, I[:-1])
			return x,y
		else:
			raise RuntimeError('give CRIS or IGISOL for string_nuclear_state_class')

	def get_HF_moments_xy(self, mass: int, I: float, param: str, scale_I: bool, n_proton: int, string_nuclear_state_class: str) -> tuple:
		'''getter function for HF moments, returning a tuple with the x eg mass or neutron number and y the values + errors of the selected parameter

        Parameters
        ----------
        mass: int
            mass number of the isotope
        I: float
            spin of the state
        param: str
            parameter to plot
        scale_I: bool
        	whether to scale the parameter with I or not
        n_proton: int
        	the amount of protons to subtract from the x-axis. eg to go from mass number to neutron number
        string_nuclear_state_class: str
			the facility at which the measurement was taken

        Returns
        -------
        tuple
        '''
		denominator = 1
		factor = 1
		if param == 'g-factor':
			param = 'Mu'
			if I[:-1] not in ['0','1']:
				denominator = float(I[:-1])
		if param == 'Q':
			factor = 100
		if string_nuclear_state_class == 'IGISOL':
			try:
				y = IG_NS(mass, float(I[:-1])).Nuc_Prop_dict[param] 
			except:
				y = IG_NS(mass, int(I[:-1])).Nuc_Prop_dict[param]
			x = self.adapt_x_nproton(mass, n_proton)
			# if scale_I:
			# 	y = self.scale_y_to_spin(y, I[:-1])
		elif string_nuclear_state_class == 'CRIS':
			try:
				y = NS(mass, float(I[:-1])).Nuc_Prop_dict[param] 
			except:
				y = NS(mass, int(I[:-1])).Nuc_Prop_dict[param]
			x = self.adapt_x_nproton(mass, n_proton)
			# if scale_I:
			# 	y = self.scale_y_to_spin(y, I[:-1])
		elif string_nuclear_state_class == 'Literature':
			x = self.adapt_x_nproton(mass, n_proton)
			y = self._lit_moment[str(mass) + '_' + str(I[:-1])][param]
		else:
			raise RuntimeError('give CRIS or IGISOL for string_nuclear_state_class')
		if scale_I:
			y = self.scale_y_to_spin(y, I[:-1])
		return x, factor*(y/denominator)

	def ax_parameters(self, ax: plt.Axes, x_value: Union[uf,unumpy.uarray], y_value: Union[uf,unumpy.uarray], I: float, label: str, **plot_kwargs) -> plt.Axes:
		'''Plots the given x and y value with all miscellaneous plot features. Also svaes legend handles in a list

        Parameters
        ----------
        ax: plt.Axes
            axes object to plot the data in
        x_value: unumpy.uarray or ufloat
            value of x
        y_value: Unumpy.uarray or ufloat
            value of y
        I: float
        	spin of the state
        label: str
        	label indicating from which facility/literature

        Returns
        -------
        plt.Axes
        '''
		if I + label not in self._legend_entries_done:
			self._legend_handles.append(ax.errorbar(x = x_value, y = unumpy.nominal_values(y_value), yerr = unumpy.std_devs(y_value), 
				markersize = self._markersize, fmt = plot_kwargs.get('fmt_func', self.get_shape_label)(label_I = [label,I[:-1]]), 
				color = plot_kwargs.get('color_func', self.get_color_I)(label_I = [label,I[:-1]]), 
				fillstyle = plot_kwargs.get('fillstyle_func', self.get_fillstyle_parity)(label_Ip = [label,I[:-1],I[-1]]),
				label = plot_kwargs.get('label_addon', self.convert_decimalstring_fractionstring(float(I[:-1])) + r'$^{' + I[-1] + '}$ ') + label))
		else:
			ax.errorbar(x = x_value, y = unumpy.nominal_values(y_value), yerr = unumpy.std_devs(y_value), 
				markersize = self._markersize, fmt = plot_kwargs.get('fmt_func', self.get_shape_label)(label_I = [label,I[:-1]]), 
				color = plot_kwargs.get('color_func', self.get_color_I)(label_I = [label,I[:-1]]), 
				fillstyle = plot_kwargs.get('fillstyle_func', self.get_fillstyle_parity)(label_Ip = [label,I[:-1],I[-1]]),
				label = plot_kwargs.get('label_addon', self.convert_decimalstring_fractionstring(float(I[:-1])) + r'$^{' + I[-1] + '}$ ') + label)
		self._legend_entries_done.append(I + label)
		return ax

	def get_fillstyle_parity(self, label_Ip: list) -> str:
		'''return fillstyle depending on the parity

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		parity = label_Ip[2]
		if parity == '+':
			return 'full'
		if parity == '-':
			return 'none'
		return 'none'

	def get_fillstyle_label(self, label_Ip: list) -> str:
		'''return fillstyle depending on the label

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		label = label_Ip[0]
		if label == 'CRIS':
			return 'full'
		if label == 'IGISOL':
			return 'none'
		if label == 'Literature':
			return 'bottom'

	def get_shape_label(self, label_I: list) -> str:
		'''return shape depending on the label

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		label = label_I[0]
		if label == 'CRIS':
			return '.'
		if label == 'IGISOL':
			return 'X'
		if label == 'Literature':
			return '*'

	def get_shape_I(self, label_I: list) -> str:
		'''return shape depending on the spin

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		I = label_I[1]
		if I == '0':
			return 's'
		if I == '0.5':
			return '*'
		if I == '1':
			return 'X'
		if I == '2':
			return '.'
		if I == '3':
			return 'X'
		if I == '3.5':
			return '.'
		if I == '4':
			return '*'
		if I == '5':
			return '.'
		if I == '6':
			return 'X'
		if I == '7':
			return '*'
		if I == '8':
			return 's'

	def get_color_label(self, label_I: list) -> str:
		'''return color depending on the label

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		label = label_I[0]
		if label == 'CRIS':
			return 'black'
		if label == 'IGISOL':
			return 'orange'
		if label == 'Literature':
			return 'blue'

	def get_color_I(self, label_I: list) -> str:
		'''return color depending on the spin

		Parameters
		----------
		label_Ip: list
			list with label, spin and parity

		Returns
		-------
		str
		'''
		I = label_I[1]
		if I == '0':
			return 'blue'
		if I == '0.5':
			return 'black'
		if I == '1':
			return 'red'
		if I == '2':
			return 'orange'
		if I == '3':
			return 'green'
		if I == '3.5':
			return 'red'
		if I == '4':
			return 'black'
		if I == '5':
			return 'red'
		if I == '6':
			return 'black'
		if I == '7':
			return 'magenta'
		if I == '8':
			return 'blue'

	def plot_layout(self, ax: plt.Axes, param: str, n_proton: int, legend_handles: list, **plot_kwargs) -> plt.Axes:
		'''add all miscellaneous features to plot such as xlabel, ylabel, etc.

		Parameters
        ----------
        ax: plt.Axes
            axes object to plot the data in
        param: str
            parameter to plot
        n_proton: int
        	the amount of protons to subtract from the x-axis. eg to go from mass number to neutron number
        legend_handles: list
        	all the legend handles to add to the plot

		Returns
		-------
		plt.Axes
		'''
		# handles, labels = ax.get_legend_handles_labels()
		# ar = np.array([labels,handles])
		# new_dict = dict()
		# for i in range(0,len(ar[0])):
		# 	if labels[i] not in new_dict.keys():
		# 		new_dict[labels[i]] = handles[i]
		# labels, handles = list(new_dict.keys())[::-1],list(new_dict.values())[::-1]
		# handles = handles, labels = labels, put this in the ax.legend(...) to 'sort' it. But this fkin sucks so better dont xd
		# you are better of by asking the class for all legend_handles and then sort them manually however the fk u want it
		# smt like this: ax.legend(handles = [P._legend_handles[0],P._legend_handles[2]]]), note tho that you reset all the **plot_kwargs for the legend so...
		ax.legend(fontsize = plot_kwargs.get('legend_fontsize', self._legend_fontsize), ncol = plot_kwargs.get('legend_ncols', 1), handles = legend_handles, 
			title = plot_kwargs.get('legend_title', None), title_fontsize = plot_kwargs.get('legend_title_fontsize', self._legend_title_fontsize))
		if n_proton == 0:
			ax.set_xlabel('Mass number', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		else:
			ax.set_xlabel('Neutron number', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		if param == 'Mu':
			ax.set_ylabel(r'$\mu$ [$\mu_{N}$]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'g-factor':
			ax.set_ylabel(r'g-factor [$\mu$/I]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Q':
			ax.set_ylabel(r'Q [fm$^{2}$]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Charge_radius':
			ax.set_ylabel(r'$\delta$<r$^2$>$^{109,A}$ [fm$^2$]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Al':
			ax.set_ylabel(r'A$_{lower}$ [MHz]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Au':
			ax.set_ylabel(r'A$_{upper}$ [MHz]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Bl':
			ax.set_ylabel(r'B$_{lower}$ [MHz]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Bu':
			ax.set_ylabel(r'B$_{upper}$ [MHz]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		elif param == 'Centroid':
			ax.set_ylabel(r'Centroid [MHz]', fontsize = plot_kwargs.get('ax_fontsize', self._axtitle_fontsize))
		ax.tick_params(axis='both', which='major', labelsize=plot_kwargs.get('ticker_fontsize', 20))
		ax.xaxis.set_major_locator(ticker.MultipleLocator(plot_kwargs.get('ax_major_ticker_distance', 4)))
		if plot_kwargs.get('ax_grid', True):
			ax.grid()
		return ax

	def adapt_x_nproton(self, mass: int, n_proton: int) -> int:
		'''Change the x to go from mass number to neutron number

		Parameters
        ----------
        mass: int
            mass number of isotope
        n_proton: int
            number of protons of isotope

		Returns
		-------
		int
		'''
		if n_proton != 0:
			return mass - n_proton
		return mass

	def scale_y_to_spin(self, y: ArrayLike, I: float) -> float:
		'''Scale y with the inverse of the spin

		Parameters
        ----------
        y: ArrayLike
            y values to plot
        I: float
            spin of the nuclear state

		Returns
		-------
		float
		'''
		I = float(I)
		if I == 0:
			denominator = 1
		else:
			denominator = I
		return y/denominator

	def convert_decimalstring_fractionstring(self, stringfloat: str) -> str:
		'''converts decimal strings to fractions strings

		Parameters
        ----------
        stringfloat: str
            decimal string

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
				stringfraction = str(int(stringfloat))
			else:
				temptop = str(stringfloat.as_integer_ratio()[0])
				tempbottom = str(stringfloat.as_integer_ratio()[1])
				stringfraction = r'$\frac{' + temptop + '}{' + tempbottom + '}$'
		return stringfraction