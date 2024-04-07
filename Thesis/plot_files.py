import os
import sys
import math
import shutil
import numpy as np
import pylab as P
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def folder_create(name,force):
	if (not os.path.exists(name)):
		os.mkdir(name)
	else:
		if (force):
			shutil.rmtree(name)
			os.mkdir(name)
def plot_files(proc, ReWriteQ):

	if (proc == "muonP"):
		process = r"$e^{-} \mu^{+} \to e^{-} \mu^{+}$"
	elif (proc == "muonM"):
		process = r"$e^{-} \mu^{-} \to e^{-} \mu^{-}$"

	process = process+ " | " +r"Soft photon emission"

# path to data (correct to fit your own) :
	path_data='/home/ldninh/hahn/working/eMu_eMu/driver/data_dist/'

# result folder
	folder_target = "./plots"

	if (not ReWriteQ): 
		if (os.path.exists(folder_target)): 
			return

# true = remove the folder and create an empty one
	folder_create(folder_target,False)

# input a list of distributions here:
	list_dist = ["cos_CM_plot","t_inv_plot"]
	list_label_x = [r"$\cos(\theta)$",r"$t$ [$GeV^2$]"]
	list_output = [r"cos_CM_plot",r"t_inv_plot"]
	list_label_y = [r"$d\sigma/d\cos(\theta)$ [$\mu$b]",r"$d\sigma/dt$ [$\mu$b/$GeV^2$]"]
#	y_plot_scale = ["linear","linear"]
	y_plot_scale = ["log","log"]
	Xsection_unit = r"$\mu$b"

	pi = math.pi
	alpha=1/137.035999074306
	e=np.sqrt(alpha*4*pi)
	fac = e**2/(16*pi**2)

	i = -1
	for i_dist in list_dist:
		i = i + 1

		treeQED = P.loadtxt(path_data + "bornQED/" + i_dist)
		virtQED_FC = P.loadtxt(path_data + "virtQED/" + i_dist)
		virtQED_TR = P.loadtxt(path_data + "virtQED_Truyen/" + i_dist)
		virtQED_TR1 = P.loadtxt(path_data + "virtQED_Truyen_soft1E4/" + i_dist)

		print("Name of distribution =",i_dist)
#		print("dist =",tree)

		x = treeQED[:,0]
		dx = x[1] - x[0]
		y_treeQED = treeQED[:,1]
		y_virtQED_FC = virtQED_FC[:,1]
		y_virtQED_TR = virtQED_TR[:,1]
		y_virtQED_TR1 = virtQED_TR1[:,1]

		y_NLOsoftQED_FC = y_treeQED + y_virtQED_FC
		y_NLOsoftQED_TR = y_treeQED + y_virtQED_TR
		y_NLOsoftQED_TR1 = y_treeQED + y_virtQED_TR1

		nbin = len(x)
# for later plotting (remove the edge bins because of numerical fluctuations)
		bin_min = 2
		bin_max = nbin-3
#		if i_dist=="cos_CM_plot":
#			bin_min = 2

### Cross sections:
#		Xsection = np.sum(y_tree)*dx
#		Xsection_mad = np.sum(y_mad)*dx
#		print("Cross section Ninh =",Xsection,Xsection_unit)
#		print("Cross section Mad  =",Xsection_mad,Xsection_unit)

#		y_min = min(y_treeQED)
#		y_max = max(y_treeQED)

# two subplots
		figL = plt.figure(1)
		gridspec.GridSpec(3,3)
# large subplot, first grid (0,0) at top left corner
		ax1=plt.subplot2grid((3,3), (0,0), rowspan=2, colspan=3)

		if (y_plot_scale[i] == "log"):
			plt.semilogy(x, y_treeQED, color = 'black', label = r'LO', drawstyle = 'default', linestyle = ':')
			plt.semilogy(x, y_NLOsoftQED_FC, color = 'blue', label = r'NLO FormCalc ($\Delta E = 10^{-3}\sqrt{s}/2$)', drawstyle = 'default')
			plt.semilogy(x, y_NLOsoftQED_TR, color = 'red', label = r'NLO This work ($\Delta E = 10^{-3}\sqrt{s}/2$)', drawstyle = 'default', linestyle = '--')
			plt.semilogy(x, y_NLOsoftQED_TR1, color = 'darkviolet', label = r'NLO This work ($\Delta E = 10^{-4}\sqrt{s}/2$)', drawstyle = 'default', linestyle = '-.')
		else:
			plt.plot(x, y_treeQED, color = 'black', label = r'LO', linestyle = ':')
			plt.plot(x, y_NLOsoftQED_FC, color = 'blue', label = r'NLO FormCalc ($\Delta E = 10^{-3}\sqrt{s}/2$)')
			plt.plot(x, y_NLOsoftQED_TR, color = 'red', label = r'NLO This work ($\Delta E = 10^{-3}\sqrt{s}/2$)', linestyle = '--')
			plt.plot(x, y_NLOsoftQED_TR1, color = 'darkviolet', label = r'NLO This work ($\Delta E = 10^{-4}\sqrt{s}/2$)', linestyle = '-.')

		plt.title(process)
		plt.ylabel(list_label_y[i])
		plt.legend(loc='upper center')

#		y_step = (y_max - y_min)/100
#		y_min = y_min - y_step*2
#		y_max = y_max + y_step*2
#		plt.ylim(y_min,y_max)
# small subplot
		fac_diff = 100
		diff_L = [0 for j in range(0,nbin)]
		diff_L1 = [0 for j in range(0,nbin)]
		for j in range(0,nbin):
			if y_treeQED[j]!=0:
				diff_L[j]=(y_NLOsoftQED_TR[j]-y_treeQED[j])/y_treeQED[j]*fac_diff
				diff_L1[j]=(y_NLOsoftQED_TR1[j]-y_treeQED[j])/y_treeQED[j]*fac_diff
#				diff_L[j]=(y_NLOsoftQED_TR[j]-y_NLOsoftQED_FC[j])/(y_treeQED[j]*fac)
			else:
				diff_L[j]=0.
				diff_L1[j]=0.
		ax2=plt.subplot2grid((3,3), (2,0), rowspan=1, colspan=3, sharex=ax1)
		plt.plot(x, diff_L,  color = 'red', linestyle = '--', label = r'$\Delta E = 10^{-3}\sqrt{s}/2$')
		plt.plot(x, diff_L1, color = 'darkviolet', linestyle = '-.', label = r'$\Delta E = 10^{-4}\sqrt{s}/2$')
		line_zero = [0 for j in range(0,nbin)]
#		plt.plot(x, line_zero, 'g')
		plt.xlabel(list_label_x[i])
		if fac_diff == 100:
			plt.ylabel(r"$\delta$ [%]")
		else:
			plt.ylabel(r"$\delta$")
###
# yrange limits for Diff plots:
		y_diff_max = -0.20*fac_diff
		y_diff_min = -0.50*fac_diff
#		ymax = max(diff_L)
#		ymin = min(diff_L)
#		ymax = min([ymax,y_diff_max])
#		ymin = max([ymin,y_diff_min])
###
		axes = plt.gca()
#		axes.set_ylim([ymin,ymax])

		axes.set_ylim([y_diff_min,y_diff_max])
		axes.set_xlim([x[bin_min],x[bin_max]])
		plt.legend()
# fit subplots and save figL
		figL.tight_layout()
		figL.subplots_adjust(hspace=0.1)

		for ax in [ax1,ax2]:
			if ax == ax1:
				plt.setp(ax.get_xticklabels(), visible=False)
			ax.tick_params(axis='both',which='both',top='on',right='on',direction='in')
			ax.tick_params(axis='both',which='major',length=6.0)
			ax.tick_params(axis='both',which='minor',length=4.0)

		figL.set_size_inches(w=7,h=7)
		suffix = "_NLOsoftQED_" + y_plot_scale[i]
                tmp = folder_target + "/" + list_output[i] + suffix + ".pdf"
		figL.savefig(tmp)
		plt.close(figL)

	print("all results are in " + folder_target)
### do some actions
#def plot_files(proc_in, ReWriteQ):
# ReWriteQ: True (force rewrite), False (no rewrite)
for proc in ["muonM"]:
 plot_files(proc,True)
