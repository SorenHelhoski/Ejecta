import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin

dirname = 'Panels'
psp.mkdir_p(dirname)

landing_bin = 50 # number of bins of the landing position of ejecta
pres_bin = 40 # number of bins in the pressure histogram
temp_bin = 60 # number of bins in the temperature histogram

pres_range = [0,200] # range plotted GPa
temp_range = [0,4000] # range plotted K

y_lim_pres = 250 
y_lim_temp = 250 

print('Opening data file...')
# open tracer file
file2 = open('data.txt','r')
x_list = list(eval(file2.readline().replace('\n','')))
y_list = list(eval(file2.readline().replace('\n','')))
v_list = list(eval(file2.readline().replace('\n','')))
b_list = list(eval(file2.readline().replace('\n','')))
p_list = list(eval(file2.readline().replace('\n','')))
T_list = list(eval(file2.readline().replace('\n','')))
m_list = list(eval(file2.readline().replace('\n','')))
t_list = list(eval(file2.readline().replace('\n','')))
r_list = list(eval(file2.readline().replace('\n','')))
used =   list(eval(file2.readline().replace('\n','')))
jumps, start_time, end_time, save_step, R, grid_spacing, g, a, v_0, layers  = eval(file2.readline().replace('\n',''))
print('DONE\n')

landing_bounds = [.9*R,3.1*R] 
sep_size = (landing_bounds[1] - landing_bounds[0])/landing_bin
Norm = grid_spacing**2/(sep_size) # used to change tracer freq to height

Loc = [R, 6*R/5, 7*R/5, 8*R/5, 9*R/5, 10*R/5, 11*R/5, 12*R/5, 13*R/5, 14*R/5]

pres_mid = (pres_range[1]+pres_range[0])/2
temp_mid = (temp_range[1]+temp_range[0])/2

#------------------------------------------------------------
#                     Histogram
#------------------------------------------------------------

print('Seperating data for 10 panel graphs...')

#pressures in each range
p1,p2,p3,p4,p5,p6,p7,p8,p9,p0 = [],[],[],[],[],[],[],[],[],[]

#temp in each range
T1,T2,T3,T4,T5,T6,T7,T8,T9,T0 = [],[],[],[],[],[],[],[],[],[]

for i in range(len(r_list)):
    if Loc[0] < r_list[i] < Loc[1]:
        p1.append(p_list[i])
        T1.append(T_list[i])
    if Loc[1] < r_list[i] < Loc[2]:
        p2.append(p_list[i])
        T2.append(T_list[i])
    if Loc[2] < r_list[i] < Loc[3]:
        p3.append(p_list[i])
        T3.append(T_list[i])
    if Loc[3] < r_list[i] < Loc[4]:
        p4.append(p_list[i])
        T4.append(T_list[i])
    if Loc[4] < r_list[i] < Loc[5]:
        p5.append(p_list[i])
        T5.append(T_list[i])
    if Loc[5] < r_list[i] < Loc[6]:
        p6.append(p_list[i])
        T6.append(T_list[i])
    if Loc[6] < r_list[i] < Loc[7]:
        p7.append(p_list[i])
        T7.append(T_list[i])
    if Loc[7] < r_list[i] < Loc[8]:
        p8.append(p_list[i])
        T8.append(T_list[i])
    if Loc[8] < r_list[i] < Loc[9]:
        p9.append(p_list[i])
        T9.append(T_list[i])
    if Loc[9] < r_list[i] < Loc[9]+R/5:
        p0.append(p_list[i])
        T0.append(T_list[i])

p_all = p1,p2,p3,p4,p5,p6,p7,p8,p9,p0
T_all = T1,T2,T3,T4,T5,T6,T7,T8,T9,T0

print('DONE\n')

#------------------------------------------------------------
#                          Density 
#------------------------------------------------------------

print('Creating KDEs and Curve Fits...')

from aux_library import Density
from scipy.optimize import curve_fit as fit
from math import sqrt

def model_func(x,A,B,C):
    return float(A)* np.exp((-0.5)*((x-float(B))/(float(C)))**2 )

KDE_P = [[None,None,None,None]]

NKDE_P = [[None,None, None,None]]

#guess = [[530, 6, 10],[360, 6, 10],[260, 8, 11],[190, 10, 12],[150, 11, 14],[100, 14, 14],[70, 15, 12],[50, 15, 11],[40, 15, 11],[30, 15, 6]]

param = []

counter = 0
for each in p_all:
    pdf = Density(each, bounds = pres_range)
    guess = [max(pdf.get_y()), pres_mid, (pres_range[1]-pres_range[0])/4]
    popt, pcov = fit(model_func, pdf.get_x(), pdf.get_y(), p0 = guess, bounds = [[1,pres_range[0],0],[2*max(pdf.get_y()),pres_range[1],pres_range[1]]], method = 'trf')
    A_fit, B_fit, C_fit = popt
    param.append(popt)
    gaussian = []

    for i in pdf.get_x():
        gaussian.append(model_func(i, A_fit, B_fit, C_fit))

    KDE_P.append([pdf.get_x(), pdf.get_y(), gaussian])
    NKDE_P.append([pdf.get_x(), pdf.get_norm_y(), None, None])
    counter += 1

KDE_T = [[None,None]]

NKDE_T = [[None,None]]

for each in T_all:
    temperature = Bin(each, bins = temp_bin, bounds = temp_range)
    pdf = Density(each, bounds = temp_range)
    KDE_T.append([pdf.get_x(), pdf.get_y()])
    NKDE_T.append([pdf.get_x(), pdf.get_norm_y()])

print('DONE\n')


#--------------------------------------------------------------
#                       10 Panel Graphs
#--------------------------------------------------------------

#Peak Pressure over each range
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = pres_range   
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.hist(p_all[each], bins = pres_bin, range = pres_range)
    axis.plot(KDE_P[each+1][0], KDE_P[each+1][1])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    axis.set_xlim(lower,upper)
    axis.set_ylim(0,y_lim_pres)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Pressure [GPa]')
ax7.set_xlabel('Peak Pressure [GPa]')
ax8.set_xlabel('Peak Pressure [GPa]')
ax9.set_xlabel('Peak Pressure [GPa]')
ax10.set_xlabel('Peak Pressure [GPa]')

fig.savefig('{}/Peak Pressure(bin={}).png'.format(dirname,pres_bin))
print('Saved: Peak Pressure(bin={}).png'.format(pres_bin))

#Peak Pressure over each range
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = pres_range   
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.hist(p_all[each], bins = pres_bin, range = pres_range)
    axis.plot(KDE_P[each+1][0], KDE_P[each+1][2])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    aa,bb,cc = param[each]
    axis.text(pres_mid,.7*y_lim_pres, "ampd : {0:.0f} \n mean: {1:.2f} \n stdv: {2:.2f}".format(aa,bb,cc))
    axis.set_xlim(lower,upper)
    axis.set_ylim(0,y_lim_pres)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Pressure [GPa]')
ax7.set_xlabel('Peak Pressure [GPa]')
ax8.set_xlabel('Peak Pressure [GPa]')
ax9.set_xlabel('Peak Pressure [GPa]')
ax10.set_xlabel('Peak Pressure [GPa]')

fig.savefig('{}/Peak Pressure (Gaussian) (bin={}).png'.format(dirname,pres_bin))
print('Saved: Peak Pressure (Gaussian) (bin={}).png'.format(pres_bin))

#Peak Pressure over each range (pdf gaussian)
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = pres_range   
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.plot(KDE_P[each+1][0], KDE_P[each+1][1])
    axis.plot(KDE_P[each+1][0], KDE_P[each+1][2])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    aa,bb,cc = param[each]
    axis.text(pres_mid,.7*y_lim_pres, "ampd : {0:.0f} \n mean: {1:.2f} \n stdv: {2:.2f}".format(aa,bb,cc))
    axis.set_xlim(lower,upper)
    axis.set_ylim(0,y_lim_pres)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Pressure [GPa]')
ax7.set_xlabel('Peak Pressure [GPa]')
ax8.set_xlabel('Peak Pressure [GPa]')
ax9.set_xlabel('Peak Pressure [GPa]')
ax10.set_xlabel('Peak Pressure [GPa]')

fig.savefig('{}/Peak Pressure (Gaussian and KDE).png'.format(dirname))
print('Saved: Peak Pressure (Gaussian and KDE).png')

#Peak Pressure over each range
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = pres_range   
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.plot(NKDE_P[each+1][0], NKDE_P[each+1][1])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    axis.set_xlim(lower,upper)
    top = 10/float(pres_range[1]-pres_range[0])
    axis.set_ylim(0,top)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Pressure [GPa]')
ax7.set_xlabel('Peak Pressure [GPa]')
ax8.set_xlabel('Peak Pressure [GPa]')
ax9.set_xlabel('Peak Pressure [GPa]')
ax10.set_xlabel('Peak Pressure [GPa]')

fig.savefig('{}/Normed Pressure(bin={}).png'.format(dirname,pres_bin))
print('Saved: Normed Pressure(bin={}).png'.format(pres_bin))

#-------------------------------------------------

#Peak Temperature over each range
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = temp_range  
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.hist(T_all[each], bins = temp_bin, range = temp_range)
    axis.plot(KDE_T[each+1][0], KDE_T[each+1][1])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    axis.set_xlim(lower,upper)
    axis.set_xlim(lower,upper)
    axis.set_ylim(0,y_lim_temp)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Temp [K]')
ax7.set_xlabel('Peak Temp [K]')
ax8.set_xlabel('Peak Temp [K]')
ax9.set_xlabel('Peak Temp [K]')
ax10.set_xlabel('Peak Temp [K]')

fig.savefig('{}/Peak Temperature(bin={}).png'.format(dirname,temp_bin))
print('Saved: Peak Temperature(bin={}).png'.format(temp_bin))

#Peak Temperature over each range
fig = plt.figure(figsize=(16, 8)) # peak pressure histograms
ax1=fig.add_subplot(251)
ax2=fig.add_subplot(252)
ax3=fig.add_subplot(253)
ax4=fig.add_subplot(254)
ax5=fig.add_subplot(255)
ax6=fig.add_subplot(256)
ax7=fig.add_subplot(257)
ax8=fig.add_subplot(258)
ax9=fig.add_subplot(259)
ax10=fig.add_subplot(2,5,10)

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
lower,upper = temp_range   
loc = R
for each in range(len(axes)):
    axis = axes[each]
    axis.plot(NKDE_T[each+1][0], NKDE_T[each+1][1])
    axis.set_title('{} km - {} km'.format(loc,loc+R/5))
    axis.set_xlim(lower,upper)
    top = 10/float(temp_range[1]-temp_range[0])
    axis.set_ylim(0,top)
    loc += R/5

ax1.set_ylabel('Frequency of Tracers')
ax6.set_ylabel('Frequency of Tracers')
ax6.set_xlabel('Peak Temp [K]')
ax7.set_xlabel('Peak Temp [K]')
ax8.set_xlabel('Peak Temp [K]')
ax9.set_xlabel('Peak Temp [K]')
ax10.set_xlabel('Peak Temp [K]')

fig.savefig('{}/Normed Temperature(bin={}).png'.format(dirname,temp_bin))
print('Saved: Normed Temperature(bin={}).png\n'.format(temp_bin))


# Changes in the gaussian fit over the 10 panels
locations = []
for each in Loc:
    locations.append(each+R/10)

amps = []
means = []
stdv = []
for i in param:
    amps.append(i[0])
    means.append(i[1])
    stdv.append(i[2])

fig = plt.figure(figsize=(16,8))
ax=fig.add_subplot(111)
ax.scatter(locations,amps)
ax.set_title('Amplitude in local 100km')
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Amplitude')
fig.savefig('{}/10 Amp.png'.format(dirname))

fig = plt.figure(figsize=(16,8))
ax=fig.add_subplot(111)
ax.scatter(locations,means)
ax.set_title('Mean in local 100km')
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Mean [GPa]')
fig.savefig('{}/10 Mean.png'.format(dirname))

fig = plt.figure(figsize=(16,8))
ax=fig.add_subplot(111)
ax.scatter(locations,stdv)
ax.set_title('Standard Deviation in local 100km')
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Standard Deviation [GPa]')
fig.savefig('{}/10 Stdv.png'.format(dirname))

