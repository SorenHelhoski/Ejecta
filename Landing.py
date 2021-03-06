import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin
from scipy.optimize import curve_fit as fit

dirname = 'Landing'
psp.mkdir_p(dirname)

landing_bin = 20 # number of bins of the landing position of ejecta
landing_bin_all = 60 # number of bins of the full landing position of ejecta
pres_bin = 20 # number of bins in the pressure histogram
temp_bin = 30 # number of bins in the temperature histogram

print('Opening data file...')
# open tracer file
file2 = open('Results.dat','r')
x_list = list(eval(file2.readline().replace('\n','')))
y_list = list(eval(file2.readline().replace('\n','')))
w_list = list(eval(file2.readline().replace('\n','')))
v_list = list(eval(file2.readline().replace('\n','')))
b_list = list(eval(file2.readline().replace('\n','')))
p_list = list(eval(file2.readline().replace('\n','')))
T_list = list(eval(file2.readline().replace('\n','')))
m_list = list(eval(file2.readline().replace('\n','')))
t_list = list(eval(file2.readline().replace('\n','')))
r_list = list(eval(file2.readline().replace('\n','')))
used =   list(eval(file2.readline().replace('\n','')))
jumps, start_time, end_time, save_step, R, grid_spacing, g, a, v_0, layers, datafile  = eval(file2.readline().replace('\n',''))
print('DONE\n')

R_list = []
for each in r_list:
    R_list.append(each/R)

landing_bounds = [.9,3.1]
sep_size = (landing_bounds[1] - landing_bounds[0])/landing_bin

guesses = [.14*R**(.74), -3.000] # guesses for the landing loc curve A*x^B
#(must have three sig figs for formatting)

#------------------------------------------------------------
#                          Graphing
#------------------------------------------------------------

# Peak Pressure
fig = plt.figure(figsize=(12, 6)) # peak pressure
ax=fig.add_subplot(111) # ejecta peak pressure
ax.set_xlabel('Landing Location [km]')
ax.set_ylabel('Peak Pressure [GPa]')
ax.set_title('Peak Pressure and Ejecta Landing Location')

ax.scatter(r_list,p_list,s=5,linewidths=0.05)

fig.savefig('{}/Pressure.png'.format(dirname))
print('Saved: Pressure.png')

fig = plt.figure(figsize=(12, 6)) # peak pressure
ax=fig.add_subplot(111) # ejecta peak pressure
ax.set_ylabel('Frequency')
ax.set_xlabel('Pressure [GPa]')
ax.set_yscale('log')
ax.set_title('Pressure Histogram')

ax.hist(p_list, bins = pres_bin)

fig.savefig('{}/Pressure_Hist(bin={}).png'.format(dirname,pres_bin))
print('Saved: Pressure_Hist(bin={}).png'.format(pres_bin))

# Peak Temperature
fig = plt.figure(figsize=(12, 6)) # peak temperature
ax=fig.add_subplot(111) # ejecta peak temperature
ax.set_xlabel('Landing Location [km]')
ax.set_ylabel('Peak Temperature [K]')
ax.set_title('Peak Temperature and Ejecta Landing Location')

ax.scatter(r_list,T_list,s=5,linewidths=0.05)

fig.savefig('{}/Temperature.png'.format(dirname))
print('Saved: Temperature.png')

fig = plt.figure(figsize=(12, 6)) # peak temp
ax=fig.add_subplot(111) # ejecta peak temp
ax.set_xlabel('Frequency')
ax.set_xlabel('Temperature [K]')
ax.set_yscale('log')
ax.set_title('Temperature Histogram')

ax.hist(T_list, bins = temp_bin)

fig.savefig('{}/Temperature_Hist(bin={}).png'.format(dirname,temp_bin))
print('Saved: Temperature_Hist(bin={}).png\n'.format(temp_bin))


#------------------------------------------------------------
#                       Landing Position
#------------------------------------------------------------

print('Binning the Landing Positions every {} km...'.format(int(R*sep_size)))

# order the landing position
R_list.sort()
location = Bin(R_list, weight = w_list, bins = landing_bin, bounds = landing_bounds)

r_bin  = location.get_x()
r_err  = location.get_xerr()
w      = location.get_xerr(factor=2)
vol    = location.get_y()
v_err  = location.get_yerr()

height, h_err = [],[]
width = (landing_bounds[1]-landing_bounds[0])/(landing_bin)
for i in range(len(vol)):
    height.append(vol[i]/(2*np.pi*r_bin[i]*R*width*R))
    h_err.append(v_err[i]/(2*np.pi*r_bin[i]*R*width*R))

location_all = Bin(R_list, weight = w_list, bins = landing_bin_all)

r_bin_all  = location_all.get_x()
r_err_all  = location_all.get_xerr()
w_all      = location_all.get_xerr(factor=2)
vol_all    = location_all.get_y()
v_err_all  = location_all.get_yerr()

print(vol_all)

height_all, h_err_all = [],[]
width = (max(R_list)-min(R_list))/(landing_bin_all)
for i in range(len(vol_all)):
    height_all.append(vol_all[i]/(2*np.pi*r_bin_all[i]*R*width*R))
    h_err_all.append(v_err_all[i]/(2*np.pi*r_bin_all[i]*R*width*R))

cumulative = location_all.get_cmlt()
cmlt_right = location_all.get_cmltr()

def model_func(x,A,B):
    return A * x ** B 

try:
    popt, pcov = fit(model_func, r_bin, height, p0 = guesses, bounds = [[.001,-6],[R/10,-1]], method = 'dogbox')
    A_fit, B_fit = popt 
except:
    A_fit, B_fit = guesses
    print('Could not find appropriate fit')

# fit for Data
height_expected = []
formula_height = []
for each in R_list:
    height_expected.append(model_func(each, A_fit, B_fit))
    formula_height.append(model_func(each, guesses[0], guesses[1]))

print('DONE\n')

# Landing Position
fig = plt.figure(figsize=(12, 6)) # ejecta landing position graph
ax=fig.add_subplot(111) # ejecta landing position
ax.set_xlabel('Normed Landing Location [{}km]'.format(R))
ax.set_ylabel('Ejecta Height [km]')
ax.set_xlim(landing_bounds)
ax.set_yscale('log')
ax.set_title('Ejecta Landing Location')
ax.text((landing_bounds[1]-landing_bounds[0])/2, .7*max(height_expected), 'Curve Fit: {0:0.3}r^({1:0.3f})'.format(A_fit,B_fit))
ax.text((landing_bounds[1]-landing_bounds[0])/2, .5*max(height_expected), 'Theoretical: {0:0.3}r^({1:0.3f})'.format(guesses[0],guesses[1]))

ax.plot(R_list, height_expected)
ax.plot(R_list, formula_height, linestyle = ':')
ax.bar(r_bin, height, width = w, alpha =.1, align = 'center', xerr = r_err, yerr = h_err)

fig.savefig('{}/Landing Location (bin={}).png'.format(dirname,landing_bin))
print('Saved: Landing Location (bin={}).png'.format(landing_bin))

# Landing Position All
fig = plt.figure(figsize=(12, 6)) # ejecta landing position graph
ax=fig.add_subplot(111) # ejecta landing position
ax.set_yscale('log')
#ax.set_xlim(.9, 3.1)
#ax.set_ylim(1, 1.1*max(height_all))
ax.set_xlabel('Normed Landing Location [{}km]'.format(R))
ax.set_ylabel('Ejecta Height [km]')
ax.set_title('Ejecta Landing Location')

ax.bar(r_bin_all, height_all, width = w_all, alpha =.1, align = 'center', xerr = r_err_all, yerr = h_err_all)

fig.savefig('{}/Landing Location (All) (bin={}).png'.format(dirname,landing_bin_all))
print('Saved: Landing Location (All) (bin={}).png'.format(landing_bin_all))

# Landing Position Cumulative
fig = plt.figure(figsize=(12, 6))
ax=fig.add_subplot(111)
ax.set_xlabel('Normed Landing Location [{}km]'.format(R))
ax.set_ylabel('Cumulative Ejecta Fraction')
ax.set_title('Ejecta Landing Location')
ax.plot(r_bin_all, cumulative)
ax.plot([0,max(R_list)], [cumulative[-1],cumulative[-1]], linestyle = ':')
ax.plot([0,max(R_list)], [0,0], linestyle = ':')

fig.savefig('{}/Cumulative Landing Location (bin={}).png'.format(dirname,landing_bin_all))
print('Saved: Cumulative Landing Location (bin={}).png'.format(landing_bin_all))

# Landing Position Cumulative from Right
fig = plt.figure(figsize=(12, 6)) 
ax=fig.add_subplot(111)
ax.set_xlabel('Normed Landing Location [{}km]'.format(R))
ax.set_ylabel('Cumulative Ejecta Fraction')
ax.set_title('Ejecta Landing Location')
ax.plot(r_bin_all, cmlt_right)
ax.plot([0,max(R_list)], [cmlt_right[0],cmlt_right[0]], linestyle = ':')
ax.plot([0,max(R_list)], [0,0], linestyle = ':')

fig.savefig('{}/Cumulative (R) Landing Location (bin={}).png'.format(dirname,landing_bin_all))
print('Saved: Cumulative (R) Landing Location (bin={}).png'.format(landing_bin_all))

