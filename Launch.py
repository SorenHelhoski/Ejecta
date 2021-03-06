import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin
from scipy.optimize import curve_fit as fit

dirname = 'Launch'
psp.mkdir_p(dirname)

# Changeable Parameters
launch_bin = 500 # number of tracers binned together in the launch graphs
mass_bins = 100  # number of bins in the mass density graphs

#(must have three sig figs for formatting)
guesses_velocity = [14000, -2.00] # initial guesses

#densities of materials in kg/m^3
dens0 = 3.25*10**3 # density of impactor
dens1 = 2.5*10**3  # density of crust
dens2 = 3.25*10**3 # density of mantle
dens3 = 0           # density of core

#change to kg/km^3
dens0 = dens0*10**9 
dens1 = dens1*10**9
dens2 = dens2*10**9
dens3 = dens3*10**9           

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

#------------------------------------------------------------
#                     Launch Values
#------------------------------------------------------------

print('Binning every {} tracers for the launch graphs...'.format(launch_bin))
# lists used in the launch graphs
x_bin = [] # average of binned x
x_err = [] # std of each in x_bin
v_bin = [] # average of binned v
v_err = [] # std of each in v_bin
b_bin = [] # average of binned b
b_err = [] # std of each in b_bin

# loops over each bin and gets the average and std error
for each in range(0,len(x_list),launch_bin):
    point = each # start of tracer loop
    # temporary lists of each value pre averaging
    x_ = []
    v_ = []
    b_ = []
    while point < each + launch_bin and point < len(x_list):
        x_.append(x_list[point])
        v_.append(v_list[point])
        b_.append(b_list[point])
        point += 1
    x_bin.append(np.average(x_))
    x_err.append(np.std(x_))
    v_bin.append(np.average(v_))
    v_err.append(np.std(v_))
    b_bin.append(np.average(b_))
    b_err.append(np.std(b_))

print('DONE\n')

#------------------------------------------------------------
#                          Heights
#------------------------------------------------------------

height = []

def max_height(y0,v0,b0,g0):
    return y0+((v0*sin(b0))**2)/(2*g0)

for i in range(len(y_list)):
    height.append(max_height(y_list[i],v_list[i],b_list[i],g))

heights = Bin(height, weight=w_list, bins = 100)

h_x = heights.get_x()
height_spread = heights.get_cmltr()

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.set_title('Launch Height Reached by Tracers')
ax.set_xlabel('Max Parabolic Height')
ax.set_ylabel('Fraction that Reach Height Threshold')
ax.plot([0,max(height)], [0,0], linestyle=':')
ax.plot([0,max(height)], [1,1], linestyle=':')
ax.plot(h_x, height_spread)
fig.savefig('{}/Tracer Height.png'.format(dirname))
print('Saved: Tracer Height.png\n')

#------------------------------------------------------------
#                          Launch Graphs
#------------------------------------------------------------

# Launch Graph Large
fig0 = plt.figure(figsize=(12, 6)) # launch graph
ax0=fig0.add_subplot(111) # scatter velocity
ax0.set_ylabel('Speed [km/s]')
ax0.set_xlabel('Launch Location [km]')
ax0.set_yscale('log')
ax0.set_xscale('log')

from scipy.optimize import curve_fit as fit

def model_func(x,A,B):
    return A* (x) ** B 
try:
    popt, pcov = fit(model_func, x_bin, v_bin, p0 = guesses_velocity, sigma=v_err)
    A_fit, B_fit = popt
except:
    A_fit, B_fit = guesses_velocity
    print('Could not find appropriate fit; Change initial guesses')

ax0.set_title('Launch Speed with curve fit : {0:0.03f}r^{1:0.03f} : {2:0.03f}(r/a)^{1:0.03f}'.format(A_fit,B_fit, A_fit/v_0*a**B_fit))

third_var = t_list # list for the colormap
third_var_name = 'Launch Time'
plt.plot([1,1000],[model_func(1,A_fit,B_fit),model_func(1000,A_fit,B_fit)])
plt.scatter(x_list, v_list,s=1, c = third_var, cmap = plt.cm.gist_rainbow,linewidths=0.01)
plt.colorbar(label = third_var_name)
fig0.savefig('{}/Tracer Launch Speed ({}).png'.format(dirname,third_var_name))

print('Saved: Tracer Launch Speed ({}).png'.format(third_var_name))


# Launch Graphs 2x2
fig = plt.figure(figsize=(12, 6)) # launch graphs
ax1=fig.add_subplot(221) # scatter velocity
ax2=fig.add_subplot(222) # scatter angles
ax3=fig.add_subplot(223) # errorbar velocity
ax4=fig.add_subplot(224) # errorbar angles
ax1.set_title('Launch Speed')
ax2.set_title('Launch Angle')
ax1.set_ylabel('Speed [km/s]')
ax2.set_ylabel('Launch Angle [radians]')
ax3.set_xlabel('Launch Location [km]')
ax3.set_ylabel('Speed [km/s]')
ax4.set_xlabel('Launch Location [km]')
ax4.set_ylabel('Launch Angle [radians]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax4.set_xscale('log')

ax1.scatter(x_list, v_list,s=1,linewidths=0.01)
ax2.scatter(x_list, b_list,s=1,linewidths=0.01)
ax3.errorbar(x_bin, v_bin, xerr = x_err, yerr = v_err)
ax4.errorbar(x_bin, b_bin, xerr = x_err, yerr = b_err)

fig.savefig('{}/Tracer Launch (bin={}).png'.format(dirname, launch_bin))
print('Saved: Tracer Launch (bin={}).png\n'.format(launch_bin))

trajectory = [] #trajector array of all paths
def parabolic_trajectory(x0,y0,v0,b,t0, t_step = 10,t_max = 1050):
    t,y = 0,0
    x_pos = []
    y_pos = []
    vx = v0*np.cos(b)
    vy = v0*np.sin(b)
    while t < t0:
        x_pos.append(x0)
        y_pos.append(0)
        t += t_step
    while y >= 0 and t < t_max:
        x = x0 + vx*(t-t0)
        y = y0-0.5*g*(t-t0)**2+vy*(t-t0)
        x_pos.append(x)
        y_pos.append(y)
        t += t_step
    while t < t_max:
        x_pos.append(x)
        y_pos.append(y)
        t += t_step
    return x_pos,y_pos

print('\nCreating Trajectories...')

for each in range(len(x_list)):
    ptraj = parabolic_trajectory(x_list[each],y_list[each],v_list[each],b_list[each],t_list[each])
    trajectory.append(ptraj)

print('DONE\n')

# Trajectories
fig = plt.figure(figsize=(12, 6))
ax=fig.add_subplot(111)
ax.set_ylabel('Height [km]')
ax.set_xlabel('Distance [km]')
ax.set_aspect('equal')
'''for traj in trajectory:
    ax.plot(traj[0],traj[1], c ='gray', linewidth = .1)
fig.savefig('{}/Tracer Trajectory.png'.format(dirname, launch_bin))'''
print('Saved: Tracer Trajectory.png\n'.format(launch_bin))

#------------------------------------------------------------
#                     Mass Density Dist
#------------------------------------------------------------


x1, x2, x3, x0 = [],[],[],[] # launch pos from each material
w1, w2, w3, w0 = [],[],[],[] # volume of each tracer
for i in range(len(m_list)):
    if m_list[i] == 'impactor':
        x0.append(x_list[i]/R)
        w0.append(w_list[i])
    elif m_list[i] == 'crust':
        x1.append(x_list[i]/R)
        w1.append(w_list[i])
    elif m_list[i] == 'mantle':
        x2.append(x_list[i]/R)
        w2.append(w_list[i])
    else:
        x3.append(x_list[i]/R)
        w3.append(w_list[i])

bins = mass_bins

vol0 = Bin(x0, weight=w0, bins=bins)
vol1 = Bin(x1, weight=w1, bins=bins)
vol2 = Bin(x2, weight=w2, bins=bins)
vol3 = Bin(x3, weight=w3, bins=bins)

freq0 = vol0.get_y(factor = dens0)
freq1 = vol1.get_y(factor = dens1)
freq2 = vol2.get_y(factor = dens2)
freq3 = vol3.get_y(factor = dens3)

xs = vol1.get_x()
mass = []

total_mass = 0
width = (max(xs)-min(xs))/(bins)
for i in range(len(xs)):
    mass_in_bin = freq0[i]+freq1[i]+freq2[i]+freq3[i]
    total_mass += mass_in_bin
    mass.append(10**-6*(mass_in_bin)/(2*np.pi*xs[i]*R*width*R))

xs_fit = []
mass_fit = []
for i in range(len(xs)):
    if xs[i] > 1:
        xs_fit.append(xs[i])
        mass_fit.append(mass[i])

def model_func(x,A,B):
    return A * x ** B 

guesses_mass = [max(mass_fit),-4.30] # initial guesses

try:
    popt, pcov = fit(model_func, xs_fit, mass_fit, p0 = mass_guesses)
    A_fit, B_fit = popt
except:
    A_fit, B_fit = guesses_mass
    print('Could not find appropriate fit; Change initial guesses') 

# fit for Data
curve = []
for i in xs_fit:
    curve.append(model_func(i,A_fit,B_fit))

fig = plt.figure(figsize=(12, 6)) 
ax=fig.add_subplot(111) 
ax.set_xlabel('Normed Distance [{} km]'.format(R))
ax.set_ylabel('Mass per Area [kg m^-2]')
ax.scatter(xs,mass)
ax.scatter(xs_fit,mass_fit)
ax.plot(xs_fit,curve)
ax.text(R,.5*max(mass),'{0:.03}*r^{1:.03} \n\nTotal Mass : {2:.03} kg'.format(A_fit,B_fit,total_mass))
ax.set_yscale('log')
ax.set_title('Ejecta Mass Density and Launch Location')
fig.savefig('{}/Mass Density(bin={}).png'.format(dirname,mass_bins))
print('Saved: Mass Density(bin={}).png\n'.format(mass_bins))
