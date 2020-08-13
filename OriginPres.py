import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin

# timesteps considered (1/skips)
skips = 20

min_pres = 0   # minimum pressure in the colormap [GPa]
max_pres = 200 # maximum pressure in the colormap [GPa]

print('Opening data file...')
# open tracer file
file2 = open('data.txt','r')
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
jumps, start_time, end_time, save_step, R, grid_spacing, g, a, v_0, layers  = eval(file2.readline().replace('\n',''))
print('DONE\n')
    

# Make an output directory
dirname='Origin_Pressure'
psp.mkdir_p(dirname)

print('Extracting Peak Pressures...')
# read the last timestep
model1=psp.opendatfile('../Chicxulub/jdata.dat')
model1.setScale('km')
step1 = model1.readStep('TrP', model1.nsteps-1)
peak = np.append(np.zeros(0), step1.TrP/10**9)
print('DONE\n')

# Open the datafile
model=psp.opendatfile('../Chicxulub/jdata.dat')
model.setScale('km')

print('Saving Images...')
# Loop over timesteps
for time in range(0,end_time,skips):
    fig=plt.figure(figsize=(8,4))
    ax=fig.add_subplot(111,aspect='equal')
    
    # Read the step:
    step = model.readStep('TrP', time)
    xx = np.append(np.zeros(0), step.xmark)
    yy = np.append(np.zeros(0), step.ymark)

    xx1,yy1 = [],[]
    for i in used:
        xx1.append(xx[i])
        yy1.append(yy[i])
               
    plt.scatter(xx1, yy1, s = .1, c = p_list, alpha = 1, cmap = plt.cm.gist_rainbow, vmin = 0, vmax = max_pres)
    plt.scatter(-xx, yy, s = .1, c = peak, alpha = 1, cmap = plt.cm.gist_rainbow, vmin = 0, vmax = max_pres)

    plt.colorbar(label='Peak Pressure [GPa]')

    ax.set_xlabel('[km]')
    ax.set_ylabel('[km]')    
    ax.set_xlim(-2*R,2*R)
    ax.set_ylim(-R,R)
       
    ax.set_title('{: 5.2f} s'.format(step.time))
    
    # Save the figure
    fig.savefig('{}/Pressure Step-{:05d}.png'.format(dirname,time), dpi = 300)
