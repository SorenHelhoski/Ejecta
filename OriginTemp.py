import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin

# timesteps considered (1/skips)
skips = 20

min_temp = 0    # minimum temperature in the colormap [K]
max_temp = 4000 # maximum temperature in the colormap [K]

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
jumps, start_time, end_time, save_step, R, grid_spacing, g, a, v_0, layers, datafile = eval(file2.readline().replace('\n',''))
print('DONE\n')
    
# Make an output directory
dirname='Origin_Temp'
psp.mkdir_p(dirname)

print('Extracting Peak Temperatures...')
# read the last timestep
model1=psp.opendatfile(datafile)
model1.setScale('km')
step1 = model1.readStep('TrT', model1.nsteps-1)
peak = np.append(np.zeros(0), step1.TrT)
print('DONE\n')

# Open the datafile
model=psp.opendatfile(datafile)
model.setScale('km')

# Loop over timesteps
for time in range(0,end_time,skips):
    fig=plt.figure(figsize=(8,4))
    ax=fig.add_subplot(111,aspect='equal')
    
    # Read the step:
    step = model.readStep('TrT', time)
    xx = np.append(np.zeros(0), step.xmark)
    yy = np.append(np.zeros(0), step.ymark)
        
    xx1,yy1 = [],[]
    for i in used:
        xx1.append(xx[i])
        yy1.append(yy[i])
               
    plt.scatter(xx1, yy1, s = .1, c = T_list, alpha = 1, cmap = plt.cm.gist_rainbow, vmin = min_temp, vmax = max_temp)
    plt.scatter(-xx, yy, s = .1, c = peak, alpha = 1, cmap = plt.cm.gist_rainbow, vmin = min_temp, vmax = max_temp)

    plt.colorbar(label='Peak Temperature [K]')

    ax.set_xlabel('[km]')
    ax.set_ylabel('[km]')
    ax.set_xlim(-2*R,2*R)
    ax.set_ylim(-R,R)
       
    ax.set_title('{: 5.2f} s'.format(step.time))
    
    # Save the figure
    fig.savefig('{}/Temp Step-{:05d}.png'.format(dirname,time), dpi = 300)
