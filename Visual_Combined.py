import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from math import atan, pi, sqrt, sin
from aux_library import Bin

psp.mkdir_p('Combined/Landing')
psp.mkdir_p('Combined/Material')
psp.mkdir_p('Combined/Pressure')
psp.mkdir_p('Combined/Temperature')

pres_lim = [0,200]
temp_lim = [0,4000]

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

tmax = int(end_time*save_step)
    
for each in range(4): 
    try:
        answer = raw_input('Material, Pressure, Temperature, Landing ? ')
        answer = answer.lower()
        answer = answer[0]
        break
    except:
        pass

if answer == 'm':
    h_map = 'Material'
elif answer == 'p':
    h_map = 'Pressure'
elif answer == 't':
    h_map = 'Temperature'
elif answer == 'l':
    h_map = 'Landing'
else:  
    h_map = 'Material'
    print('Defaulting to Material')

print('Extracting Materials...')
# Open the datafile
model=psp.opendatfile(datafile)
model.setScale('km')
step = model.readStep('TrP', 0)
yy = np.append(np.zeros(0), step.ymark)
col = []
for i in range(len(yy)): # find materials based on starting height
    if yy[i] > layers[0]:
        co = 'black'
    elif yy[i] > layers[1]:
        co = 'darkgrey'
    elif yy[i] > layers[2]:
        co = 'y'
    else:
        co = 'r'
    col.append(co)
print('DONE\n')

col1 = []
for i in range(len(m_list)): # change materials to their colors
    if m_list[i] == 'impactor':
        color = 'black'
    elif m_list[i] == 'crust':
        color = 'darkgrey'
    elif m_list[i] == 'mantle':
        color = 'y'
    else:
        color = 'r'
    col1.append(color)

print('Extracting Peak Pressures...')
# read the last timestep
model1=psp.opendatfile(datafile)
model1.setScale('km')
step1 = model1.readStep('TrP', model1.nsteps-1)
peak = np.append(np.zeros(0), step1.TrP/10**9)
print('DONE\n')

# Open the datafile
model=psp.opendatfile(datafile)
model.setScale('km')

#------------------------------------------------------------
#                     Trajectories
#------------------------------------------------------------

trajectory = [] #trajector array of all paths
def parabolic_trajectory(x0,y0,v0,b,t0, t_step = 1,t_max = 1000):
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
    ptraj = parabolic_trajectory(x_list[each],y_list[each],v_list[each],b_list[each],t_list[each],t_max = tmax)
    trajectory.append(ptraj)

print('DONE\n')

#------------------------------------------------------------
#                          Graphing
#------------------------------------------------------------

print('Saving Images...')

t_s = 10
for t in range(0,tmax-t_s,t_s):

    fig=plt.figure(figsize=(8,4))
    ax=fig.add_subplot(111,aspect='equal')
 
    # Read the step:
    if h_map == 'Temperature':
        step = model.readStep('TrT', int(t/save_step))
        value = np.append(np.zeros(0), step.TrT)
    else:
        step = model.readStep('TrP', int(t/save_step))
        value = np.append(np.zeros(0), step.TrP)

    xx = np.append(np.zeros(0), step.xmark)
    yy = np.append(np.zeros(0), step.ymark)

    if h_map == 'Pressure':
        ax.scatter(-xx, yy, c = value, cmap = plt.cm.gist_rainbow, s = .1, vmin = pres_lim[0], vmax = pres_lim[1])
    elif h_map == 'Temperature':
        ax.scatter(-xx, yy, c = value, cmap = plt.cm.gist_rainbow, s = .1, vmin = temp_lim[0], vmax = temp_lim[1])
    else:
        ax.scatter(-xx, yy, c = col, s = .1)

    ax.set_xlabel('[km]')
    ax.set_ylabel('[km]')    
    ax.set_xlim(-2*R,2*R)
    ax.set_ylim(-R,R)
       
    ax.set_title('{: 5.2f} s'.format(step.time))

    Xs, Ys = [],[]
    for traj in trajectory:
        Xs.append(traj[0][t])
        Ys.append(traj[1][t])

    if h_map == 'Material':
        ax.scatter(Xs, Ys, c = col1, s = .1)
    elif h_map == 'Pressure':
        ax.scatter(Xs, Ys, c = p_list, cmap = plt.cm.gist_rainbow, s = .1, vmin = pres_lim[0], vmax = pres_lim[1])
    elif h_map == 'Temperature':
        ax.scatter(Xs, Ys, c = T_list, cmap = plt.cm.gist_rainbow, s = .1, vmin = temp_lim[0], vmax = temp_lim[1])
    elif h_map == 'Landing':
        ax.scatter(Xs, Ys, c = r_list, cmap = plt.cm.gist_rainbow, s = .1)
    else:
        print('Nothing is being plotted')
 
    fig.savefig('Combined/{0}/{1:04d}.png'.format(h_map,t))


