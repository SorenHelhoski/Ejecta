import numpy as np
from math import atan, pi, sqrt, sin, cos, floor
import pySALEPlot as psp

#============================================================
#                      Custom Parameters 
#============================================================

input_file = open('./Ejecta.inp','r')
param = {'STS':None,'FTS':None,'FCS':None,'ESC':None,'CRT':None,'SIM':None}
for line in input_file:
    row = line.split(' ')
    if row[0] in param.keys():
        param[row[0]] = eval(row[-1])

datafile = param['SIM']

start_time = param['STS'] # Starting time STEP (must be larger than 0)
end_time = param['FTS'] # final time STEP
R = param['FCS'] # Rough final size of crater (km) appeal to plots of the impact

# exclusion parameters (not stored anywhere but here)
esc_vel = param['FCS'] # escape velocity (km/s)
crit_angle = param['CRT'] # minimum angle to be considered uplift (rad)

jumps = 1 # number of skips between searches. (use 1, unless troubleshooting)
# terminal message
if jumps != 1:
      print('1 out of every {} tracers considered'.format(jumps)) 
print('Running ejecta_search.py from timesteps {}-{}\n'.format(start_time, end_time))

peak_file = datafile # data file that is read to get peak pressures (at single timestep)
data_file = datafile # data file used to extract all other information

#============================================================
#                Parameters from asteroid.inp)
#============================================================

print('Opening asteroid.inp...')
# open input file
ast_input = open('../asteroid.inp','r')
keywords = ['GRIDSPC', 'GRAV_V', 'OBJRESH', 'OBJVEL', 'DTSAVE']
ast_dict = {}
for line in ast_input:
    word  = line[0:16].replace(' ','')
    value = (line[54:-1].replace(' ','').replace(':',',')).replace('D','*10**')
    if word == 'LAYPOS':
        layers0 = eval('['+value+']')
    if word == 'TR_SPCH' or word == "TR_SPCV":
        tr_spacing = eval('['+value+']')
        for each in tr_spacing:
            if each != -1:
                print('Tracer Spacing must be set to -1')
                raise
    if word in keywords:
        ast_dict[word] = eval(value)

save_step = ast_dict['DTSAVE']             # save interval of sim (s)
grid_spacing = ast_dict['GRIDSPC'] *.001   # (km)
g   = ast_dict['GRAV_V'] * -0.001          # gravity (km/s^2)
a   = ast_dict['OBJRESH'] * grid_spacing   # radius of impactor (km)
v_0 = ast_dict['OBJVEL'] * -0.001          # impact velcity (km/s)

layers = []
base = layers0[-1]*grid_spacing
for i in layers0:
    layers.append(i*grid_spacing - base)
layers.reverse()
layers.append(-1*10**16)
layers.append(-1*10**16)

print('DONE\n')   

# make sure Results.dat exists
open('Results.dat','a').close()
filetest = open('Results.dat','r')    
filetest.close() 
   
#============================================================
#                 Peak Values Search Starts Here
#============================================================


#Find the peak pressures of all tracers
print('Extracting Peak Pressure')
model1=psp.opendatfile(peak_file)
model1.setScale('km')
step1 = model1.readStep('TrP', model1.nsteps-1)
pres = np.append(np.zeros(0), step1.TrP/10**9)
print('Using timestep {0} of {1} with {2} tracers'.format(model1.nsteps-1, peak_file, len(pres)))
print('DONE\n')

# extract the materials and weights
print('Extracting Materials...')
model=psp.opendatfile(data_file)
model.setScale('km')
step = model.readStep('TrT', 0)
xx_0 = np.append(np.zeros(0), step.xmark)
yy_0 = np.append(np.zeros(0), step.ymark)

#get the materials:
mat = []    # materials of all tracers
weight = [] # weights of all the tracers
for i in range(len(yy_0)): # sort based on initial location
    if yy_0[i] > layers[0]:
        mat.append('impactor')
    elif yy_0[i] > layers[1]:  # bottom of first layer
        mat.append('crust')
    elif yy_0[i] > layers[2]: # bottom of second layer
        mat.append('mantle') 
    else:                # deepest layer
        mat.append('core')
    weight.append(2*np.pi*xx_0[i]*grid_spacing**2)
        
print('Using timestep {0} of {1} with {2} tracers'.format(0, data_file, len(mat)))
print('DONE\n')
        
        
print('==============================')
print('Starting Ejecta Search...')
print('==============================')
# Start the ejecta search
model=psp.opendatfile(data_file)
model.setScale('km')
print('Using {} with {} tracers in each timestep'.format(data_file, len(mat)))

estimate_time = 0.00000012 * len(mat)/jumps * (end_time - start_time) 
# constant times tracers times timesteps
hours = int(floor(estimate_time))
minutes = int(60*(estimate_time % 1))

print('Estimated Time : {} hours, {} minutes'.format(hours, minutes))
step = model.readStep('TrT', start_time-1)
xx_0 = np.append(np.zeros(0), step.xmark)
yy_0 = np.append(np.zeros(0), step.ymark)

# inintialize lists
x_list = [] # x position (km)
y_list = [] # y position (km)
w_list = [] # volume of material for that tracer [km3]
v_list = [] # velocity (normalized to v_0)
b_list = [] # angles (radians)
p_list = [] # peak pressure (GPa)
T_list = [] # peak temperature (K)
m_list = [] # material list
t_list = [] # launch times (s)
r_list = [] # landing position (km)
used = []   # tracers that have already been used

excluded = [0,0,0] # number of rejected tracers
                   # both, above escape, above angle

# Loop over Tracers for Considered Timesteps

for time in range(start_time , end_time+1):
    # Read the step
    try:
        step = model.readStep('TrT', time)
    except:
        print('Timestep {} OUT OF RANGE'.format(time))
        end_time = time
        break
    temp = np.append(np.zeros(0), step.TrT)
    xx_1 = np.append(np.zeros(0), step.xmark)
    yy_1 = np.append(np.zeros(0), step.ymark)
    if time == 0:
        print(len(xx_1))
    for each in range(0,len(xx_1),jumps):
        # Compare with previous step
        if (each not in used) and yy_1[each]> 0:
            D_x = xx_1[each] - xx_0[each] # Change in x
            D_y = yy_1[each] - yy_0[each] # Change in y
            angle = atan(D_y / D_x) # launch angle
            distance = sqrt((D_x)**2 + (D_y)**2) # change in distance
            velocity = distance/save_step # launch velocity
            if D_x <= 0 or D_y <= 0:
                # pass over tracers that are traveling down or in
                continue

            used.append(each) # indexes of previously used tracers

            if D_y/save_step > esc_vel and angle > crit_angle:
                excluded[0] += 1
            elif D_y/save_step > esc_vel:
                excluded[1] += 1
            elif angle > crit_angle:
                excluded[2] += 1
            else:
                # rejected tracers of this statement are excluded forever
                # tracers above escape velocity and above critical angle
                x_list.append(xx_1[each]) # x position normalized
                y_list.append(yy_1[each]) # y position
                v_list.append(velocity) # velocity normalized
                b_list.append(angle) # angle in radians 
                t_list.append(time*save_step) # ejection time
                T_list.append(temp[each]) # peak temperature          
    xx_0, yy_0 = xx_1, yy_1 # makes the current step the previous step

for each in used:
    p_list.append(pres[each])
    m_list.append(mat[each])
    w_list.append(weight[each])

print('==============================')
print('Finalizing Ejecta Search...')
print('==============================')

# order the list in same fashion based on x_list
x_list, y_list, w_list, v_list, b_list, p_list, T_list,m_list, t_list, used = zip(*sorted(zip(x_list, y_list, w_list, v_list, b_list, p_list, T_list, m_list, t_list, used)))

print('Total Number of Ejecta Tracers: {}'.format(len(used))) 
print('Final Ejecta Launch at Time   : {}'.format(max(t_list)))
print('                      Timestep: {}'.format(int(max(t_list)/save_step)))
print('Above escape velocity     : {}'.format(excluded[0]+excluded[1])) 
print('Above critical angle      : {}'.format(excluded[0]+excluded[2])) 
print('overlap of velocity+angle : {}'.format(excluded[0])) 
print('DONE\n')
      
print('Calculating Landing Positions...')
for i in range(len(x_list)):
    vx = v_list[i]*np.cos(b_list[i])
    vy = v_list[i]*np.sin(b_list[i])
    time = vy/g + sqrt(vy**2+2*g*y_list[i])/g
    pos = x_list[i] + vx*time
    r_list.append(pos)
print('DONE\n')     
      
print('Writing Output File...')
try:
    file0 = open('Results.dat','x')
except:
    file0 = open('Results.dat','r+')
    file0.truncate(0)

prop = [jumps, start_time, end_time, save_step, R, grid_spacing, g, a, v_0, layers,datafile]

file0.write(str(x_list)+'\n')
file0.write(str(y_list)+'\n')
file0.write(str(w_list)+'\n')
file0.write(str(v_list)+'\n')
file0.write(str(b_list)+'\n')
file0.write(str(p_list)+'\n')
file0.write(str(T_list)+'\n')
file0.write(str(m_list)+'\n')
file0.write(str(t_list)+'\n')
file0.write(str(r_list)+'\n')
file0.write(str(used)+'\n')
file0.write(str(prop)+'\n')

file0.close()
print('DONE\n')
      
