import numpy as np
from math import atan, pi, sqrt, sin, cos, floor
import pySALEPlot as psp

# Changeable Parameters
start_time = 1 # must be larger than 0  
end_time = 600 #
jumps = 1 # 1 over jumps is the fraction of tracers looped over
# Parameters from the simulation (from asteroid.inp)
save_step = 1 # save interval of the simulation
grid_spacing = 1 # km
g = 1.63 *10**-3 # gravity (km/s^2)
a   = 50  # radius of impactor (km)
v_0 = 12  # impact velcity (km/s)
# Rough final size of crater
R = 500 # basin size (km)

if jumps != 1:
      print('1 out of every {} tracers considered'.format(jumps)) 
print('Running ejecta.py from timesteps {}-{}\n'.format(start_time, end_time))
   

#Find the peak pressures of all tracers
print('Extracting Peak Pressure')
peak_file = 'jdata.dat'
model1=psp.opendatfile('../Chicxulub/{}'.format(peak_file))
model1.setScale('km')
step1 = model1.readStep('TrP', model1.nsteps-1)
pres = np.append(np.zeros(0), step1.TrP/10**9)
print('Using timestep {0} of {1} with {2} tracers'.format(model1.nsteps-1, peak_file, len(pres)))
print('DONE\n')

# extract the materials
print('Extracting Materials...')
data_file = 'jdata.dat'
model=psp.opendatfile('../Chicxulub/{}'.format(data_file))
model.setScale('km')
step = model.readStep('TrT', 0)
yy_0 = np.append(np.zeros(0), step.ymark)

#get the materials:
mat = []    # materials of all tracers
for i in range(len(yy_0)): # sort based on initial location
    if yy_0[i] > 0:
        mat.append('impactor')
    elif yy_0[i] > -30:  # bottom of first layer
        mat.append('crust')
    elif yy_0[i] > -500: # bottom of second layer
        mat.append('mantle') 
    else:                # deepest layer
        mat.append('core')
        
print('Using timestep {0} of {1} with {2} tracers'.format(0, data_file, len(mat)))
print('DONE\n')
        
        
print('==============================')
print('Starting Ejecta Search...')
print('==============================')
# Start the ejecta search
model=psp.opendatfile('../Chicxulub/{}'.format(data_file))
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
v_list = [] # velocity (normalized to v_0)
b_list = [] # angles (radians)
p_list = [] # peak pressure (GPa)
T_list = [] # peak temperature (K)
m_list = [] # material list
t_list = [] # launch times (s)
r_list = [] # landing position (km)
used = []   # tracers that have already been used

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
            #if velocity < 15000*(xx_1[each])**-2:
                # added in log regime
                #continue

            used.append(each) # indexes of previously used tracers

            if D_y/save_step < 2.4:
                # rejected tracers of this statement are excluded forever
                # tracers above escape velocity
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

print('==============================')
print('Finalizing Ejecta Search...')
print('==============================')

# order the list in same fashion based on x_list
x_list, y_list, v_list, b_list, p_list, T_list,m_list, t_list, used = zip(*sorted(zip(x_list, y_list, v_list, b_list, p_list, T_list, m_list, t_list, used)))

print('Total Number of Ejecta Tracers: {}'.format(len(used))) 
print('Final Ejecta Launch at Time   : {}'.format(max(t_list)))
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
    file = open('data.txt','x')
except:
    file = open('data.txt','r+')
    file.truncate(0)

prop = [jumps, start_time, end_time, save_step, grid_spacing, g, a, v_0, R]

file.write(str(x_list)+'\n')
file.write(str(y_list)+'\n')
file.write(str(v_list)+'\n')
file.write(str(b_list)+'\n')
file.write(str(p_list)+'\n')
file.write(str(T_list)+'\n')
file.write(str(m_list)+'\n')
file.write(str(t_list)+'\n')
file.write(str(r_list)+'\n')
file.write(str(used)+'\n')
file.write(str(prop)+'\n')

file.close()
print('DONE\n')
      
