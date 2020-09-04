import numpy as np
from scipy import integrate
from math import sqrt, pi

# Manual Histogram Class

'''
This class takes a x_list and then bins it into a histogram,
and stores the non-zero bins and the respective frequencies in each

bins is the number of bins

bounds determine the upper and lower limit of counting frequencies
'''

class Bin:
    def __init__(self,x_list, weight=[], bins = 10, bounds = []):

        #set the weights if no value is given
        if weight == []:
            for i in range(len(x_list)):
                weight.append(1)
        
        if x_list != []:
            x_list, weight = zip(*sorted(zip(x_list, weight)))
            if bounds == []:
                bounds = [min(x_list),max(x_list)]
        else:
            bounds = [0,1] # dummy values to keep from crashing when empty

        lower, upper = bounds
        bin_size = (upper - lower) / bins

        # inintialize certain parameters
        bottom = lower # lower bound of bin
        x_bin = [] # average of each
        x_err = [] # std of each bin
        x_ = [] # temporary list of each bin
        y_ = [] # temporary list of weighted frequencies
        freq = [] # frequency of tracers in each bin
        freq_err = []
    
        i = 0
        final = len(x_list)
        for each in range(len(x_list)):
            if x_list[each] < lower:
                i += 1
            if x_list[each] > upper:
                final = each
                break
        # loops over data
        while i < final:
            # append to temporary bin list if within the bounds
            if bottom <= x_list[i] <= bottom + bin_size:
                x_.append(x_list[i])
                y_.append(weight[i])
                i += 1
            else:
                x_bin.append(bottom+bin_size/2)
                freq.append(sum(y_))
                freq_err.append((sum(y_))**(1/2))
                x_err.append(np.std(x_))
                x_ = []
                y_ = []
                bottom += bin_size
        x_bin.append(bottom+bin_size/2)
        freq.append(sum(y_))
        freq_err.append((sum(y_))**(1/2))
        x_err.append(np.std(x_))

        while len(x_bin) <= bins:
            x_bin.append(bottom+bin_size/2)
            freq.append(0)
            freq_err.append(0)
            x_err.append(0)
            bottom += bin_size
        
        c = 0 # total in each bin
        cmlt = [] # cumulative total in each bin
        for each in freq:
            cmlt.append(c)
            c += each

        c = 0 # total in each bin
        cmltr = [] # cumulative total in each bin from the right
        for i in range(len(freq),0,-1):
            cmltr.append(c)
            c += freq[i-1]
        cmltr.reverse()

        self.x_bin    = x_bin
        self.x_err    = x_err
        self.freq     = freq
        self.freq_err = freq_err
        self.cmlt     = cmlt
        self.cmltr    = cmltr
        self.bin_size = bin_size

    def get_cmlt(self): # returns the cumulative percent dist (from left)
        cmlt_ = []
        _max = self.cmlt[-1]
        for each in self.cmlt:
            cmlt_.append(float(each)/_max)
        return cmlt_

    def get_cmltr(self): # returns the cumulative percent dist (from right)
        cmlt_ = []
        _max = self.cmltr[0]
        for each in self.cmltr:
            cmlt_.append(float(each)/_max)
        return cmlt_

    def get_x(self, factor=1): # returns the average value of each bin
        dummy = []
        for each in self.x_bin:
            dummy.append(float(each)*factor)
        return dummy

    def get_xerr(self, factor=1): # returns standard deviation of points in each bin
        dummy = []
        for each in self.x_err:
            dummy.append(float(each)*factor)
        return dummy

    def get_y(self, factor=1): # returns the freq in each bin
        dummy = []
        for each in self.freq:
            dummy.append(float(each)*factor)
        return dummy

    def get_yerr(self, factor=1): # returns approx error of each frequency (Poisson)
        dummy = []
        for each in self.freq_err:
            dummy.append(float(each)*factor)
        return dummy

    def get_norm_y(self, factor=1): # returns the freq in each bin
        dummy = []
        for each in self.freq:
            dummy.append(float(each)*factor)
        dummy0 = []
        for each in dummy:
            dummy0.append(float(each)/sum(dummy))
        return dummy0

#--------------------------------------------------------------------------------

# KDE class

''' 
Input an x_list to be turned into a Kernal Density Estimation.
Each point in the x_list becomes a gaussian centered at that point with amplitude = 1
Then every point is added up
The result is like a continuous histogram
Each gaussian has the same stdv value -> the optimal bandwidth of the x_list

The bounds determine the range over which the KDE is calculated
it will still consider points in x_list outside this range

The resolution determines how many points of the KDE are calculated
This does not alter the KDE, only how many points of it are returned
'''

class Density:
    def __init__(self, x_list, weight=[], bounds = [0,1000], resolution = 100):
        
        #set the weights if no value is given
        if weight == []:
            for i in range(len(x_list)):
                weight.append(1)
        
        lower, upper = bounds
        step_size = (upper-lower)/resolution

        i = lower    
        x_positions = []
        prob_dens = []
        # use Silverman's Rule
        n = float(len(x_list))
        optimal_bandwidth = 1.06*np.std(x_list)*n**(-.2)

        while lower <= i <= upper:
            pdf = 0
            for j in range(len(x_list)):
                pdf += weight[j]*np.exp( (-0.5)*((i-x_list[j])/(optimal_bandwidth))**2 )/(sqrt(2*np.pi))
    
            x_positions.append(i)
            prob_dens.append(pdf)
            i += step_size

        Normal = 0
        for each in prob_dens:
            Normal += step_size * each
        norm_prob_dens = []
        for each in prob_dens:
            norm_prob_dens.append(each/Normal)

        self.xx = x_positions
        self.yy = prob_dens
        self.ny = norm_prob_dens
        
    def get_x(self): # returns the x values of the KDE
        return self.xx
    
    def get_y(self): # returns the y values of the KDE
        return self.yy
    
    def get_norm_y(self, norm = 1): # returns normalized y values of KDE
        return norm*self.ny

        








