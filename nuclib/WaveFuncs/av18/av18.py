#!/usr/bin/env python
import numpy as np
from scipy.interpolate import griddata

#--u and w are S and D state wave functions in momentum space with units fm^{3/2}
#--normalization integral(dk k^2 (u^2 + w^2)) = 1
#--Ref: Fantoni and Pandharipande, Nucl. Phys. A427 (1984) 473 
#--k in units of 1/fm
def av18(k):
    #--load data for (k,u,w)
    data = np.load('WaveFuncs/av18/av18deut.npy')
    k_data   = np.array(data[0])
    u_data   = np.array(data[1])
    w_data   = np.array(data[2])

    #--interpolate to necessary k value
    u = griddata(k_data,u_data,k,method='cubic',fill_value=0)
    w = griddata(k_data,w_data,k,method='cubic',fill_value=0)

    return u,w



