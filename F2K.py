import sys,os,time
import numpy as np
import copy

import lhapdf

os.environ['LHAPDF_DATA_PATH']='./:'+os.environ['LHAPDF_DATA_PATH']

if __name__=='__main__':

    """
    Initializing the LHAPDF grid and flavors
    """
    lhagrid = lhapdf.mkPDFs('JAM25kaon_nlonll_F2K')
    flavs = [21,-5,-4,-3,-2,-1,1,2,3,4,5]

    """
    Define the replica index (in this example, index = 165), X, and Q2 kinematics
    """
    index = 165
    X = np.linspace(0,1,100)
    Q2 = 10

    """
    Calculate F2K as a sum over all flavors.
    No factors needed as the factors are baked into the grids already.
    """
    F2K = np.array([np.sum([lha[index_check].xfxQ2(f, _,Q2) for f in flavs]) for _ in X])