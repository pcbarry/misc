
import sys,os
import numpy as np
path = "/work/JAM/barryp/JAM/fitpack2testing/"
sys.path.insert(0, path)
os.environ['FITPACK']=path

from tools.config import conf,load_config2
from tools.tools  import load, save, checkdir
from fitlib.maxlike import MAXLIKE

conf['verbose']   = 100
conf['bootstrap'] = True
conf['ftol']      = 1e-8

load_config2("analysis-LCS/%s/%s/input.npy"%(sys.argv[3],sys.argv[2]))

data=np.load("analysis-LCS/ccLCS/DM/mcrun/%s"%sys.argv[1],allow_pickle=True).item(0)
#--parameters from numpy data
for k in conf['params']:
    if k in data['params']:
        for kk in conf['params'][k]:
             if kk in data['params'][k]:
                conf['params'][k][kk]['value']=data['params'][k][kk]['value']

#--normalizations from numpy data
for k in conf['datasets']:
    if k in data['datasets']:
        for kk in conf['datasets'][k]['norm']:
            if kk in data['datasets'][k]['norm']:
                conf['datasets'][k]['norm'][kk]['value']=data['datasets'][k]['norm'][kk]['value']

ml=MAXLIKE(verbose=1000)
ml.run2()

data={"params":conf["params"],"datasets":conf["datasets"]}
checkdir("analysis-LCS/%s/%s/mcrun"%(sys.argv[3],sys.argv[2]))
np.save("analysis-LCS/%s/%s/mcrun/%s"%(sys.argv[3],sys.argv[2],sys.argv[1]),data)
