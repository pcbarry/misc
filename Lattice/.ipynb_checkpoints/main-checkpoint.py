
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

load_config2("analysis-LCS/ccLCS/stat_u_100percent/input.npy")
data=load("/work/JAM/barryp/JAM/analysis-LCS/NLOpy3/msr-inspected/%s"%sys.argv[1]) #--for prior distributions
#--parameters
for k in conf['params']:
    for kk in conf['params'][k]:
        for i in range(len(data['order'][10])):
            if data['order'][10][i][1]=='pdf' and data['order'][10][i][2]==kk:
                conf['params'][k][kk]['value']=data['params'][10][i]
        for i in range(len(data['order'][15])):
            if data['order'][15][i][0]==1:
                if data['order'][15][i][1]==k and data['order'][15][i][2]==kk:
                    conf['params'][k][kk]['value']=data['params'][15][i]
#--normalizations
for k in conf['datasets']:
    for kk in conf['datasets'][k]['norm']:
        for i in range(len(data['order'][15])):
            if data['order'][15][i][0]==2:
                if data['order'][15][i][1]==k and data['order'][15][i][2]==kk:
                    conf['datasets'][k]['norm'][kk]['value']=data['params'][15][i]

ml=MAXLIKE()
ml.run2(hybrid=['pdf','pdf-pion','p->pi,n','dy-pion','ln'])

data={"params":conf["params"],"datasets":conf["datasets"]}
checkdir("analysis-LCS/ccLCS/stat_u_100percent/mcrun")
np.save("analysis-LCS/ccLCS/stat_u_100percent/mcrun/%s"%(sys.argv[1]),data)
