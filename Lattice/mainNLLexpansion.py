
import sys,os
import numpy as np
path = "/work/JAM/barryp/JAM/fitpack2/"
sys.path.insert(0, path)
os.environ['FITPACK']=path

from tools.config import conf,load_config2
from tools.tools  import load, save, checkdir
from fitlib.maxlike import MAXLIKE

conf['verbose']   = 100
conf['bootstrap'] = True
conf['ftol']      = 1e-8


load_config2("analysis-LCS/%s/%s/input.npy"%(sys.argv[3],sys.argv[2]))
data=load("/work/JAM/barryp/Notebooks-Jupyterhub/expansioncopy/step18pos/msr-inspected/%s"%sys.argv[1]) #--for prior distributions
#--parameters
for k in conf['params']:
    for kk in conf['params'][k]:
        for i in range(len(data['order'][10])):
            if data['order'][10][i][1]=='pdf' and data['order'][10][i][2]==kk:
                conf['params'][k][kk]['value']=data['params'][10][i]
        for i in range(len(data['order'][18])):
            if data['order'][18][i][0]==1:
                if data['order'][18][i][1]==k and data['order'][18][i][2]==kk:
                    conf['params'][k][kk]['value']=data['params'][18][i]
                    print(data['order'][18][i],data['params'][18][i])
#--normalizations
for k in conf['datasets']:
    for kk in conf['datasets'][k]['norm']:
        for i in range(len(data['order'][18])):
            if data['order'][18][i][0]==2:
                if data['order'][18][i][1]==k and data['order'][18][i][2]==kk:
                    conf['datasets'][k]['norm'][kk]['value']=data['params'][18][i]

ml=MAXLIKE()
ml.run2()

data={"params":conf["params"],"datasets":conf["datasets"]}
checkdir("analysis-LCS/%s/%s/mcrun"%(sys.argv[3],sys.argv[2]))
np.save("analysis-LCS/%s/%s/mcrun/%s"%(sys.argv[3],sys.argv[2],sys.argv[1]),data)
