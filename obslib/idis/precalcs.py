#!/usr/bin/env python
import sys,os
from nuclib import deuterium,helium
from qcdlib import aux, eweak, pdf2, ht0, ht1, alphaS, mellin
from obslib.idis.reader import READER
from obslib.idis.theory import STFUNCS
from tools.config import conf

#--general setup
conf['order'] = 'NLO'
conf['Q20']   = 1.0
conf['dglap mode']='truncated'

conf['aux']   = aux.AUX()
conf['mellin']= mellin.MELLIN(npts=4)
conf['dmellin']= mellin.DMELLIN(nptsN=4,nptsM=4)
conf['alphaS']= alphaS.ALPHAS()
conf['eweak'] = eweak.EWEAK()
conf['pdf']   = pdf2.PDF()
conf['ht']    = False
conf['dsmf']  = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])
conf['hsmf']  = None#helium.HELIUM(root='%s/nuclib/grids/helium'%os.environ['FITPACK'])


tables=[]

def deuteron_datasets():
    tables.append([10011,'idis/expdata/10011.xlsx'])  # deuteron  F2        SLAC    
    tables.append([10017,'idis/expdata/10017.xlsx'])  # deuteron  F2        BCDMS    
    #tables.append([10021,'idis/expdata/10021.xlsx'])  # d/p       F2d/F2p   NMC   
    #tables.append([10002,'idis/expdata/10002.xlsx'])  # deuteron  F2        JLab Hall C (E00-106)
    #tables.append([10033,'idis/expdata/10033.xlsx'])  # n/d       F2        Bonus
    #tables.append([10006,'idis/expdata/10006.xlsx'])  # deuteron  F2        HERMES          
    #tables.append([10050,'idis/expdata/10050.xlsx'])  # plotting data, Q2 = 1.6129          
    #tables.append([10051,'idis/expdata/10051.xlsx'])  # plotting data, Q2 = 10          
    #tables.append([10052,'idis/expdata/10052.xlsx'])  # plotting data, Q2 = 5          

def gen_table(tmc=False,offshell=False,conj=False,ncores=2):
    """
    here we just need to process the deuteron data to add the nuclear effects
    we add npts to the end to remind that is npts mellin points per z interval 
    """

    conf['nuc']=True
    conf['tmc']=tmc
    conf['offshell']=offshell

    npts=4
    conf['mellin']= mellin.MELLIN(npts=npts)
    conf['dmellin']= mellin.DMELLIN(nptsN=npts,nptsM=npts)

    if   tmc==False:   
        conf['path2idistab']  = '%s/grids/grids-idis/notmc%d'%(os.environ['FITPACK'],npts)
    elif tmc=='AOT':
        conf['path2idistab']  = '%s/grids/grids-idis/aot%d'%(os.environ['FITPACK'],npts)
    elif tmc=='GP':
        conf['path2idistab']  = '%s/grids/grids-idis/gp%d'%(os.environ['FITPACK'],npts)

    deuteron_datasets()

    conf['datasets']={}
    conf['datasets']['idis']={}
    conf['datasets']['idis']['xlsx']={}

    for table in tables:
        idx,path=table
        conf['datasets']['idis']['xlsx'][idx]=path

    conf['datasets']['idis']['filters']=[]
    conf['datasets']['idis']['filters'].append("Q2>1.0") 
    conf['datasets']['idis']['filters'].append("W2>4.0")
    conf['stfuncs']=STFUNCS()
    conf['idis tabs']=READER().load_data_sets('idis')
    conf['stfuncs'].gen_tables('idis tabs',tmc,offshell,conj,ncores=ncores)

if __name__=="__main__":
    ncores=20

    #gen_table(tmc=False,offshell=False)
    #gen_table(tmc='GP',offshell=False)
    #gen_table(tmc='AOT',offshell=False)

    #gen_table(tmc=False,offshell=True,conj=False,ncores=ncores)
    #gen_table(tmc='GP',offshell=True,conj=False,ncores=ncores)
    #gen_table(tmc='AOT',offshell=True,conj=False,ncores=ncores)

    #gen_table(tmc=False,offshell=True,conj=True,ncores=ncores)
    #gen_table(tmc='GP',offshell=True,conj=True,ncores=ncores)
    gen_table(tmc='AOT',offshell=True,conj=True,ncores=ncores)




