#!/usr/bin/env python
import sys,os
import time
import numpy as np
from subprocess import Popen, STDOUT

#--from nuclib
import nuclib.nuclib as nuclib

#--from qcdlib
from qcdlib import pdf0
from qcdlib import ppdf0
from qcdlib import ff0
from qcdlib import aux, eweak,ht, alphaS, mellin

#--from fitlib
from fitlib.parman import PARMAN
from fitlib.resman import RESMAN

#--from tools 
from tools.tools    import checkdir
from tools.config   import conf,load_config


def main():
 
    def timing(): 
        nworkers=20
    
        Sets=[]
        Sets.append('all')
        Sets.append('idis')
        Sets.append('pidis')
        Sets.append('sidis')
        Sets.append('psidis')
        Sets.append('dy')
        Sets.append('sia')
    
        data={}
    
        for keep in Sets:
            #if keep=='all': continue
            load_config('input.py')
            if keep!='all':
              for reaction in Sets[1:]:
                  if reaction==keep: continue
                  del conf['datasets'][reaction]
                  print reaction,conf['datasets'].keys()
            resman=RESMAN(nworkers)
            time=resman.test()
            resman.shutdown()
    
            data[keep]=time
    
    
        print '\nSummary:\n'
        fmt='%10s time=%5.2f  relative=%5.2f'
        for _ in data: print fmt%(_,data[_],data[_]/data['all']*100)
    
    def data_points(): 
        nworkers=20
    
        load_config('input.py')
        resman=RESMAN(nworkers)
    
        data={}    
    
        #--compute residuals
        if 'idis' in conf['datasets']:
            res=resman.idisres.get_residuals(calc=False)[0]
            data['idis']=len(res)
        if 'pidis' in conf['datasets']:
            res=resman.pidisres.get_residuals(calc=False)[0]
            data['pidis']=len(res)
        if 'sidis' in conf['datasets']:
            res=resman.sidisres.get_residuals(calc=False)[0]
            data['sidis']=len(res)
        if 'psidis' in conf['datasets']:
            res=resman.psidisres.get_residuals(calc=False)[0]
            data['psidis']=len(res)
        if 'dy'   in conf['datasets']:
            res=resman.dyres.get_residuals(calc=False)[0]
            data['dy']=len(res)
        if 'sia' in conf['datasets']:
            res=resman.siares.get_residuals(calc=False)[0]
            data['sia']=len(res)
    
    
        total=0
        for _ in data: total+=data[_]
    
        fmt='%10s num points=%8d relative=%5.2f'
        print '\nnumber of data points:\n'
        for _ in data: print fmt%(_,data[_],float(data[_])/total*100)
        print 
        print fmt%('total',total,100.)
        resman.shutdown()
    
    def parameters():
        nworkers=20
    
        load_config('input.py')
        resman=RESMAN(nworkers)
        parman=resman.parman
        par=parman.par
        order=parman.order
    
        data1={}
        data2={}
    
        for _ in order:
    
            Type,Kind,dum=_
            if Type==1 and Kind not in data1: data1[Kind]=0
            if Type==2 and Kind not in data2: data2[Kind]=0
            if Type==1:data1[Kind]+=1
            if Type==2:data2[Kind]+=1
    
    
        print '\nSummary\n'
        print 'total number of parameters: ',len(par) 
        fmt='kind = %10s  num.par= %5d  relative=%5.2f'
        for _ in data1: print fmt%(_,data1[_],data1[_]/float(len(par))*100)
        for _ in data2: print fmt%(_,data2[_],data2[_]/float(len(par))*100)
        resman.shutdown()
    
    def test(): 
        nworkers=32
        load_config('input.py')
        resman=RESMAN(nworkers)
        resman.test(10)
        resman.shutdown()
    
    def cpu_usage():
        DEVNULL = open(os.devnull, 'wb', 0)
        p = Popen(["./metrics.py", "0"], stdout=DEVNULL, stderr=STDOUT)
        ru = os.wait4(p.pid, 0)[2]
        cpu_usage=ru.ru_utime+ru.ru_stime
        print cpu_usage
   
def main(): 
    load_config('input.py')
    nworkers=1
    resman=RESMAN(nworkers)
    resman.test()
    resman.shutdown()

if __name__=="__main__":

    main()



