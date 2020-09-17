#!/usr/bin/env python
import sys,os
import numpy as np
import time
import copy

#--from scipy
from scipy.optimize  import minimize,leastsq
from scipy.optimize  import least_squares

#--from tools
from tools.tools     import checkdir,save,load
import tools.config
from tools.config    import load_config, conf, options
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator

#--from fitlib
from fitlib.resman import RESMAN

class MAXLIKE:
  
    def __init__(self,inputfile,ncores=20):
        self.nworkers=ncores
        self.inputfile=inputfile
    
    def set_counters(self):
        self.chi2tot=1e1000
        self.dchi2=0
        self.t0 = time.time()
        self.cnt=0
  
    def print_status(self,res,rres,nres):

        #--update status parameters  
        shifts=self.parman.shifts
        etime = (time.time()-self.t0)/60
        npts=res.size
        chi2=np.sum(res**2)
        rchi2=np.sum(rres**2)
        nchi2=np.sum(nres**2)
        chi2tot=chi2+rchi2+nchi2
        dchi2=chi2tot-self.chi2tot
        if  shifts>2: 
            if  chi2tot<self.chi2tot:
                self.dchi2=self.chi2tot-chi2tot
                self.chi2tot=chi2tot
  
        #--build header  
        status=[]
        status.append('JAM FITTER')
        status.append('count = %d'%self.cnt)
        status.append('elapsed time(mins)=%f'%etime)
        status.append('shifts  = %d'%shifts)
        status.append('npts    = %d'%npts)
        status.append('chi2    = %f'%chi2)
        status.append('rchi2   = %f'%rchi2)
        status.append('nchi2   = %f'%nchi2)
        status.append('chi2tot = %f'%(chi2tot))
        status.append('dchi2(iter)  = %f'%self.dchi2)
        status.append('dchi2(local) = %f'%dchi2)

        #--special output for pdfs
        if 'pdf'  in conf['params']:
            for _ in conf['pdf'].sr:
                status.append('pdf %s:%f'%(_,conf['pdf'].sr[_]))

            #status.append('proton uvsr = %f'%conf['pdf'].sr['uvsr'])
            #status.append('proton dvsr = %f'%conf['pdf'].sr['dvsr'])
            #status.append('proton msr  = %f'%conf['pdf'].sr['msr'])
            #if 'svsr' in conf['pdf'].sr: status.append('proton svsr = %f'%conf['pdf'].sr['svsr'])

        #--report from resman 
        status.append('')
        status.extend(self.resman.gen_report())

        #--report from parman
        parstatus = self.parman.gen_report()
 
        #--print into screen 
        nstatus=len(status)
        nparstatus=len(parstatus)
        os.system('clear') 
        for i in range(max([nstatus,nparstatus])):
            data=[]
            if i<nstatus: data.append(status[i])
            else: data.append('')
            if i<nparstatus: data.append(parstatus[i])
            else: data.append('')
            print '%-120s  | %s'%tuple(data)
        return status,parstatus
  
    def get_residuals(self,par):
        res,rres,nres=self.resman.get_residuals(par)
        self.cnt+=1
        self.print_status(res,rres,nres)
        if len(rres)!=0: res=np.append(res,rres)
        if len(nres)!=0: res=np.append(res,nres)
        return res

    def checklimits(self):

        for k in conf['params']:
            for kk in conf['params'][k]:
                if conf['params'][k][kk]['fixed']!=False: continue
                p=conf['params'][k][kk]['value']
                pmin=conf['params'][k][kk]['min']
                pmax=conf['params'][k][kk]['max']
                if  p<pmin or p>pmax:
                    print '%s-%s out of limits. '%(k,kk) 
                    sys.exit()

        for k in conf['datasets']:
            for kk in conf['datasets'][k]['norm']:
                p=conf['datasets'][k]['norm'][kk]['value']
                pmin=conf['datasets'][k]['norm'][kk]['min']
                pmax=conf['datasets'][k]['norm'][kk]['max']
                if  p<pmin or p>pmax: 
                    print '%s-%s out of limits. '%(k,kk) 
                    sys.exit()

    def get_bounds(self):
        order=self.parman.order

        bounds_min=[]
        bounds_max=[]
        for entry in order:
            i,k,kk=entry
            if  i==1:
                p=conf['params'][k][kk]['value']
                pmin=conf['params'][k][kk]['min']
                pmax=conf['params'][k][kk]['max']
                if p<pmin or p>pmax: 
                    msg='%s/%s outsize the limits %f %f %f'%(k,kk,p,pmin,pmax)
                    raise ValueError(msg)
                bounds_min.append(conf['params'][k][kk]['min'])
                bounds_max.append(conf['params'][k][kk]['max'])
            elif i==2:
                p=conf['datasets'][k]['norm'][kk]['value']
                pmin=conf['datasets'][k]['norm'][kk]['min']
                pmax=conf['datasets'][k]['norm'][kk]['max']
                if p<pmin or p>pmax: 
                    msg='%s/%s outsize the limits %f %f %f'%(k,kk,p,pmin,pmax)
                    raise ValueError(msg)
                bounds_min.append(conf['datasets'][k]['norm'][kk]['min'])
                bounds_max.append(conf['datasets'][k]['norm'][kk]['max'])

        return (bounds_min, bounds_max)

    def get_guess(self):
        guess=self.parman.par
        return guess

    def gen_summary(self,par):
    
        res,rres,nres    = self.resman.get_residuals(par)
        status,parstatus = self.print_status(res,rres,nres)
        status.extend(parstatus)
        status=[l+'\n' for l in status]
        fname='summary'
        F=open(fname,'w')
        F.writelines(status)
        F.close()

    def gen_output(self):
        """
        modification of the input is done using a dedicated scrip
        at tools/inputmod.py
        """
        inputmod=INPUTMOD(self.inputfile)
  
        for kind in conf['params']: 
            for par in conf['params'][kind]:
                value=conf['params'][kind][par]['value']
                inputmod.mod_par(kind,par,'value',value) 
  
        for reaction in conf['datasets']:
            for idx in conf['datasets'][reaction]['norm']:
                value=conf['datasets'][reaction]['norm'][idx]['value']
                inputmod.mod_norm(reaction,idx,'value',value) 
  
        fname='output.py'
        inputmod.gen_input(fname)

    def run(self):

        global conf
        load_config(self.inputfile)

        #--current value must be within min and max
        self.checklimits()

        #--initialize resman
        self.resman=RESMAN(self.nworkers)
        self.parman=self.resman.parman

        #--setups  
        guess=self.get_guess()
        bounds=self.get_bounds()
        self.set_counters()
        self.parman.set_new_params(guess,initial=True)

        #--run fit         
        fit = least_squares(self.get_residuals, guess,bounds=bounds,method='trf',ftol=conf['ftol'])

        #--final results 
        self.gen_summary(fit.x)
        self.gen_output()

        #--close resman     
        self.resman.shutdown()

if __name__=='__main__':
    
    MAXLIKE('input.py').run()







