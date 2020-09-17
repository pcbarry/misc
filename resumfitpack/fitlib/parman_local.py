import sys
from tools.config import load_config, conf
from numpy.random import uniform
import numpy as np

class PARMAN:

    def __init__(self):
        self.get_ordered_free_params()

    def get_ordered_free_params(self):
        self.par=[]
        self.order=[]
        self.pmin=[]
        self.pmax=[]

        if 'check lims' not in conf: conf['check lims']=False

        for k in conf['params']:
            for kk in conf['params'][k]:
                if  conf['params'][k][kk]['fixed']==False:
                    p=conf['params'][k][kk]['value']
                    pmin=conf['params'][k][kk]['min']
                    pmax=conf['params'][k][kk]['max']
                    self.pmin.append(pmin)
                    self.pmax.append(pmax)
                    if p<pmin or p>pmax:
                       if conf['check lims']: raise ValueError('par limits are not consistend with central: %s %s'%(k,kk))

                    self.par.append(p)
                    self.order.append([1,k,kk])

        if 'datasets' in conf:
            for k in conf['datasets']:
                for kk in conf['datasets'][k]['norm']:
                    if  conf['datasets'][k]['norm'][kk]['fixed']==False:
                        p=conf['datasets'][k]['norm'][kk]['value']
                        pmin=conf['datasets'][k]['norm'][kk]['min']
                        pmax=conf['datasets'][k]['norm'][kk]['max']
                        self.pmin.append(pmin)
                        self.pmax.append(pmax)
                        if p<pmin or p>pmax:
                           if conf['check lims']: raise ValueError('par limits are not consistend with central: %s %s'%(k,kk))
                        self.par.append(p)
                        self.order.append([2,k,kk])

        self.pmin=np.array(self.pmin)
        self.pmax=np.array(self.pmax)
        self.par=np.array(self.par)
        self.set_new_params(self.par,initial=True)

    def gen_flat(self,setup=True):
        r=uniform(0,1,len(self.par))
        par=self.pmin + r * (self.pmax-self.pmin)
        if setup: self.set_new_params(par,initial=True)
        return par
        #cnt=0
        #while 1:
        #    cnt+=1
        #    #flag=True
        #    #if  conf['params']['pdf']['g1 N']['value']<0: flag=False
        #    #if flag: break
        #    #if cnt%100==0: print('\n have not found par after %d tries\n'%cnt)
        #return par

    def check_lims(self):
        flag=True
        for k in conf['params']:
            for kk in conf['params'][k]:
                if  conf['params'][k][kk]['fixed']==False:
                    p=conf['params'][k][kk]['value']
                    pmin=conf['params'][k][kk]['min']
                    pmax=conf['params'][k][kk]['max']
                    if  p<pmin or p>pmax:
                        print k,kk, p,pmin,pmax
                        flag=False

        if  'datasets' in conf:
            for k in conf['datasets']:
                for kk in conf['datasets'][k]['norm']:
                    if  conf['datasets'][k]['norm'][kk]['fixed']==False:
                        p=conf['datasets'][k]['norm'][kk]['value']
                        pmin=conf['datasets'][k]['norm'][kk]['min']
                        pmax=conf['datasets'][k]['norm'][kk]['max']
                        if p<pmin or p>pmax:
                          flag=False
                          print k,kk, p,pmin,pmax

        return flag

    def set_new_params(self,parnew,initial=False):
        self.par=parnew
        self.shifts=0
        semaphore={}

        for i in range(len(self.order)):
            ii,k,kk=self.order[i]
            if  ii==1:
                if k not in semaphore: semaphore[k]=0
                if  conf['params'][k][kk]['value']!=parnew[i]:
                    conf['params'][k][kk]['value']=parnew[i]
                    semaphore[k]=1
                    self.shifts+=1
            elif ii==2:
                if kk in conf['datasets'][k]['norm']:
                    if  conf['datasets'][k]['norm'][kk]['value']!=parnew[i]:
                        conf['datasets'][k]['norm'][kk]['value']=parnew[i]
                        self.shifts+=1

        if  initial:
            for k in conf['params']: semaphore[k]=1

        self.propagate_params(semaphore)

    def gen_report(self):
        L=[]
        cnt=0
        for k in conf['params']:
            for kk in sorted(conf['params'][k]):
                if  conf['params'][k][kk]['fixed']==False:
                    cnt+=1
                    if  conf['params'][k][kk]['value']<0:
                        L.append('%d %10s  %10s  %10.5e'%(cnt,k,kk,conf['params'][k][kk]['value']))
                    else:
                        L.append('%d %10s  %10s   %10.5e'%(cnt,k,kk,conf['params'][k][kk]['value']))

        for k in conf['datasets']:
            for kk in conf['datasets'][k]['norm']:
                if  conf['datasets'][k]['norm'][kk]['fixed']==False:
                    cnt+=1
                    L.append('%d %10s %10s %10d  %10.5e'%(cnt,'norm',k,kk,conf['datasets'][k]['norm'][kk]['value']))
        return L

    def propagate_params(self,semaphore):
        flag=False

        #--leading power collinear distributions
        if 'pdf'     in semaphore and semaphore['pdf']==1          : self.set_pdf_params()
        if 'ppdf'    in semaphore and semaphore['ppdf']==1         : self.set_ppdf_params()
        if 'ffpion'  in semaphore and semaphore['ffpion']==1       : self.set_ffpion_params()
        if 'ffkaon'  in semaphore and semaphore['ffkaon']==1       : self.set_ffkaon_params()
        if 'ffhadron'  in semaphore and semaphore['ffhadron']==1   : self.set_ffhadron_params()

        #--power correction terms for idis and pidis
        if 't4F2' in semaphore and semaphore['t4F2']==1: 
            self.set_t4_params('t4F2')
            #--the following lines are needed to fix t4FL/t4F3 = t4F2 as they do not appear in semaphore if they are fixed to another distribution
            if 'same ht' in conf and conf['same ht']==True:
              self.set_t4_params('t4FL')
              self.set_t4_params('t4F3')
        if 't4FL' in semaphore and semaphore['t4FL']==1            : self.set_t4_params('t4FL')
        if 't4F3' in semaphore and semaphore['t4F3']==1            : self.set_t4_params('t4F3')
        if 't4g1' in semaphore and semaphore['t4g1']==1            : self.set_t4_params('t4g1')

        #--offshell correction terms for idis and pidis
        if 'F2off' in semaphore and semaphore['F2off']==1            : self.set_off_params('F2off')
        if 'g1off' in semaphore and semaphore['g1off']==1            : self.set_off_params('g1off')

        #--pion pdfs and spitting functions
        if 'pdf-pion' in semaphore and semaphore['pdf-pion']==1    : self.set_pionpdf_params()
        if 'p->pi,n' in semaphore and semaphore['p->pi,n']==1      : self.set_p_to_pi_n_params()

        #--special pdfs for lattice observables
        if 'pdf(qpdf)' in semaphore and semaphore['pdf(qpdf)']==1  : self.set_qpdf_params()

    def set_params(self,dist,FLAV,PAR):

        #--setup the constraints
        for flav in FLAV:
            for par in PAR:
                if flav+' '+par not in conf['params'][dist]: continue
                if conf['params'][dist][flav+' '+par]['fixed']==True: continue
                if conf['params'][dist][flav+' '+par]['fixed']==False: continue
                reference_flav=conf['params'][dist][flav+' '+par]['fixed']

                if len(reference_flav.split())==2:
                    conf['params'][dist][flav+' '+par]['value']=conf['params'][dist][reference_flav]['value']
                    #print('Fixing %s %s %s to be equal to %s %s'%(dist,flav,par,dist,reference_flav))
                elif len(reference_flav.split())==3:  #allows one to reference from another distribution
                    reference_dist=reference_flav.split()[0]
                    reference_flav=reference_flav.split()[1] + ' ' + reference_flav.split()[2] 
                    conf['params'][dist][flav+' '+par]['value']=conf['params'][reference_dist][reference_flav]['value']
                    #print('Fixing %s %s %s to be equal to %s %s'%(dist,flav,par,reference_dist,reference_flav))
        #--update values at the class
        for flav in FLAV:
            idx=0
            for par in PAR:
                if  flav+' '+par in conf['params'][dist]:
                    conf[dist].params[flav][idx]=conf['params'][dist][flav+' '+par]['value']
                else:
                    conf[dist].params[flav][idx]=0
                idx+=1
        conf[dist].setup()

        #--update values at conf
        for flav in FLAV:
            idx=0
            for par in PAR:
                if  flav+' '+par in conf['params'][dist]:
                    conf['params'][dist][flav+' '+par]['value']= conf[dist].params[flav][idx]
                idx+=1

    def set_pdf_params(self):
        if  conf['pdf parametrization']==0:
            FLAV=['g1','uv1','dv1','ub1','db1','s1','sb1']
            FLAV.extend(['g2','uv2','dv2','ub2','db2','s2','sb2'])
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)
        if  conf['pdf parametrization']==1:
            FLAV=['g1','uv1','dv1','dbmub1','dbpub1','sm1','sp1']
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)
        if  conf['pdf parametrization']==2:
            FLAV=['g1','uv1','dv1','sea1','sea2','db1','ub1','s1','sb1']
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)
        if  conf['pdf parametrization']==3:
            FLAV=['g+1','uv+1','dv+1','sea1+1','sea2+1','db+1','ub+1','s+1','sb+1']
            FLAV.extend(['g-1','uv-1','dv-1','sea1-1','sea2-1','db-1','ub-1','s-1','sb-1'])
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)
            if 'ppdf' in conf:
                conf['ppdf'].params = conf['pdf'].params
                conf['ppdf'].setup()
        if  conf['pdf parametrization']==4:
            FLAV=['g1','uv1','dv1','sea1','sea2','ub1','db1','s1','sb1']
            FLAV.extend(['g2','uv2','dv2','ub2','db2','s2','sb2'])
            FLAV.extend(['g3','uv3','dv3','ub3','db3','s3','sb3'])
            # FLAV.extend(['g2','g3'])
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)
        if  conf['pdf parametrization']==10:
            FLAV=['g1','uv1','dv1','sea1','sea2','ub1','db1','s1','sb1']
            FLAV.extend(['mix'])
            PAR=['N','a','b','c','d']
            self.set_params('pdf',FLAV,PAR)

    def set_ppdf_params(self):
        if  conf['ppdf parametrization']==0:
            FLAV=['g1','up1','dp1','sp1','ub1','db1','sb1']
            PAR=['N','a','b','c','d']
            self.set_params('ppdf',FLAV,PAR)

    def set_ffpion_params(self):
        if  conf['ffpion parametrization']==0:
            FLAV=['g1','u1','d1','s1','c1','b1']
            FLAV.extend(['ub1','db1','sb1','cb1','bb1'])
            PAR=['N','a','b','c','d']
            self.set_params('ffpion',FLAV,PAR)
        if  conf['ffpion parametrization']==1:
            FLAV=['g1','up1','dp1','sp1','cp1','bp1']
            FLAV.extend(['up2','dp2'])
            PAR=['N','a','b','c','d']
            self.set_params('ffpion',FLAV,PAR)

    def set_ffkaon_params(self):
        if  conf['ffkaon parametrization']==0:
            FLAV=['g1','u1','d1','s1','c1','b1']
            FLAV.extend(['ub1','db1','sb1','cb1','bb1'])
            PAR=['N','a','b','c','d']
            self.set_params('ffkaon',FLAV,PAR)
        if  conf['ffkaon parametrization']==1:
            FLAV=['g1','up1','dp1','sp1','cp1','bp1']
            FLAV.extend(['up2','sp2'])
            PAR=['N','a','b','c','d']
            self.set_params('ffkaon',FLAV,PAR)

    def set_ffhadron_params(self):
        if  conf['ffhadron parametrization']=='0h':
            FLAV=['g1','u1','d1','s1','c1','b1']
            FLAV.extend(['ub1','db1','sb1','cb1','bb1'])
            PAR=['N','a','b','c','d']
            self.set_params('ffhadron',FLAV,PAR)
        if  conf['ffhadron parametrization']=='sum2':
            FLAV=['g1','u1','d1','s1','c1','b1']
            FLAV.extend(['ub1','db1','sb1','cb1','bb1'])
            PAR=['N','a','b','c','d']
            self.set_params('ffhadron',FLAV,PAR)

    def set_pionpdf_params(self):
        if conf['pdf-pion parametrization']==0:
            FLAV=['g1','g2','ubv1','ubv2','dv1','dv2','u1','u2','db1','db2','s1','s2','sb1','sb2']
            PAR=['N','a','b','c','d']
            self.set_params('pdf-pion',FLAV,PAR)
        if conf['pdf-pion parametrization']==1:
            FLAV=['g1','g2','g3','ubv1','ubv2','ubv3','dv1','dv2','dv3','u1','u2','u3','db1','db2','db3','s1','s2','s3','sb1','sb2','sb3']
            PAR=['N','a','b','c','d']
            self.set_params('pdf-pion',FLAV,PAR)
        if conf['pdf-pion parametrization']=='sea' or conf['pdf-pion parametrization']=='resum':
            FLAV=['g1','g2','g3','ubv1','ubv2','ubv3','dv1','dv2','dv3','u1','u2','u3','db1','db2','db3','s1','s2','s3','sb1','sb2','sb3']
            PAR=['N','a','b','c','d']
            self.set_params('pdf-pion',FLAV,PAR)

    def set_p_to_pi_n_params(self):
        if 'ln' in conf:
            conf['ln'].L_p2pin=conf['params']['p->pi,n']['lambda']['value']

    def set_t4_params(self,dist):
        if    conf['ht parametrization']==0:                                  FLAV=['proton','neutron']
        elif  conf['ht parametrization']==1 or conf['ht parametrization']==2: FLAV=['proton','deuteron']
        PAR=['N','a','b','c','d']
        self.set_params(dist,FLAV,PAR)

    def set_off_params(self,dist):
        FLAV=['proton','neutron']
        if conf['offshell parametrization']==0 or conf['offshell parametrization']==1: PAR=['N','a','b','c','d']
        if conf['offshell parametrization']==2 or conf['offshell parametrization']==3 or conf['offshell parametrization']==4: PAR=['N','x0','x1']
        self.set_params(dist,FLAV,PAR)

    def check_residual_penalty(self):
        res=np.array([0])
        if 'dy-pion' in conf['datasets']:
            if conf['pdf-pion'].params['u1'][0]<0: res[0]=10000.0
        return res






