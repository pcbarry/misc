import sys
import numpy as np
from obslib.idis.reader import READER
from obslib.idis.theory import STFUNCS
from tools.residuals    import _RESIDUALS
from tools.config       import conf

class RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='idis'
        self.tabs=conf['idis tabs']
        self.stfuncs=conf['idis stfuncs']
        self.setup()

        self.tmc = False
        self.hq  = False
        self.ht  = False
        self.Nf  = None
        self.nuc = False
        self.offshell = False
        if 'hq'   in conf: self.hq  = conf['hq']
        if 'tmc'  in conf: self.tmc = conf['tmc']
        if 'ht'   in conf: self.ht  = conf['ht']
        if 'nuc'  in conf: self.nuc = conf['nuc']
        if 'offshell' in conf: self.offshell = conf['offshell']
        if self.hq : self.Nf  = 3 
 
    def _get_sigma_r(self,k,i,F2,FL,F3):
  
        if   '_minus' in self.tabs[k]['lepton beam'][i]: sign = 1
        elif '_plus'  in self.tabs[k]['lepton beam'][i]: sign = -1
        elif 'none'   in self.tabs[k]['lepton beam'][i]: sign = 0
  
        else: 
          print 'ERR in data set %d: lepton id is not identified'%(k)
          sys.exit()
  
        x=self.tabs[k]['X'][i]
        Q2=self.tabs[k]['Q2'][i]
        M2=conf['aux'].M2
  
        if 'Elab' in self.tabs[k]:
            ELab=self.tabs[k]['Elab'][i]
            s=M2+2*ELab*M2**0.5
            y=(Q2/2/x)/((s-M2)/2)
        elif 'Y' in self.tabs[k]:
            y=self.tabs[k]['Y'][i]
        YP=1+(1-y)**2
        YM=1-(1-y)**2
  
        thy=F2-y**2/YP*FL+sign*YM/YP*x*F3
        return thy
  
    def _get_sigma_rcc(self,k,i,W2,WL,W3):
  
        if    '_plus'  in self.tabs[k]['lepton beam'][i]: sign=-1
        elif  '_minus' in self.tabs[k]['lepton beam'][i]: sign=1
        else: 
            print 'ERR in data set %d: lepton id is not identified'%(k)
            sys.exit()
  
        x=self.tabs[k]['X'][i]
        Q2=self.tabs[k]['Q2'][i]
        M2=conf['aux'].M2
  
        if 'Elab' in self.tabs[k]:
            ELab=self.tabs[k]['Elab'][i]
            s=M2+2*ELab*M2**0.5
            y=(Q2/2/x)/((s-M2)/2)
        elif 'Y' in self.tabs[k]:
            y=self.tabs[k]['Y'][i]

        YP=1+(1-y)**2
        YM=1-(1-y)**2
  
        thy=YP/2*W2-y**2/2*WL+sign*YM/2*x*W3
        return thy
  
    def _get_sigma_dxdy(self,k,i,F2,FL,F3):
  
        if    '_plus'  in self.tabs[k]['lepton beam'][i]: sign=-1
        elif  '_minus' in self.tabs[k]['lepton beam'][i]: sign=1
        else: 
          print 'ERR in data set %d: lepton id is not identified'%(k)
          sys.exit()
  
        x=self.tabs[k]['X'][i]
        y=self.tabs[k]['Y'][i]
        Q2=self.tabs[k]['Q2'][i]
        M2=conf['aux'].M2
        YP=1+(1-y)**2
        YM=1-(1-y)**2
  
        thy=2*np.pi*conf['aux'].alfa**2/x/y/Q2
        thy*=(YP+2*x**2*y**2*M2/Q2)*F2-y**2*FL+sign*YM*x*F3
  
        return thy

    def _get_sigma_dxdQ2(self,k,i,F2,FL,F3):
        thy=self._get_sigma_dxdy(k,i,F2,FL,F3)
        x=self.tabs[k]['X'][i]
        rs=self.tabs[k]['RS'][i]
        M2=conf['aux'].M2
        return thy/x/(rs**2-M2)
 
    def get_F2p(self,x,Q2):

        thy = 0

        #--leading twist

        thy+= self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)

        #if  self.hq:
        #    thy+= self.hqstfuncs.get_FX(x,Q2,'F2','c',TMC=False,evolve=False)\
        #        + self.hqstfuncs.get_FX(x,Q2,'F2','b',TMC=False,evolve=False)

        #--higher twist
 
        if  self.ht==True:

            thy+=self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=4)

        return thy

    def get_FLp(self,x,Q2):

        thy = 0

        #--leading twist

        thy+= self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)
  
        #if  self.hq:
        #    thy+= self.hqstfuncs.get_FX(x,Q2,'F2','c',TMC=False,evolve=False)\
        #        + self.hqstfuncs.get_FX(x,Q2,'F2','b',TMC=False,evolve=False)

        #--higher twist
  
        if  self.ht==True: 

            thy+=self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=4)

        return thy

    def get_F3p(self,x,Q2):

        thy = 0

        #--leading twist

        thy+= self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)
  
        #if  self.hq:
        #    thy+= self.hqstfuncs.get_FX(x,Q2,'F2','c',TMC=False,evolve=False)\
        #        + self.hqstfuncs.get_FX(x,Q2,'F2','b',TMC=False,evolve=False)

        #--higher twist
  
        if  self.ht==True: 

            thy+=self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=4)

        return thy

    def get_F2n(self,x,Q2):

        thy = 0

        #--leading twist

        thy+= self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)
  
        #if  self.hq:
        #    thy+= self.hqstfuncs.get_FX(x,Q2,'F2','c',TMC=False,evolve=False)\
        #        + self.hqstfuncs.get_FX(x,Q2,'F2','b',TMC=False,evolve=False)

        #--higher twist
  
        if  self.ht==True: 

            thy+=self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=4)

        return thy

    def get_F2d(self,x,Q2):
        thy=0

        #--leading twist
        if self.nuc==False:

            thy+= self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)\
                + self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)

            #if  self.hq:
            #    thy+= 2.0*self.hqstfuncs.get_FX(x,Q2,'F2','c',tmc=False,evolve=False)\
            #        + 2.0*self.hqstfuncs.get_FX(x,Q2,'F2','b',tmc=False,evolve=False)

        else: 

            thy+=  self.stfuncs.get_F2D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)

            if self.offshell==True:

                thy+= self.stfuncs.get_F2D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,offshell=True,Nf=self.Nf,evolve=False)

            #if  self.hq:
            #    thy+= self.hqstfuncs.get_F2D(x,Q2,'c',tmc=False,evolve=False)\
            #        + self.hqstfuncs.get_F2D(x,Q2,'b',tmc=False,evolve=False)

        #--higher twist

        if self.ht==True:

            if self.nuc==False:

                thy+= self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=4)\
                     +self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=4)
            else:

                thy+=  self.stfuncs.get_F2D(x,Q2,twist=4,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)

        return 0.5*thy

    def get_sig_r_p(self,x,Q2,k,i):

        current=self.tabs[k]['current'][i].strip()
  
        if current=='NC':
  
            F2 = self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)
            FL = self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)
            F3 = self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=2,Nf=self.Nf,tmc=self.tmc,evolve=False)

            if  self.hq: 

                F2 += self.hqstfuncs.get_FX(x,Q2,'F2','c',TMC=False,evolve=False)\
                    + self.hqstfuncs.get_FX(x,Q2,'F2','b',TMC=False,evolve=False)
                FL += self.hqstfuncs.get_FX(x,Q2,'FL','c',TMC=False,evolve=False)\
                    + self.hqstfuncs.get_FX(x,Q2,'FL','b',TMC=False,evolve=False)

            #--ignore power corrections for FL/F3
            if  self.ht:
                F2 += self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=4)
                #FL += self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=4)
                #F3 += self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=4)
  
            thy = self._get_sigma_r(k,i,F2,FL,F3)
  
        if current=='CC':
  
            if   '_plus'  in self.tabs[k]['lepton beam'][i]: sign=+1
            elif '_minus' in self.tabs[k]['lepton beam'][i]: sign=-1
  
            W2 = self.stfuncs.get_WXN(x,Q2,'W2',sign=sign,nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
            WL = self.stfuncs.get_WXN(x,Q2,'WL',sign=sign,nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
            W3 = self.stfuncs.get_WXN(x,Q2,'W3',sign=sign,nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
  
            #--ignore power corrections for WX. Cannot be fit at the moment.
            #if  self.ht: 

                #W2 += self.stfuncs.get_WXN(x,Q2,'W2',sign,nucleon='proton',twist=4)
                #WL += self.stfuncs.get_WXN(x,Q2,'WL',sign,nucleon='proton',twist=4)
                #W3 += self.stfuncs.get_WXN(x,Q2,'W3',sign,nucleon='proton',twist=4)
  
            thy = self._get_sigma_rcc(k,i,W2,WL,W3)

        return thy
    
    def get_sig_r_d(self,x,Q2,k,i):

        current=self.tabs[k]['current'][i].strip()
  
        if  current=='NC':

            if self.nuc==False:
  
                F2 = self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
                FL = self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
                F3 = self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
  
                F2 += self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=2,tmc=self.tmc,evolve=False)
                FL += self.stfuncs.get_FXN(x,Q2,'FL',nucleon='neutron',twist=2,tmc=self.tmc,evolve=False)
                F3 += self.stfuncs.get_FXN(x,Q2,'F3',nucleon='neutron',twist=2,tmc=self.tmc,evolve=False)
  
                if  self.ht: 

                    F2 += self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=4)
                    FL += self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=4)
                    F3 += self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=4)
  
                    F2 += self.stfuncs.get_FXN(x,Q2,'F2',nucleon='neutron',twist=4)
                    FL += self.stfuncs.get_FXN(x,Q2,'FL',nucleon='neutron',twist=4)
                    F3 += self.stfuncs.get_FXN(x,Q2,'F3',nucleon='neutron',twist=4)
  
            if self.nuc==True:

                F2 = self.stfuncs.get_F2D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)
                FL = self.stfuncs.get_FLD(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)
                F3 = self.stfuncs.get_F3D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)

                if  self.ht:

                    F2 += self.stfuncs.get_F2D(x,Q2,twist=4,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)
                    FL += self.stfuncs.get_FLD(x,Q2,twist=4,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)
                    F3 += self.stfuncs.get_F3D(x,Q2,twist=4,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,evolve=False)

                if self.offshell:

                    #--no offshell corrections for FL or F3 at the moment
                    F2 += self.stfuncs.get_F2D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,offshell=True,evolve=False)
                    #FL += self.stfuncs.get_FLD(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,offshell=True,evolve=False)
                    #F3 += self.stfuncs.get_F3D(x,Q2,twist=2,tmc=self.tmc,nuc=self.nuc,Nf=self.Nf,offshell=True,evolve=False)

            thy = 0.5*self.get_sigma_r(k,i,F2,FL,F3)

        else:

            raise ValueError('%s not supported'%current)

        return thy

    def get_F2d_over_F2p(self,x,Q2):

        F2p=self.get_F2p(x,Q2)
        F2d=self.get_F2d(x,Q2)

        return F2d/F2p

    def get_F2n_over_F2d(self,x,Q2):

        F2n=self.get_F2n(x,Q2)
        F2d=2*self.get_F2d(x,Q2)
        return F2n/F2d
 
    def get_sigma_dxdy(self,x,Q2,k,i):
        F2 = self.stfuncs.get_FXN(x,Q2,'F2',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
        FL = self.stfuncs.get_FXN(x,Q2,'FL',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
        F3 = self.stfuncs.get_FXN(x,Q2,'F3',nucleon='proton',twist=2,tmc=self.tmc,evolve=False)
        return self._get_sigma_dxdy(k,i,F2,FL,F3)

    def _get_theory(self,entry):

        idx,i=entry
        x=self.tabs[idx]['X'][i]
        Q2=self.tabs[idx]['Q2'][i]
        target=self.tabs[idx]['target'][i]
        obs=self.tabs[idx]['obs'][i].strip()

        conf['pdf'].evolve(Q2)
  
        if   obs=='F2'    and target=='p':  thy=self.get_F2p(x,Q2)
        elif obs=='F2'    and target=='d':  thy=self.get_F2d(x,Q2)
        elif obs=='FL'    and target=='p':  thy=self.get_FLp(x,Q2)
        elif obs=='F3'    and target=='p':  thy=self.get_F3p(x,Q2)
        elif obs=='sig_r' and target=='p':  thy=self.get_sig_r_p(x,Q2,idx,i)
        elif obs=='sig_r' and target=='d':  thy=self.get_sig_r_d(x,Q2,idx,i)
        elif obs=='F2d/F2p'              :  thy=self.get_F2d_over_F2p(x,Q2)
        elif obs=='F2n/F2d'              :  thy=self.get_F2n_over_F2d(x,Q2)
        elif obs=='dsig/dxdy'            :  thy=get_sigma_dxdy(x,Q2,idx,i)
        elif obs=='dsig/dxdQ2'           :  thy=get_sigma_dxdQ2(x,Q2,idx,i)
        else:
            msg='exp=%d obs=%s and target=%s not implemented'
            raise ValueError(msg%(k,obs,target))  

        return thy
  
    def gen_report(self,verb=1,level=1):
      """
      verb = 0: Do not print on screen. Only return list of strings
      verv = 1: print on screen the report
      level= 0: only the total chi2s
      level= 1: include point by point 
      """

      L=[]

      if len(self.tabs.keys())!=0:
        L.append('reaction: unpol DIS')
        for f in conf['datasets']['idis']['filters']:
          L.append('filters: %s'%f)

        L.append('%7s %3s %10s %5s %10s %10s %10s %10s'%('idx','tar','col','npts','chi2','chi2-npts','rchi2','nchi2'))
        for k in self.tabs:
          if len(self.tabs[k])==0: continue 
          res=self.tabs[k]['residuals']

          rres=[]
          for c in conf['rparams']['idis'][k]:
            rres.append(conf['rparams']['idis'][k][c]['value'])
          rres=np.array(rres)

          if k in conf['datasets']['idis']['norm']:
            norm=conf['datasets']['idis']['norm'][k]
            nres=(norm['value']-1)/norm['dN']
          else:
            nres=0

          
          chi2=np.sum(res**2)
          rchi2=np.sum(rres**2)
          nchi2=nres**2
          tar=self.tabs[k]['target'][0]
          col=self.tabs[k]['col'][0].split()[0]
          npts=res.size
          L.append('%7d %3s %10s %5d %10.2f %10.2f %10.2f %10.2f'%(k,tar,col,npts,chi2,chi2-npts,rchi2,nchi2))

        if level==1:
          L.append('-'*100)  
          for k in self.tabs:
            if len(self.tabs[k]['value'])==0: continue 
            if k in conf['datasets']['idis']['norm']:
              norm=conf['datasets']['idis']['norm'][k]
              nres=(norm['value']-1)/norm['dN']
              norm=norm['value']
            else:
              norm=1.0
              nres=0
            for i in range(len(self.tabs[k]['value'])):
              x     = self.tabs[k]['X'][i]
              Q2    = self.tabs[k]['Q2'][i]
              res   = self.tabs[k]['residuals'][i]
              thy   = self.tabs[k]['thy'][i]
              exp   = self.tabs[k]['value'][i]
              alpha = self.tabs[k]['alpha'][i]
              rres  = self.tabs[k]['r-residuals'][i]
              col   = self.tabs[k]['col'][i]
              shift = self.tabs[k]['shift'][i]
              tar   = self.tabs[k]['target'][i]
              msg='col=%7s, tar=%5s, x=%10.3e, Q2=%10.3e, exp=%10.3e, alpha=%10.3e, thy=%10.3e, shift=%10.3e, chi2=%10.3e, res=%10.3e, norm=%10.3e, '
              L.append(msg%(col,tar,x,Q2,exp,alpha,thy,shift,res**2,res,norm))

      if verb==0:
        return L
      elif verb==1:
        for l in L: print l
        return L


