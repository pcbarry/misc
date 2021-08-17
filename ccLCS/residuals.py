#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from tools.config       import conf
from tools.residuals    import _RESIDUALS
from obslib.ccLCS.reader import READER
from obslib.ccLCS.theory import CCLCS

class RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='ccLCS'
        self.tabs=conf['ccLCS tabs']
        self.thy=conf['ccLCS']
        self.setup()

    def get_theory(self):
        for k in self.tabs:
            z   = self.tabs[k]['z']
            nu  = self.tabs[k]['nu']
            obs = self.tabs[k]['obs']
            mpi = self.tabs[k]['mpi']
            L   = self.tabs[k]['L']
            a   = self.tabs[k]['a']
            if obs[0] == 'sigma_VA':
                thy =np.array([self.thy.get_ccLCS(nu[i],z[i]**2,mpi[i],L[i],a[i]) for i in range(len(obs))])
                #thy =np.array([self.thy.current_current(nu[i],z[i]**2,mpi[i],L[i],a[i]) for i in range(len(obs))])
            else:
                print('ERR: obs not implemented')
                sys.exit()
            self.tabs[k]['thy']=thy

    def gen_report(self,verb=1,level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point 
        """
        L=[]

        L.append('reaction: ccLCS')
        msg ='%7s'%'idx'
        msg+='%10s'%'tar'
        msg+='%5s'%'npts'
        msg+='%10s'%'chi2'
        msg+='%10s'%'chi2-npts'
        msg+='%10s'%'rchi2'
        msg+='%10s'%'nchi2'
        L.append(msg)
        for idx in self.tabs:
            res =self._get_residuals(idx)
            rres=self._get_rres(idx)
            nres=self._get_nres(idx)

            chi2=np.sum(res**2)
            rchi2=np.sum(rres**2)
            nchi2=nres**2
            npts=len(self.tabs[idx]['value'])
            msg ='%7d'%idx
            msg+='%10s'%'pion'
            msg+='%5d'%npts
            msg+='%10.2f'%chi2
            msg+='%10.2f'%(chi2-npts)
            msg+='%10.2f'%rchi2
            msg+='%10.2f'%nchi2
            L.append(msg)

        """
        if level==1:
            L.append('-'*100)

            msg='obs=%7s,  '
            msg='obs=%7s,  '
            msg='obs=%7s,  '
            msg+='exp=%10.3e,  '
            msg+='thy=%10.3e,  '
            msg+='chi2=%10.3f,  '

            for k in self.tabs:
              if k in conf['datasets']['pITD']['norm']:
                norm=conf['datasets']['pITD']['norm'][k]['value']
              else:
                norm=1
              for i in range(len(self.tabs[k]['value'])):
                  z    = self.tabs[idx]['z'][i]
                  nu   = self.tabs[idx]['nu'][i]
                  a    = self.tabs[idx]['a'][i]
                  mpi  = self.tabs[idx]['mpi'][i]
                  LatL = self.tabs[idx]['L'][i]
                  obs  = self.tabs[idx]['obs'][i]
                  Nf   = self.tabs[idx]['Nf'][i]

              res   = self.tabs[k]['residuals'][i]
              thy   = self.tabs[k]['thy'][i]
              exp   = self.tabs[k]['value'][i]
              #alpha = self.tabs[k]['alpha'][i]
              #rres  = self.tabs[k]['r-residuals'][i]
              #col   = self.tabs[k]['col'][i]
              #shift = self.tabs[k]['shift'][i]
              #reaction = self.tabs[k]['reaction'][i]
              if res<0: chi2=-res**2
              else: chi2=res**2

              L.append(msg%(col,reaction,obs,xFmin,xFmax,pT,exp,alpha,thy,shift,chi2,norm))
        """

        if verb==0:
            return L
        elif verb==1:
            for l in L: print(l)
            return L

if __name__=='__main__':

    from qcdlib import alphaS,aux,eweak,mellin,pdfpion

    conf['Q20']   = 1.0
    conf['order']='NLO'
    conf['alphaSmode']='backward'
    conf['dglap mode']='truncated'
    conf['scheme']='ZMVFS'
    conf['aux']=aux.AUX()
    conf['alphaS']=alphaS.ALPHAS()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)
    conf['pdf-pion']=pdfpion.PDF()
    conf['LCS scale']=3
    
    conf['datasets']={}
    conf['datasets']['ccLCS']={}
    conf['datasets']['ccLCS']['xlsx']={}
    conf['datasets']['ccLCS']['xlsx'][12781]='current-current-pion/expdata/12781.xlsx'
    conf['datasets']['ccLCS']['xlsx'][13581]='current-current-pion/expdata/13581.xlsx'
    conf['datasets']['ccLCS']['norm']={}
    conf['datasets']['ccLCS']['filters']=[]
    conf['ccLCS tabs']=READER().load_data_sets('ccLCS')
    conf['ccLCS']=CCLCS()
    conf['bootstrap']=False

    residuals=RESIDUALS()
    residuals.get_residuals()
    print(residuals.tabs[12781]['thy'])
    print(residuals.tabs[13581]['thy'])



