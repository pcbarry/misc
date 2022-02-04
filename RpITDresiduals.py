#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from tools.config       import conf
from tools.residuals    import _RESIDUALS
from obslib.pITD.reader import READER
from obslib.pITD.theory import PITD

class RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='pITD'
        self.tabs=conf['pITD tabs']
        self.thy=conf['pITD']
        self.setup_pITD()

    def get_theory(self):
        for k in self.tabs:
            z   = self.tabs[k]['z']
            nu  = self.tabs[k]['nu']
            obs = self.tabs[k]['obs']
            mpi = self.tabs[k]['mpi']
            L   = self.tabs[k]['L']
            a   = self.tabs[k]['a']
            if obs[0] == 'Re(rpITD)':
                thy =np.array([self.thy.get_ReITD_values(nu[i],z[i]**2,mpi[i],L[i],a[i]) for i in range(len(obs))])
            elif obs[0] == 'Im(rpITD)':
                thy =np.array([self.thy.ImITD(nu[i],z[i]**2) for i in range(len(obs))]) #--needs work still!!
            else:
                print('ERR: obs not implemented')
                sys.exit()
            self.tabs[k]['thy']=thy

    def resample_pITD(self):
        self.tabs=copy.deepcopy(self.original)
        for k in self.tabs:
            npts=len(self.tabs[k]['value'])
            c=[x for x in self.tabs[k] if 'cov' in x]
            cov=np.zeros((len(c),len(c)))
            for i in range(len(c)):
                cov[i]=self.tabs[k][c[i]]
            w,v=np.linalg.eig(cov)
            v=v.T
            pseudo=np.einsum('ij,i->j',v,np.sqrt(w)*np.random.randn(len(w)))
            self.tabs[k]['value']+=pseudo
        
    def _get_pITD_residuals(self,k):

        npts=len(self.tabs[k]['value'])
        exp=self.tabs[k]['value']

        thy=self.tabs[k]['thy'] #--assuming no norm
        self.tabs[k]['prediction']=thy
        N=np.ones(exp.size)

        c=[x for x in self.tabs[k] if 'cov' in x]
        cov=np.zeros((len(c),len(c)))
        for i in range(len(c)):
            cov[i]=self.tabs[k][c[i]]

        covinv=np.linalg.inv(cov)
        vect=exp-thy
        chi2=np.einsum('i,ij,j',vect,covinv,vect)
        res=chi2**0.5

        self.tabs[k]['residuals']=np.zeros(len(self.tabs[k]['value']))
        self.tabs[k]['residuals'][0]=res
        return res

    def get_pITD_residuals(self):

        res,rres,nres=[],[],[] #--currently, no need for rres or nres
        self.get_theory()

        for k in self.tabs:
            res=np.append(res,self._get_pITD_residuals(k))

        return res,rres,nres

    def setup_pITD(self):
        self.percent_to_absolute()
        self.add_columns()
        self.get_alpha()
        self.retrieve_norm_uncertainty()
        self.setup_rparams()
        self.original=copy.deepcopy(self.tabs)
        if 'bootstrap' in conf and conf['bootstrap']: self.resample_pITD()

    def gen_report(self,verb=1,level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point 
        """
        L=[]

        L.append('reaction: pITD')
        msg ='%7s'%'idx'
        msg+='%10s'%'tar'
        msg+='%5s'%'npts'
        msg+='%10s'%'chi2'
        msg+='%10s'%'chi2-npts'
        msg+='%10s'%'rchi2'
        msg+='%10s'%'nchi2'
        L.append(msg)
        for idx in self.tabs:
            res =self._get_pITD_residuals(idx)
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
    conf['imell']=mellin.IMELLIN()
    conf['pdf-pion-int']=pdfpion.PDF(mellin=conf['imell'])
    conf['pITD scale']=3
    
    conf['datasets']={}
    conf['datasets']['pITD']={}
    conf['datasets']['pITD']['xlsx']={}
    conf['datasets']['pITD']['xlsx'][1001]='pITD-pion/expdata/1001.xlsx'
    conf['datasets']['pITD']['xlsx'][1002]='pITD-pion/expdata/1002.xlsx'
    conf['datasets']['pITD']['norm']={}
    conf['datasets']['pITD']['filters']=[]
    conf['datasets']['pITD']['filters'].append('z>0.3')
    conf['datasets']['pITD']['filters'].append('z<0.6')
    conf['pITD tabs']=READER().load_data_sets('pITD')
    conf['pITD']=PITD()
    conf['bootstrap']=True

    residuals=RESIDUALS()
    residuals.get_pITD_residuals()
    print(conf['pITD tabs'][1001]['thy'])
    print(conf['pITD tabs'][1002]['thy'])



