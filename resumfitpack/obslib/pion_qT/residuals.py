#!/usr/bin/env python
import numpy as np
from tools.config       import conf
from tools.residuals    import _RESIDUALS
from obslib.pion_qT.reader   import READER
from obslib.pion_qT.theory   import PION_QT
import lhapdf
import time

class RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='pion_qT'
        self.tabs=conf['pion_qT tabs']
        self.setup()

    def _get_theory(self,entry):
        idx,i=entry
        pT   = self.tabs[idx]['pT'][i]
        s    = self.tabs[idx]['s'][i]
        mu = pT/2
        mu2=mu**2

        if self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx':

            xF    = self.tabs[idx]['xF'][i]
            Qmin  = self.tabs[idx]['Qmin'][i]
            Qmax  = self.tabs[idx]['Qmax'][i]
            thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx (valence)':

            xF    = self.tabs[idx]['xF'][i]
            Qmin  = self.tabs[idx]['Qmin'][i]
            Qmax  = self.tabs[idx]['Qmax'][i]
            conf['pdf-pion'].evolve(mu2)
            for f in ['g','u','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy1 = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
            del conf['pdf-pion'].storage[mu2]
            conf['pdf-pion'].evolve(mu2)
            for f in ['g','ub','d','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy2 = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
            del conf['pdf-pion'].storage[mu2]

            thy = thy1 - thy2

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx (sea)':

            xF    = self.tabs[idx]['xF'][i]
            Qmin  = self.tabs[idx]['Qmin'][i]
            Qmax  = self.tabs[idx]['Qmax'][i]
            conf['pdf-pion'].evolve(mu2)
            conf['pdf-pion'].storage[mu2]['g']=np.zeros(len(conf['pdf-pion'].storage[mu2]['g']))
            conf['pdf-pion'].storage[mu2]['ub']=conf['pdf-pion'].storage[mu2]['u']
            conf['pdf-pion'].storage[mu2]['d']=conf['pdf-pion'].storage[mu2]['db']
            thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
            del conf['pdf-pion'].storage[mu2]

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx (gluon)':

            xF    = self.tabs[idx]['xF'][i]
            Qmin  = self.tabs[idx]['Qmin'][i]
            Qmax  = self.tabs[idx]['Qmax'][i]
            conf['pdf-pion'].evolve(mu2)
            for f in ['u','ub','d','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
            del conf['pdf-pion'].storage[mu2]

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm':

            Q    = self.tabs[idx]['Q'][i]
            ymax = self.tabs[idx]['ymax'][i]
            thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm (valence)':

            Q    = self.tabs[idx]['Q'][i]
            ymax = self.tabs[idx]['ymax'][i]

            conf['pdf-pion'].evolve(mu2)
            for f in ['g','u','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy1 = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
            del conf['pdf-pion'].storage[mu2]

            conf['pdf-pion'].evolve(mu2)
            for f in ['g','ub','d','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy2 = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
            del conf['pdf-pion'].storage[mu2]
            thy = thy1 - thy2

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm (sea)':

            Q    = self.tabs[idx]['Q'][i]
            ymax = self.tabs[idx]['ymax'][i]

            conf['pdf-pion'].evolve(mu2)
            conf['pdf-pion'].storage[mu2]['g']=np.zeros(len(conf['pdf-pion'].storage[mu2]['g']))
            conf['pdf-pion'].storage[mu2]['ub']=conf['pdf-pion'].storage[mu2]['u']
            conf['pdf-pion'].storage[mu2]['d']=conf['pdf-pion'].storage[mu2]['db']
            thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
            del conf['pdf-pion'].storage[mu2]

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm (gluon)':

            Q    = self.tabs[idx]['Q'][i]
            ymax = self.tabs[idx]['ymax'][i]

            conf['pdf-pion'].evolve(mu2)
            for f in ['u','ub','d','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
            del conf['pdf-pion'].storage[mu2]

        elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm (asy)':

            Q    = self.tabs[idx]['Q'][i]
            ymax = self.tabs[idx]['ymax'][i]

            Qmin = self.tabs[idx]['m_min'][i]
            Qmax = self.tabs[idx]['m_max'][i]

            thy = conf['pion_qT-asy'].get_asy_Q_term(Q,pT,s,ymax)

        elif self.tabs[idx]['obs'][0] == 'd3sigma/dpTdmdx':

            Q   = self.tabs[idx]['Q'][i]
            y   = self.tabs[idx]['y'][i]
            thy = conf['pion_qT_3'].get_xsec(pT,s,Q,y)

        elif self.tabs[idx]['obs'][0] == 'd3sigma/dpTdmdx (valence)':

            Q   = self.tabs[idx]['Q'][i]
            y   = self.tabs[idx]['y'][i]

            conf['pdf-pion'].evolve(mu2)
            for f in ['g','u','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy1 = conf['pion_qT_3'].get_xsec(pT,s,Q,y)
            del conf['pdf-pion'].storage[mu2]
            conf['pdf-pion'].evolve(mu2)
            for f in ['g','ub','d','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy2 = conf['pion_qT_3'].get_xsec(pT,s,Q,y)
            del conf['pdf-pion'].storage[mu2]

            thy = thy1 - thy2

        elif self.tabs[idx]['obs'][0] == 'd3sigma/dpTdmdx (sea)':

            Q   = self.tabs[idx]['Q'][i]
            y   = self.tabs[idx]['y'][i]

            conf['pdf-pion'].evolve(mu2)
            conf['pdf-pion'].storage[mu2]['g']=np.zeros(len(conf['pdf-pion'].storage[mu2]['g']))
            conf['pdf-pion'].storage[mu2]['ub']=conf['pdf-pion'].storage[mu2]['u']
            conf['pdf-pion'].storage[mu2]['d']=conf['pdf-pion'].storage[mu2]['db']
            thy = conf['pion_qT_3'].get_xsec(pT,s,Q,y)
            del conf['pdf-pion'].storage[mu2]

        elif self.tabs[idx]['obs'][0] == 'd3sigma/dpTdmdx (gluon)':

            Q  = self.tabs[idx]['Q'][i]
            y  = self.tabs[idx]['y'][i]

            conf['pdf-pion'].evolve(mu2)
            for f in ['u','ub','d','db','s','sb','c','cb','b','bb']:
                conf['pdf-pion'].storage[mu2][f]=np.zeros(len(conf['pdf-pion'].storage[mu2][f]))
            thy = conf['pion_qT_3'].get_xsec(pT,s,Q,y)
            del conf['pdf-pion'].storage[mu2]

        return thy

    def gen_report(self,verb=1,level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point 
        """
        L=[]

        L.append('reaction: pion_qT')
        msg ='%7s'%'idx'
        msg+='%10s'%'tar'
        msg+='%10s'%'col'
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
            reaction=self.tabs[idx]['reaction'][0]
            col=self.tabs[idx]['col'][0].split()[0]
            npts=res.size
            msg ='%7d'%idx
            msg+='%10s'%reaction
            msg+='%10s'%col
            msg+='%5d'%npts
            msg+='%10.2f'%chi2
            msg+='%10.2f'%(chi2-npts)
            msg+='%10.2f'%rchi2
            msg+='%10.2f'%nchi2
            L.append(msg)

        if level==1:
            L.append('-'*100)  

            msg ='col=%7s,  '
            msg+='tar=%7s,  '
            msg+='obs=%7s,  '
            msg+='xF=%10.3e,  '
            msg+='Y=%10.3e,  '
            msg+='Q2=%10.3e,  '
            msg+='exp=%10.3e,  ' 
            msg+='alpha=%10.3e,  ' 
            msg+='thy=%10.3e,  ' 
            msg+='shift=%10.3e,  ' 
            msg+='chi2=%10.3f,  '
            msg+='norm=%10.3f,  '

            for k in self.tabs:
              if k in conf['datasets']['pion_qT']['norm']:
                norm=conf['datasets']['pion_qT']['norm'][k]['value']
              else:
                norm=1
              for i in range(len(self.tabs[k]['value'])):
                  if self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx':
                      xFmin   = self.tabs[idx]['x_min'][i]
                      xFmax  = self.tabs[idx]['x_max'][i]
                      pT  = self.tabs[idx]['pT'][i]
                      obs = self.tabs[idx]['obs'][i]  
                  elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm':
                      Qmin   = self.tabs[idx]['m_min'][i]
                      Qmax  = self.tabs[idx]['m_max'][i]
                      pT  = self.tabs[idx]['pT'][i]
                      obs = self.tabs[idx]['obs'][i]
                      Q = self.tabs[idx]['Q'][i]

              res   = self.tabs[k]['residuals'][i]
              thy   = self.tabs[k]['thy'][i]
              exp   = self.tabs[k]['value'][i]
              alpha = self.tabs[k]['alpha'][i]
              rres  = self.tabs[k]['r-residuals'][i]
              col   = self.tabs[k]['col'][i]
              shift = self.tabs[k]['shift'][i]
              reaction = self.tabs[k]['reaction'][i]
              if res<0: chi2=-res**2
              else: chi2=res**2
              
              if self.tabs[idx]['obs'][0] == 'd2sigma/dpTdx':
                  L.append(msg%(col,reaction,obs,xFmin,xFmax,pT,exp,alpha,thy,shift,chi2,norm))
              elif self.tabs[idx]['obs'][0] == 'd2sigma/dpTdm':
                  L.append(msg%(col,reaction,obs,Qmin,Qmax,pT,exp,alpha,thy,shift,chi2,norm))

        if verb==0:
            return L
        elif verb==1:
            for l in L: print l
            return L

if __name__=='__main__':

    from qcdlib import alphaS,aux,eweak,mellin,pdfpion0
    from qcdlib import aux,mellin,alphaS,eweak

    conf['Q20']   = 1.0
    conf['order']='NLO'
    conf['alphaSmode']='backward'
    conf['dglap mode']='truncated'
    conf['scheme']='ZMVFS'
    conf['aux']=aux.AUX()
    conf['alphaS']=alphaS.ALPHAS()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)
    conf['pdf-pion']=pdfpion0.PDF()
    
    conf['datasets']={}
    conf['datasets']['pion_qT']={}
    conf['datasets']['pion_qT']['xlsx']={}
    conf['datasets']['pion_qT']['xlsx'][1001]='pion_qT/expdata/1001.xlsx'
    conf['datasets']['pion_qT']['xlsx'][1002]='pion_qT/expdata/1002.xlsx'
    conf['datasets']['pion_qT']['norm']={}
    conf['datasets']['pion_qT']['filters']=[]
    conf['pion_qT tabs']=READER().load_data_sets('pion_qT')
    conf['pion_qT']=PION_QT()

    residuals=RESIDUALS()
    t1=time.time()
    residuals.get_residuals()
    print(conf['pion_qT tabs'][1001]['thy'])
    t2=time.time()
    print t2-t1



