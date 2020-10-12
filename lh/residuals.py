#!/usr/bin/env python
import numpy as np
from tools.config import conf
from tools.residuals import _RESIDUALS
from obslib.lh.reader import READER
from obslib.lh.theory import LH

class RESIDUALS(_RESIDUALS):

    def __init__(self):
        self.reaction='lh'
        self.tabs=conf['lh tabs']
        self.setup()

    def _get_theory(self,entry):
        k,i=entry
        x       =  self.tabs[k]['x'][i]
        xK      =  self.tabs[k]['xK'][i]
        y       =  self.tabs[k]['y'][i]
        Q2      =  self.tabs[k]['Q2'][i]
        kT2max  =  self.tabs[k]['kT2max'][i]
        obs     =  self.tabs[k]['obs'][i]

        if obs=='F2LH (valence)':
            conf['pdf-kaon'].evolve(Q2)
            conf['pdf-kaon'].storage[Q2]['g']=np.zeros(len(conf['pdf-kaon'].storage[Q2]['g']))
            conf['pdf-kaon'].storage[Q2]['qp'][1]=-conf['pdf-kaon'].storage[Q2]['qm'][1]
            conf['pdf-kaon'].storage[Q2]['qp'][2]=conf['pdf-kaon'].storage[Q2]['qm'][2]
            for i in range(3,6):
                conf['pdf-kaon'].storage[Q2]['qp'][i]=np.zeros(len(conf['pdf-kaon'].storage[Q2]['qp'][i]))

            thy=conf['lh'].get_F2LH(xK,y,Q2,kT2max)
            del conf['pdf-kaon'].storage[Q2]

        elif obs=='F2LH (sea)':
            conf['pdf-kaon'].evolve(Q2)
            conf['pdf-kaon'].storage[Q2]['g']=np.zeros(len(conf['pdf-kaon'].storage[Q2]['g']))
            conf['pdf-kaon'].storage[Q2]['qp'][1]=2*conf['pdf-kaon'].storage[Q2]['u']
            conf['pdf-kaon'].storage[Q2]['qp'][2]=2*conf['pdf-kaon'].storage[Q2]['db']
            thy=conf['lh'].get_F2LH(xK,y,Q2,kT2max)
            del conf['pdf-kaon'].storage[Q2]

        elif obs =='F2LH (gluon)':
            conf['pdf-kaon'].evolve(Q2)
            for i in range(1,6):
                conf['pdf-kaon'].storage[Q2]['qp'][i]=np.zeros(len(conf['pdf-kaon'].storage[Q2]['qp'][i]))
            thy=conf['lh'].get_F2LN(xpi,y,Q2,kT2max)
            del conf['pdf-kaon'].storage[Q2]

        elif obs=='F2LH': 
            thy=conf['lh'].get_F2LH(xK,y,Q2,kT2max)
            #print conf['ln'].L_p2pin
        elif obs=='R':  
            deltaxL=self.tabs[k]['deltaxL'][i]
            thy=conf['lh'].get_R(x,xK,deltaxL,y,Q2,kT2max)
        elif obs=='dsig/dxdQ2dxL':
            ye=self.tabs[k]['ye'][i]
            thy=conf['lh'].get_dsigdxdQ2dxL(x,xK,y,Q2,kT2max,ye)
            units=0.389379372e12
            thy*=units
        else:
            print 'ERR: obs not implemented'
            sys.exit()
        return thy

    def gen_report(self,verb=1,level=1):
      """
      verb = 0: Do not print on screen. Only return list of strings
      verv = 1: print on screen the report
      level= 0: only the total chi2s
      level= 1: include point by point 
      """
      L=[]

      L.append('reaction: ln')
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
          col=self.tabs[idx]['col'][0].split()[0]
          tar=self.tabs[idx]['tar'][0]
          npts=res.size
          msg ='%7d'%idx
          msg+='%10s'%tar
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
          msg+='xpi=%10.3e,  '
          msg+='y=%10.3e,  '
          msg+='Q2=%10.3e,  '
          msg+='exp=%10.3e,  ' 
          msg+='alpha=%10.3e,  ' 
          msg+='thy=%10.3e,  ' 
          msg+='shift=%10.3e,  ' 
          msg+='chi2=%10.3f,  '
          msg+='norm=%10.3f,  '
          
          for k in self.tabs:
              if k in conf['datasets']['ln']['norm']:
                  norm=conf['datasets']['ln']['norm'][k]['value']
              else:
                  norm=1
              for i in range(len(self.tabs[k]['value'])):
                  obs   = self.tabs[k]['obs'][i]
                  xpi   = self.tabs[k]['xpi'][i]
                  y     = self.tabs[k]['y'][i]
                  Q2    = self.tabs[k]['Q2'][i]
                  res   = self.tabs[k]['residuals'][i]
                  thy   = self.tabs[k]['thy'][i]
                  exp   = self.tabs[k]['value'][i]
                  alpha = self.tabs[k]['alpha'][i]
                  rres  = self.tabs[k]['r-residuals'][i]
                  col   = self.tabs[k]['col'][i]
                  shift = self.tabs[k]['shift'][i]
                  if res<0: chi2=-res**2
                  else: chi2=res**2
                  L.append(msg%(col,tar,obs,xpi,y,Q2,exp,alpha,thy,shift,chi1,norm))

      if verb==0:
          return L
      elif verb==1:
          for l in L: print l
          return L


if __name__=='__main__':

  import qcdlib
  from qcdlib import aux,eweak,alphaS,mellin,pdfpion0,pdf2
  import obslib
  from obslib import idis
  from obslib.idis import theory
  from reader import READER,READER2
  from theory import STFUNCS,LN

  conf['Q20'] = 1.0
  conf['alphaSmode']='backward'
  conf['dglap mode']='truncated'
  conf['order']='NLO'
  conf['scheme']='ZMVFS'
  conf['ln mode']='model'
  conf['SPLFUNC model']='cov exp'

  conf['aux']=qcdlib.aux.AUX()
  conf['eweak']=qcdlib.eweak.EWEAK()
  conf['alphaS']=qcdlib.alphaS.ALPHAS()
  conf['mellin-pion']=qcdlib.mellin.MELLIN(extended=True,npts=8)
  conf['mellin']=qcdlib.mellin.MELLIN()
  conf['pdf-pion']=qcdlib.pdfpion0.PDF()
  conf['pion-stfuncs']=STFUNCS()
  conf['pdf']=qcdlib.pdf2.PDF()
  conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
  conf['ln']=LN()

  conf['datasets']={}
  conf['datasets']['ln']={}
  conf['datasets']['ln']['xlsx']={}
  conf['datasets']['ln']['xlsx'][1000]='ln/expdata/1000.xlsx'
  conf['datasets']['ln']['xlsx'][2000]='ln/expdata/2000.xlsx'
  conf['datasets']['ln']['filters']=[]
  conf['datasets']['ln']['filters'].append("y<0.3") 
  conf['datasets']['ln']['norm']={}
  conf['datasets']['ln']['norm'][1000]={'value':1,'fixed':True,'min':0,'max':2} 
  conf['datasets']['ln']['norm'][2000]={'value':1,'fixed':True,'min':0,'max':2} 
  conf['ln tabs']=READER().load_data_sets('ln')

  res=RESIDUALS()
  res.get_residuals()
  #res.gen_report(verb=1,level=1)

