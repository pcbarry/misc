#!/usr/bin/env python
import sys,os
import numpy as np
import time

#--from scipy
from scipy.integrate import fixed_quad
from scipy.special import gamma

#--from tools
from tools.bar import BAR
from tools.multiproc import MULTIPROC
from tools.tools import load, save, checkdir
from tools.config import conf
from tools.parallel import PARALLEL

#--from qcdlib
from qcdlib import aux,mellin,alphaS,eweak,pdf0

#--local
from fakepdf import FAKEPDF
from theory import DY
from reader import READER

#====================================
# Not finished yet! Do not use
# Not sure how to get pions to use
# conf['mellin-pion'] while other
# hadrons use conf['mellin'].
#====================================

class DYHYBRID(DY):

  def gen_SIGN(self,N,Q2,S,Y,muF2,part,flav,msg='calculating...'):
    #Nf=self.alphaS.get_Nf(Q2)
    pts=len(N)
    SIGN=np.zeros(pts,dtype=complex)
    #if flav=='c'  and Nf<=3: return SIGN
    #if flav=='cb' and Nf<=3: return SIGN 
    #if flav=='b'  and Nf<=4: return SIGN
    #if flav=='bb' and Nf<=4: return SIGN
    bar=BAR(msg,pts)
    for i in range(pts): 
      self.n=N[i]
      #self.flav=flav
      real=self.get_xsec(Q2,S,Y,muF2,ilum='hybrid-real',part=part)
      imag=self.get_xsec(Q2,S,Y,muF2,ilum='hybrid-imag',part=part)
      SIGN[i]=np.complex(real,imag)
      bar.next()
    bar.finish()
    return SIGN

  def _gen_melltab_hybrid(self,Q2,S,Y,muF2):
    N =self.mellin.N
    data={}

    data['qA,qbB']={}
    for flav in ['ub','db','sb','cb','bb']:
      data['qA,qbB'][flav]=self.gen_SIGN(N,Q2,S,Y,muF2,'qA,qbB',flav,msg='calculating %s %s'%('qA,qbB',flav))

    data['qbA,qB']={}
    for flav in ['u','d','s','c','b']:
      data['qbA,qB'][flav]=self.gen_SIGN(N,Q2,S,Y,muF2,'qbA,qB',flav,msg='calculating %s %s'%('qbA,qB',flav))

    data['qA,gB']={}
    for flav in ['g']:
      data['qA,gB'][flav]=self.gen_SIGN(N,Q2,S,Y,muF2,'qA,gB','g',msg='calculating %s %s'%('qA,gB',flav))

    data['gA,qB']={}
    for flav in ['u','d','s','c','b','ub','db','sb','cb','bb']:
      data['gA,qB'][flav]=self.gen_SIGN(N,Q2,S,Y,muF2,'gA,qB',flav,msg='calculating %s %s'%('gA,qB',flav))

    return data

  def gen_melltab_hybrid(self):
    path2dytab=conf['path2dytab-hybrid']
    for k in conf['DY tabs']:
      path2dytabK='%s/%d'%(path2dytab,k)
      checkdir(path2dytabK)
      npts=len(conf['DY tabs'][k]['values'])
      bar=BAR('hybrid-mell 4 %d'%k,npts)
      for idx in range(conf['DY tabs'][k]['idx']):
        S  = conf['dytab'][k]['S'][idx]
        Y  = conf['dytab'][k]['Y'][idx]
        Q2 = conf['dytab'][k]['Q2'][idx]
        muF2=Q2
        data=self._gen_melltab_hybrid(Q2,S,Y,muF2)
        fname='%s/%d.melltab'%(path2dytabK,idx)
        save(data,fname)
        bar.next()
      bar.finish()

  def load_melltab_hybrid(self):
    self.melltab={}
    path2dytab=conf['path2dytab-hybrid']
    for k in conf['dy-pion tabs']:
      path2dytabK='%s/%d'%(path2dytab,k)
      self.melltab[k]={}
      for idx in conf['dy-pion tabs'][k]['idx']:
        fname='%s/%d.melltab'%(path2dytabK,idx)
        self.melltab[k][idx]=load(fname)

  def get_xsec_mell_hybrid(self,k,i,Q2):
    Nf=conf['alphaS'].get_Nf(Q2)
    idx=conf['dy-pion tabs'][k]['idx'][i]
    data=self.melltab[k][idx]
    conf['pdf-pion'].evolve(Q2)
    PDFA=conf['pdf-pion'].storage[Q2]

    """
    Put in the moments for pdfA (the pion)
    through calling the NPDF class and get_moments(flav,Q2) function
    
    """
    xsec = self.eU2*PDFA['u']  * data['qA,qbB']['ub']\
          +self.eD2*PDFA['d']  * data['qA,qbB']['db']\
          +self.eD2*PDFA['s']  * data['qA,qbB']['sb']\
          +self.eU2*PDFA['ub'] * data['qbA,qB']['u']\
          +self.eD2*PDFA['db'] * data['qbA,qB']['d']\
          +self.eD2*PDFA['sb'] * data['qbA,qB']['s']\
          +self.eU2*PDFA['u']  * data['qA,gB']['g']\
          +self.eD2*PDFA['d']  * data['qA,gB']['g']\
          +self.eD2*PDFA['s']  * data['qA,gB']['g']\
          +self.eU2*PDFA['ub'] * data['qA,gB']['g']\
          +self.eD2*PDFA['db'] * data['qA,gB']['g']\
          +self.eD2*PDFA['sb'] * data['qA,gB']['g']\
          +self.eU2*PDFA['g']  * data['gA,qB']['u']\
          +self.eD2*PDFA['g']  * data['gA,qB']['d']\
          +self.eD2*PDFA['g']  * data['gA,qB']['s']\
          +self.eU2*PDFA['g']  * data['gA,qB']['ub']\
          +self.eD2*PDFA['g']  * data['gA,qB']['db']\
          +self.eD2*PDFA['g']  * data['gA,qB']['sb']

    if Nf>3: 
      xsec+= self.eU2*PDFA['c']  * data['qA,qbB']['cb']\
            +self.eU2*PDFA['cb'] * data['qbA,qB']['c']\
            +self.eU2*PDFA['c']  * data['qA,gB']['g']\
            +self.eU2*PDFA['cb'] * data['qA,gB']['g']\
            +self.eU2*PDFA['g']  * data['gA,qB']['c']\
            +self.eU2*PDFA['g']  * data['gA,qB']['cb']

    if Nf>4: 
      xsec+= self.eD2*PDFA['b']  * data['qA,qbB']['bb']\
            +self.eD2*PDFA['bb'] * data['qbA,qB']['b']\
            +self.eD2*PDFA['b']  * data['qA,gB']['g']\
            +self.eD2*PDFA['bb'] * data['qA,gB']['g']\
            +self.eD2*PDFA['g']  * data['gA,qB']['b']\
            +self.eD2*PDFA['g']  * data['gA,qB']['bb']

    return self.mellin_ext.invert(1,xsec) # x=1 because x**-N is inside data[..][..]

if __name__=='__main__':

  import os
  from qcdlib import aux,mellin,alphaS,eweak
  from fakepdf import FAKEPDF

  conf['Q20']   = 1.0
  conf['alphaSmode']='backward'
  conf['mode']='truncated'
  conf['order']='NLO'
  conf['scheme']='ZMVFS'
  conf['path2dytab']='%s/grids/dytab'%os.environ['FITPACK']
  conf['aux']=aux.AUX()
  conf['mellin']=mellin.MELLIN(npts=4)
  conf['alphaS']=alphaS.ALPHAS()
  conf['eweak']=eweak.EWEAK()
  conf['pdfA']=FAKEPDF()
  conf['pdfB']=FAKEPDF()

  ############################################
  # tests

  #Q2=10.0
  #S=32.0**2
  #Y=0.0
  #muF2=Q2
  #dy=DY()
  #print dy.get_xsec(Q2,S,Y,muF2,ilum='normal',part='full')

  #############################################
  # precalcs

  from qcdlib import aux
  from reader import READER  
  import time


  conf['aux']=aux.AUX()
  conf['datasets']={}
  conf['datasets']['dy']={}
  conf['datasets']['dy']['xlsx']={}
  conf['datasets']['dy']['xlsx'][10001]='dy/expdata/10001.xlsx'
  conf['datasets']['dy']['xlsx'][10002]='dy/expdata/10002.xlsx'
  conf['datasets']['dy']['norm']={}
  conf['datasets']['dy']['filters']=[]
  conf['datasets']['dy']['filters'].append("Q2>1") 


  conf['dy tabs']=READER().load_data_sets('dy')
  dy=DY()
  #dy.gen_melltab()
  dy.load_melltab()

  t1=time.time()
  for k in conf['dy tabs']:
    npts=len(conf['dy tabs'][k]['value'])
    for i in range(npts):
      Q2 = conf['dy tabs'][k]['Q2'][i]
      reaction = conf['dy tabs'][k]['reaction'][i]
      approx = dy.get_mxsec(k,i,Q2,reaction)
      #print i
      S  = conf['dy tabs'][k]['S'][i]
      Y  = conf['dy tabs'][k]['Y'][i]
      muF2=Q2
      exact =  dy.get_xsec(Q2,S,Y,muF2,ilum='normal',part='full')
      rel_err=abs((approx-exact)/exact)*100
      print rel_err


  t2=time.time()
  print t2-t1











