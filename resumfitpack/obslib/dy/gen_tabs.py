#!/usr/bin/env python
import os,sys
import numpy as np

import lhapdf

#--from tools
from tools.tools import load,save,checkdir
from tools.config import load_config,conf
#--from qcdlib
from qcdlib import mellin,aux,alphaS,eweak
#-- from dy
import obslib.dy.piontheory
#import obslib.dy.piontheorynlo
#import obslib.dy.piontheoryplus
import obslib.dy.reader

class TUNGSTEN():

    def __init__(self):
        col='EPPS16nlo_CT14nlo_W184'
        Set=0
        self.lha=lhapdf.mkPDF(col,Set)

    def get_pdfs(self,x,Q2):
      p=np.zeros(11)
      p[0]=self.lha.xfxQ2(21,x,Q2)
      p[1]=self.lha.xfxQ2( 2,x,Q2)
      p[2]=self.lha.xfxQ2(-2,x,Q2)
      p[3]=self.lha.xfxQ2( 1,x,Q2)
      p[4]=self.lha.xfxQ2(-1,x,Q2)
      p[5]=self.lha.xfxQ2( 3,x,Q2)
      p[6]=self.lha.xfxQ2(-3,x,Q2)
      p[7]=self.lha.xfxQ2( 4,x,Q2)
      p[8]=self.lha.xfxQ2(-4,x,Q2)
      p[9]=self.lha.xfxQ2( 5,x,Q2)
      p[10]=self.lha.xfxQ2(-5,x,Q2)
      p/=x
      return p

def gen_melltab():
    conf['path2dytab-hybrid']='%s/grids/grids-dypion'%os.environ['FITPACK']

    conf['dy-pion']=obslib.dy.piontheory.DY_PION()
    #conf['dy-pion']=obslib.dy.piontheoryplus.DY_PION_PLUS()
    conf['dy-pion'].mellin=conf['mellin-pion']
    conf['dy-pion'].gen_melltab_hybrid()

if __name__=='__main__':

    conf['datasets']={}
    conf['datasets']['dy-pion']={}
    conf['datasets']['dy-pion']['filters']=[]
    conf['datasets']['dy-pion']['filters'].append("Q2>4.16**2")
    conf['datasets']['dy-pion']['filters'].append("Q2<8.34**2")
    conf['datasets']['dy-pion']['filters'].append("xF>0") 
    conf['datasets']['dy-pion']['filters'].append("xF<0.9") 
    conf['datasets']['dy-pion']['xlsx']={}
    #conf['datasets']['dy-pion']['xlsx'][10001]='dy-pion/expdata/10001.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][10002]='dy-pion/expdata/10002.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][10003]='dy-pion/expdata/10003.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][40001]='dy-pion/expdata/40001.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][40002]='dy-pion/expdata/40002.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][40003]='dy-pion/expdata/40003.xlsx'
    #conf['datasets']['dy-pion']['xlsx'][60000]='dy-pion/expdata/60000.xlsx'
    conf['datasets']['dy-pion']['xlsx'][80000]='dy-pion/expdata/80000.xlsx'
    
    conf['aux']=aux.AUX()
    conf['Q20']=conf['aux'].mc2
    conf['eweak']=eweak.EWEAK()
    conf['order']='NLO'
    conf['alphaS']=alphaS.ALPHAS()
    conf['pdfB']=TUNGSTEN()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)

    conf['dy-pion tabs']=obslib.dy.reader.READER_PIONS().load_data_sets('dy-pion')

    gen_melltab() #--only turn on to change the mellin tabs; already calculated

