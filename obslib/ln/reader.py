#!/usr/bin/env python
import numpy as np
import pandas as pd
from tools.config import conf
from tools.reader import _READER

class READER(_READER):

  def __init__(self):
    self.aux=conf['aux']

  def get_xpi(self,tab):
    xpi=tab['x']/tab['y']
    tab['xpi']=pd.Series(xpi,index=tab.index)
    return tab

  def modify_table(self,tab):
    tab=self.get_xpi(tab)   
    tab=self.apply_cuts(tab)
    return tab

class READER2(_READER):
  '''
  for db-ub
  '''

  def __init__(self):
    self.aux=conf['aux']

  def modify_table(self,tab):
    tab=self.apply_cuts(tab)
    return tab

if __name__ == "__main__":

  from qcdlib.aux import AUX

  conf['aux']=AUX()
  conf['datasets']={}
  conf['datasets']['ln']={}
  conf['datasets']['ln']['xlsx']={}
  conf['datasets']['ln']['xlsx'][1000]='../../database/ln/expdata/1000.xlsx'
  conf['datasets']['ln']['xlsx'][2000]='../../database/ln/expdata/2000.xlsx'
  conf['datasets']['ln']['filters']=[]
  conf['datasets']['ln']['filters'].append("y<0.3") 

  TAB=READER().load_data_sets('ln')

  conf['datasets']={}
  conf['datasets']['dbub']={}
  conf['datasets']['dbub']['xlsx']={}
  conf['datasets']['dbub']['xlsx'][10001]='../../database/dbub/expdata/10001.xlsx'

  conf['aux']=AUX()
  TAB=READER2().load_data_sets('dbub')

