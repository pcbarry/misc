import numpy as np
import pandas as pd

from .tools import isnumeric,lprint
from .config import conf


class _READER:

  def __init__(self):
    pass

  def apply_cuts(self,tab):
    if 'filters' in conf['datasets'][self.reaction]:
      for f in conf['datasets'][self.reaction]['filters']:
        tab=tab.query(f)
    return tab

  def load_data_sets(self,reaction,verb=True):
    self.reaction=reaction
    if reaction not in conf['datasets']: return None
    XLSX=conf['datasets'][reaction]['xlsx']
    TAB={}
    for k in XLSX: 
      if verb: print 'loading %s data sets %d'%(reaction,k)
      fname=conf['datasets'][reaction]['xlsx'][k]
      tab=pd.read_excel(fname)
      tab=self.modify_table(tab)
      npts=tab.index.size
      if npts==0: continue
      TAB[k]=tab.to_dict(orient='list')
      for kk in TAB[k]: 
        if isnumeric(TAB[k][kk][0]):
          TAB[k][kk]=np.array(TAB[k][kk])
    return TAB

