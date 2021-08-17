#!/usr/bin/env python
import numpy as np
import pandas as pd
from tools.reader import _READER
from tools.config import conf

class READER(_READER):

    def __init__(self):
        self.aux=conf['aux']

    def get_idx(self,tab):
        tab['idx']=pd.Series(tab.index,index=tab.index)
        return tab

    def get_alpha(self,tab):
        c=[x for x in tab if 'stat_' in x]
        alpha=np.zeros(len(tab['value'].values))
        for i in range(len(tab['value'])):
            alpha[i]=tab['stat_%i'%(i+1)][i]**0.5
        tab['alpha_u']=pd.Series(alpha,index=tab.index)
        return tab

    def get_cov(self,tab):
        c=[x for x in tab if 'stat_' in x]
        for idx in tab['idx']:
            tab['cov_%i'%(idx+1)]=pd.Series(tab['stat_%i'%(idx+1)].values,index=tab.index)
        for i in range(len(c)):
            del tab[c[i]]
        return tab

    def get_pcut(self,tab):
        if 'cut p3small' in conf and conf['cut p3small']:
            if tab['L'].values[0]==24:
                d=tab.query('p<3')
                return d
            else: return tab
        else: return tab

    def modify_table(self,tab):
        tab=self.get_alpha(tab)
        tab=self.apply_cuts(tab)
        tab=self.get_pcut(tab)
        tab=self.get_idx(tab)
        tab=self.get_cov(tab)
        return tab


if __name__ == "__main__":

    from qcdlib import aux

    conf['aux']=aux.AUX()
    conf['datasets']={}
    conf['datasets']['pITD']={}
    conf['datasets']['pITD']['xlsx']={}
    conf['datasets']['pITD']['xlsx'][1001]='pITD-pion/expdata/1001.xlsx'
    conf['datasets']['pITD']['xlsx'][1002]='pITD-pion/expdata/1002.xlsx'
    conf['datasets']['pITD']['norm']={}
    conf['datasets']['pITD']['filters']=[]
    conf['datasets']['pITD']['filters'].append("z>0.3") 
    conf['datasets']['pITD']['filters'].append("z<0.6") 

    TAB=READER().load_data_sets('pITD')

    print(TAB)



