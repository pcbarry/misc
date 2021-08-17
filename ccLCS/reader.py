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
        if 'stat_u' in tab: return tab
        c=[x for x in tab if 'stat_' in x]
        alpha=np.zeros(len(tab['value'].values))
        for i in range(len(tab['value'])):
            alpha[i]=tab['stat_%i'%(i+1)][i]**0.5
        tab['alpha_u']=pd.Series(alpha,index=tab.index)
        return tab

    def get_cov(self,tab):
        if 'stat_u' in tab: return tab
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
    conf['datasets']['ccLCS']={}
    conf['datasets']['ccLCS']['xlsx']={}
    conf['datasets']['ccLCS']['xlsx'][22781]='current-current-pion/expdata/22781.xlsx'
    conf['datasets']['ccLCS']['xlsx'][23581]='current-current-pion/expdata/23581.xlsx'
    conf['datasets']['ccLCS']['xlsx'][24131]='current-current-pion/expdata/24131.xlsx'
    conf['datasets']['ccLCS']['xlsx'][24132]='current-current-pion/expdata/24132.xlsx'
    conf['datasets']['ccLCS']['xlsx'][14132]='current-current-pion/expdata/14132.xlsx'
    conf['datasets']['ccLCS']['norm']={}
    conf['datasets']['ccLCS']['filters']=[]
    #conf['datasets']['ccLCS']['filters'].append("z>0.3") 
    #conf['datasets']['ccLCS']['filters'].append("z<0.6") 

    TAB=READER().load_data_sets('ccLCS')

    print(TAB)



