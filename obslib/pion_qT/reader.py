#!/usr/bin/env python
import numpy as np
import pandas as pd
from tools.reader import _READER
from tools.config import conf

class READER(_READER):

    def __init__(self):
        self.aux=conf['aux']

    def get_s(self,tab):
        s=np.ones(len(tab.index))
        s *= 21.8**2
        tab['s']=pd.Series(s,index=tab.index)
        return tab

    def get_xF(self,tab):
        #print(tab.obs.values[0])
        #if tab.obs.values[0] == 'd2sigma/dpTdx':
        if 'd2sigma/dpTdx' in tab.obs.values[0]:
            xF = 0.5*(tab['x_min'].values+tab['x_max'].values)
            tab['xF']=pd.Series(xF,index=tab.index)
            return tab
        else:
            return tab
    
    def get_Q(self,tab):
        #print(tab.obs.values[0])
        #if tab.obs.values[0] == 'd2sigma/dpTdm':
        if 'd2sigma/dpTdm' in tab.obs.values[0]:
            Q = 0.5*(tab['m_min'].values+tab['m_max'].values)
            tab['Q']=pd.Series(Q,index=tab.index)
            #return tab
        #if tab.obs.values[0] == 'd2sigma/dpTdx':
        elif 'd2sigma/dpTdx' in tab.obs.values[0]:
            Q = (4.05+8.55)/2 * (tab['x_max'].values)**0
            tab['Q']=pd.Series(Q,index=tab.index)
        return tab

    def get_Qlims(self,tab):
        #if tab.obs.values[0] == 'd2sigma/dpTdx':
        if 'd2sigma/dpTdx' in tab.obs.values[0]:
            Qmin=np.ones(tab['value'].values.size)*4.05
            Qmax=np.ones(tab['value'].values.size)*8.55
            tab['Qmin']=pd.Series(Qmin,index=tab.index)
            tab['Qmax']=pd.Series(Qmax,index=tab.index)
            return tab
        else:
            return tab

    def get_ymax(self,tab):
        #if tab.obs.values[0] == 'd2sigma/dpTdm':    
        if 'd2sigma/dpTdm' in tab.obs.values[0]:    
            ymax = np.arcsinh(tab['s']**0.5/(2*np.sqrt(tab['Q']**2+tab['pT'].values**2)))
            tab['ymax']=pd.Series(ymax, index=tab.index)
            return tab
        else:
            return tab

    def get_y(self,tab):
        if 'd3sigma/dpTdmdx' in tab.obs.values[0]:
            y = np.arcsinh(tab['xF'].values*tab['s']**0.5/(2*np.sqrt(tab['Q'].values**2+tab['pT'].values**2)))
            tab['y']=pd.Series(y,index=tab.index)
        return tab

    def modify_table(self,tab):
        tab = self.get_s(tab)
        tab = self.get_xF(tab)
        tab = self.get_Q(tab)
        tab = self.get_Qlims(tab)
        tab = self.get_ymax(tab)
        tab = self.get_y(tab)
        tab=self.apply_cuts(tab)
        return tab

if __name__ == "__main__":

    from qcdlib import aux

    conf['aux']=aux.AUX()
    conf['datasets']={}
    conf['datasets']['pion_qT']={}
    conf['datasets']['pion_qT']['xlsx']={}
    conf['datasets']['pion_qT']['xlsx'][1001]='pion_qT/expdata/1001.xlsx'
    conf['datasets']['pion_qT']['xlsx'][1002]='pion_qT/expdata/1002.xlsx'
    conf['datasets']['pion_qT']['norm']={}
    conf['datasets']['pion_qT']['filters']=[]
    conf['datasets']['pion_qT']['filters'].append('pT>2.2')
    TAB=READER().load_data_sets('pion_qT')


    tab=TAB[1001]
    npts=len(tab['value'])
    for i in range(npts):
        msg='pT=%10.2f  m_max=%10.5f m_min=%10.5f'
        print msg%(tab['pT'][i],tab['m_max'][i],tab['m_min'][i])







