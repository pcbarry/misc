#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd

DATA=[]

DATA.append([1001 ,'pion' ,'TASSO'])
DATA.append([1002 ,'pion' ,'TASSO'])
DATA.append([1003 ,'pion' ,'TASSO'])
DATA.append([1004 ,'pion' ,'TASSO'])
DATA.append([1005 ,'pion' ,'TASSO'])
DATA.append([1006 ,'pion' ,'TASSO'])
DATA.append([1007 ,'pion' ,'TPC'])
DATA.append([1008 ,'pion' ,'TPC'])
DATA.append([1009 ,'pion' ,'TPC(uds)'])
DATA.append([1010 ,'pion' ,'TPC(c)'])
DATA.append([1011 ,'pion' ,'TPC(b)'])
DATA.append([1012 ,'pion' ,'HRS'])
DATA.append([1013 ,'pion' ,'TOPAZ'])
DATA.append([1014 ,'pion' ,'SLD'])
DATA.append([1015 ,'pion' ,'SLD(uds)'])
DATA.append([1016 ,'pion' ,'SLD(c)'])
DATA.append([1017 ,'pion' ,'SLD(b)'])
DATA.append([1018 ,'pion' ,'ALEPH'])
DATA.append([1019 ,'pion' ,'OPAL'])
DATA.append([1020 ,'pion' ,'OPAL(u)'])
DATA.append([1021 ,'pion' ,'OPAL(d)'])
DATA.append([1022 ,'pion' ,'OPAL(s)'])
DATA.append([1023 ,'pion' ,'OPAL(c)'])
DATA.append([1024 ,'pion' ,'OPAL(b)'])
DATA.append([1025 ,'pion' ,'DELPHI'])
DATA.append([1026 ,'pion' ,'DELPHI(uds)'])
DATA.append([1027 ,'pion' ,'DELPHI(b)'])
DATA.append([1028 ,'pion' ,'BABAR'])
DATA.append([1029 ,'pion' ,'BELLE'])
DATA.append([1030 ,'pion' ,'ARGUS'])

DATA.append([2001 ,'kaon' ,'TASSO'])
DATA.append([2002 ,'kaon' ,'TASSO'])
DATA.append([2003 ,'kaon' ,'TASSO'])
DATA.append([2004 ,'kaon' ,'TASSO'])
DATA.append([2005 ,'kaon' ,'TASSO'])
DATA.append([2006 ,'kaon' ,'TASSO'])
DATA.append([2007 ,'kaon' ,'TPC'])
DATA.append([2008 ,'kaon' ,'TPC'])
DATA.append([2012 ,'kaon' ,'HRS'])
DATA.append([2013 ,'kaon' ,'TOPAZ'])
DATA.append([2014 ,'kaon' ,'SLD'])
DATA.append([2015 ,'kaon' ,'SLD(uds)'])
DATA.append([2016 ,'kaon' ,'SLD(c)'])
DATA.append([2017 ,'kaon' ,'SLD(b)'])
DATA.append([2018 ,'kaon' ,'ALEPH'])
DATA.append([2019 ,'kaon' ,'OPAL'])
DATA.append([2020 ,'kaon' ,'OPAL(u)'])
DATA.append([2021 ,'kaon' ,'OPAL(d)'])
DATA.append([2022 ,'kaon' ,'OPAL(s)'])
DATA.append([2023 ,'kaon' ,'OPAL(c)'])
DATA.append([2024 ,'kaon' ,'OPAL(b)'])
DATA.append([2025 ,'kaon' ,'DELPHI'])
DATA.append([2026 ,'kaon' ,'DELPHI(uds)'])
DATA.append([2027 ,'kaon' ,'DELPHI(b)'])
DATA.append([2028 ,'kaon' ,'BABAR'])
DATA.append([2029 ,'kaon' ,'BELLE'])
DATA.append([2030 ,'kaon' ,'ARGUS'])
DATA.append([2031 ,'kaon' ,'DELPHI'])



for data in DATA:
  idx,hadron,col=data
  fname='src2/%d.xlsx'%idx
  tab=pd.read_excel(fname)
  tab['hadron']=pd.Series([hadron for i in range(len(tab.index))],index=tab.index)
  tab['col']=pd.Series([col for i in range(len(tab.index))],index=tab.index)
  writer = pd.ExcelWriter('%d.xlsx'%idx)
  tab.to_excel(writer,'Sheet1',index=False)
  writer.save()




























