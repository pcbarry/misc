#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import pandas as pd

L=open("TPC").readlines()
L=[l.strip() for l in L]
L=[l for l in L if l!='']
H=L[:6]
table=L[6:]
for l in H: print l

table=[[float(_) for _ in l.split()] for l in table]
table=np.transpose(table)

data={}
data['xpmin']=table[0]
data['xpmax']=table[1]
data['value']=table[2]
data['error_u']=table[3]
data['hadron']=['hadron' for _ in range(len(table[0]))]
data['col']=['TPC' for _ in range(len(table[0]))]
data['RS']=[29 for _ in range(len(table[0]))]
data['obs']=["1/sig dsig/dxp" for _ in range(len(table[0]))]
data=pd.DataFrame(data)
data=data[['col','hadron','obs','RS','xpmin','xpmax','value','error_u']]
writer = pd.ExcelWriter('TPC.xlsx')
data.to_excel(writer,'TPC.xlsx',index=False)
writer.save()


