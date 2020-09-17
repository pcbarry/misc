#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import pandas as pd

L=open("TASSO").readlines()
L=[l.strip() for l in L]
L=[l for l in L if l!='']
H=L[:8]
table=L[8:]
for l in H: print l
print 
for l in table: print l


table=[[float(_) for _ in l.split()] for l in table]
table=np.transpose(table)

data={}
data['xpmin']=table[0]
data['xpmax']=table[1]
data['value']=table[4]
data['error_u']=table[5]
data['norm_c']=0.045*table[4]
data['hadron']=['hadron' for _ in range(len(table[0]))]
data['col']=['TASSO' for _ in range(len(table[0]))]
data['RS']=[29 for _ in range(len(table[0]))]
data['obs']=["1/sig dsig/dxp" for _ in range(len(table[0]))]
data=pd.DataFrame(data)
data=data[['col','hadron','obs','RS','xpmin','xpmax','value','error_u','norm_c']]
writer = pd.ExcelWriter('TASSO.xlsx')
data.to_excel(writer,'TASSO.xlsx',index=False)
writer.save()


