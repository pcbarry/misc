#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import pandas as pd

L=open("ALEPH").readlines()
L=[l.strip() for l in L]
L=[l for l in L if l!='']
H=L[:8]
table=L[8:]
for l in H: print l

table=[[float(_) for _ in l.split()] for l in table]
table=np.transpose(table)

data={}
data['xpmin']=table[0]
data['xpmax']=table[1]
data['value']=table[2]
data['stat_u']=table[3]
data['syst_u']=table[4]
data['hadron']=['hadron' for _ in range(len(table[0]))]
data['col']=['ALEPH' for _ in range(len(table[0]))]
data['RS']=[91.2 for _ in range(len(table[0]))]
data['obs']=["1/sig dsig/dxp" for _ in range(len(table[0]))]
data=pd.DataFrame(data)
data=data[['col','hadron','obs','RS','xpmin','xpmax','value','stat_u','syst_u']]
writer = pd.ExcelWriter('ALEPH.xlsx')
data.to_excel(writer,'ALEPH.xlsx',index=False)
writer.save()


#numbers=[]
#for i in range(len(L)): 
#  numbers.extend(L[i].replace('{',' ').split())
#numbers=[float(x)  for x in numbers]
#
#table=[]
#row=[]
#cnt=0
#while 1:
#  if len(numbers)==0: break
#  row.append(numbers.pop(0))
#  cnt+=1
#  if cnt==10: 
#    table.append(row)
#    row=[]
#    cnt=0  
#
#for i in range(len(table)):
#  table[i]=table[i][-5:]
#
#table=np.transpose(table)
#data={}
#data['xpmin']=table[0]
#data['xpmax']=table[1]
#data['value']=table[2]
#data['error_u']=table[3]
#data['hadron']=['hadron' for _ in range(len(table[0]))]
#data['col']=['DELPHI' for _ in range(len(table[0]))]
#data['RS']=[91.2 for _ in range(len(table[0]))]
#data=pd.DataFrame(data)
#writer = pd.ExcelWriter('DELPHI.xlsx')
#data.to_excel(writer,'DELPHI.xlsx',index=False)
#writer.save()










