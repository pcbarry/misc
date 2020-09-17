#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys,os
import numpy as np
import pylab as py
import pandas as pd

L=open("SLD").readlines()

L=[l.strip() for l in L]
L=[l for l in L if l!='']
L=[l.replace('?',' ') for l in L]
L=[l.replace(r'–',' ') for l in L]
L=[[float(_) for _ in l.split()] for l in L]

numbers=[]
for i in range(len(L)): 
  numbers.extend(L[i])

table=[]
row=[]
cnt=0
while 1:
  if len(numbers)==0: break
  row.append(numbers.pop(0))
  cnt+=1
  if cnt==13: 
    table.append(row)
    row=[]
    cnt=0  

#for l in table: print l

table2=[]
for i in range(len(table)):
  row=[]
  row.append(np.exp(-table[i][3]))
  row.append(np.exp(-table[i][4]))
  row.append(table[i][9])
  row.append(table[i][10])
  table2.append(row)

table=np.transpose(table2)
data={}
data['xpmin']=table[0]
data['xpmax']=table[1]
data['value']=table[2]
data['error_u']=table[3]
data['hadron']=['hadron' for _ in range(len(table[0]))]
data['col']=['SLD' for _ in range(len(table[0]))]
data['RS']=[91.28 for _ in range(len(table[0]))]
data['obs']=["1/sig dsig/dxp" for _ in range(len(table[0]))]
data=pd.DataFrame(data)
data=data[['col','hadron','obs','RS','xpmin','xpmax','value','error_u']]
writer = pd.ExcelWriter('SLD.xlsx')
data.to_excel(writer,'SLD.xlsx',index=False)
writer.save()










