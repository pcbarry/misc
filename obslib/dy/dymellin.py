#!/usr/bin/env python
import sys,os
import numpy as np
import time

#--from scipy
from scipy.integrate import fixed_quad
from scipy.special import gamma

#--from tools
from tools.bar import BAR
from tools.multiproc import MULTIPROC
from tools.tools import load, save, checkdir
from tools.config import conf
from tools.parallel import PARALLEL 

#--from qcdlib
from qcdlib import aux,mellin,alphaS,eweak,pdf0

#--local
from fakepdf  import FAKEPDF
from theory   import DY
from reader   import READER

class DYMELLIN(DY):
  
    def gen_SIGNM(self,N,M,Q2,S,Y,muF2,part,msg='calculating...'):
        pts=len(N)
        bar=BAR(msg,pts**2)
        SIGNM=np.zeros((pts,pts),dtype=complex)
        for i in range(pts): 
            for j in range(pts): 
                self.n=N[i]
                self.m=M[j]
                real=self.get_xsec(Q2,S,Y,muF2,ilum='mell-real',part=part)
                imag=self.get_xsec(Q2,S,Y,muF2,ilum='mell-imag',part=part)
                SIGNM[i,j]=np.complex(real,imag)
                bar.next()
        bar.finish()
        return SIGNM
  
    def _gen_melltab(self,Q2,S,Y,muF2):
        N =self.mellin.N
        NS=np.conjugate(N)
        data={}
        for part in ['qA,qbB','qbA,qB','qA,gB','gA,qB']:
            data[part]={}
            data[part]['SIGNM']=self.gen_SIGNM(N ,N,Q2,S,Y,muF2,part,msg='calculating %s N M'%part)
            data[part]['SIGNCM']=self.gen_SIGNM(NS,N,Q2,S,Y,muF2,part,msg='calculating %s N*M'%part)
        return data

    def gen_melltab(self):

        #--construct requests
        requests=[]
        for k in conf['dy tabs']:
            path2dytabK='%s/%d'%(conf['path2dytab'],k)
            checkdir(path2dytabK)
            npts=len(conf['dy tabs'][k]['value'])
            for i in range(npts):
                data={}
                data['S']    = conf['dy tabs'][k]['S'][i]
                data['Y']    = conf['dy tabs'][k]['Y'][i]
                data['Q2']   = conf['dy tabs'][k]['Q2'][i]
                data['idx']  = conf['dy tabs'][k]['idx'][i]
                data['muF2'] = data['Q2']
                data['path2dytabK'] = path2dytabK
                requests.append(data)
                    
        def task(data):
            Q2  = data['Q2']
            S   = data['S']
            Y   = data['Y']
            idx = data['idx']
            muF2= data['muF2']
            mtab=self._gen_melltab(Q2,S,Y,muF2)
            fname='%s/%d.melltab'%(data['path2dytabK'],idx)
            save(mtab,fname)
    
        #--now we setup the parallel class
        nworkers=20
        parallel=PARALLEL()
        parallel.task=task
        parallel.setup_master()
        parallel.setup_workers(nworkers)
        parallel.send_tasks(requests)
        parallel.stop_workers()

    def load_melltab(self):
        self.mtab={}
        path2dytab=conf['path2dytab']
        for k in conf['dy tabs']:
            path2dytabK='%s/%d'%(path2dytab,k)
            self.mtab[k]={}
            npts=len(conf['dy tabs'][k]['value'])
            bar=BAR('loading DY tables of %d'%k,npts)
            for i in range(npts):
                idx= conf['dy tabs'][k]['idx'][i]
                fname='%s/%d.melltab'%(path2dytabK,idx)
                self.mtab[k][i]=load(fname)
                bar.next()
            bar.finish() 
  
    def invert(self,fA,fB,sigNM,sigNCM):
        fAS=np.conjugate(fA)
        W=self.mellin.W*self.mellin.JAC
        phase2=self.mellin.phase**2
        xsec1=np.einsum('i,j,i,j,ij',W,W,fA,fB,phase2*sigNM,optimize=True)
        xsec2=np.einsum('i,j,i,j,ij',W,W,fAS,fB,sigNCM,optimize=True)
        return np.real((-2/4./np.pi**2)*(xsec1-xsec2))
  
    def get_mxsec(self,k,i,muF2,reaction):
     
        conf['pdfA'].evolve(muF2)
        conf['pdfB'].evolve(muF2)
  
        gA  = conf['pdfA'].storage[muF2]['g'] 
        uA  = conf['pdfA'].storage[muF2]['u'] 
        dA  = conf['pdfA'].storage[muF2]['d'] 
        sA  = conf['pdfA'].storage[muF2]['s'] 
        cA  = conf['pdfA'].storage[muF2]['c'] 
        bA  = conf['pdfA'].storage[muF2]['b'] 
        ubA = conf['pdfA'].storage[muF2]['ub']
        dbA = conf['pdfA'].storage[muF2]['db']
        sbA = conf['pdfA'].storage[muF2]['sb']
        cbA = conf['pdfA'].storage[muF2]['cb']
        bbA = conf['pdfA'].storage[muF2]['bb']
  
        gB  = conf['pdfB'].storage[muF2]['g'] 
        uB  = conf['pdfB'].storage[muF2]['u'] 
        dB  = conf['pdfB'].storage[muF2]['d'] 
        sB  = conf['pdfB'].storage[muF2]['s'] 
        cB  = conf['pdfB'].storage[muF2]['c'] 
        bB  = conf['pdfB'].storage[muF2]['b'] 
        ubB = conf['pdfB'].storage[muF2]['ub']
        dbB = conf['pdfB'].storage[muF2]['db']
        sbB = conf['pdfB'].storage[muF2]['sb']
        cbB = conf['pdfB'].storage[muF2]['cb'] 
        bbB = conf['pdfB'].storage[muF2]['bb'] 
        
  
        if  reaction=='pd':
            uB=0.5*(uB+dB)
            dB=uB
            ubB=0.5*(ubB+dbB)
            dbB=ubB
        
        xsec = self.eU2*self.invert(uA  ,ubB ,self.mtab[k][i]['qA,qbB']['SIGNM'],self.mtab[k][i]['qA,qbB']['SIGNCM'])\
              +self.eD2*self.invert(dA  ,dbB ,self.mtab[k][i]['qA,qbB']['SIGNM'],self.mtab[k][i]['qA,qbB']['SIGNCM'])\
              +self.eD2*self.invert(sA  ,sbB ,self.mtab[k][i]['qA,qbB']['SIGNM'],self.mtab[k][i]['qA,qbB']['SIGNCM'])\
              +self.eU2*self.invert(ubA ,uB  ,self.mtab[k][i]['qbA,qB']['SIGNM'],self.mtab[k][i]['qbA,qB']['SIGNCM'])\
              +self.eD2*self.invert(dbA ,dB  ,self.mtab[k][i]['qbA,qB']['SIGNM'],self.mtab[k][i]['qbA,qB']['SIGNCM'])\
              +self.eD2*self.invert(sbA ,sB  ,self.mtab[k][i]['qbA,qB']['SIGNM'],self.mtab[k][i]['qbA,qB']['SIGNCM'])\
              +self.eU2*self.invert(uA  ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eD2*self.invert(dA  ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eD2*self.invert(sA  ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eU2*self.invert(ubA ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eD2*self.invert(dbA ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eD2*self.invert(sbA  ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
              +self.eU2*self.invert(gA  ,uB  ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
              +self.eD2*self.invert(gA  ,dB  ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
              +self.eD2*self.invert(gA  ,sB  ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
              +self.eU2*self.invert(gA  ,ubB ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
              +self.eD2*self.invert(gA  ,dbB ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
              +self.eD2*self.invert(gA  ,sbB ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )
  
  
        Nf=conf['alphaS'].get_Nf(muF2)
        if Nf>3: 
          xsec+= self.eU2*self.invert(cA  ,cbB  ,self.mtab[k][i]['qA,qbB']['SIGNM'],self.mtab[k][i]['qA,qbB']['SIGNCM'])\
                +self.eU2*self.invert(cbA ,cB   ,self.mtab[k][i]['qbA,qB']['SIGNM'],self.mtab[k][i]['qbA,qB']['SIGNCM'])\
                +self.eU2*self.invert(cA  ,gB   ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
                +self.eU2*self.invert(cbA ,gB   ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
                +self.eU2*self.invert(gA  ,cB   ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
                +self.eU2*self.invert(gA  ,cbB  ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )
                                                          
        if Nf>4:                                          
          xsec+= self.eD2*self.invert(bA  ,bbB ,self.mtab[k][i]['qA,qbB']['SIGNM'],self.mtab[k][i]['qA,qbB']['SIGNCM'])\
                +self.eD2*self.invert(bbA ,bB  ,self.mtab[k][i]['qbA,qB']['SIGNM'],self.mtab[k][i]['qbA,qB']['SIGNCM'])\
                +self.eD2*self.invert(bA  ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
                +self.eD2*self.invert(bbA ,gB  ,self.mtab[k][i]['qA,gB']['SIGNM'] ,self.mtab[k][i]['qA,gB']['SIGNCM'] )\
                +self.eD2*self.invert(gA  ,bB  ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )\
                +self.eD2*self.invert(gA  ,bbB ,self.mtab[k][i]['gA,qB']['SIGNM'] ,self.mtab[k][i]['gA,qB']['SIGNCM'] )
        return xsec

    def test(self):

        self.load_melltab()
        t=time.time()
        for k in conf['dy tabs']:
            npts=len(conf['dy tabs'][k]['value'])
            for i in range(npts):
                Q2 = conf['dy tabs'][k]['Q2'][i]
                reaction = conf['dy tabs'][k]['reaction'][i]
                approx = self.get_mxsec(k,i,Q2,reaction)
                S  = conf['dy tabs'][k]['S'][i]
                Y  = conf['dy tabs'][k]['Y'][i]
                muF2=Q2
                exact = self.get_xsec(Q2,S,Y,muF2,ilum='normal',part='full')
                rel_err=abs((approx-exact)/exact)*100
                print 'rap=%5.2f  Q2=%5.2f  rel-err=%5.2f'%(Y,Q2,rel_err)
        print time.time()-t
  
if __name__=='__main__':

    conf['Q20']   = 1.0
    conf['alphaSmode']='backward'
    conf['dglap mode']='truncated'
    conf['order']='NLO'
    conf['scheme']='ZMVFS'
    conf['path2dytab']='%s/grids/grids-dy'%os.environ['FITPACK']
    conf['aux']=aux.AUX()
    conf['mellin']=mellin.MELLIN(npts=4)
    conf['alphaS']=alphaS.ALPHAS()
    conf['eweak']=eweak.EWEAK()
    conf['pdfA']=FAKEPDF()
    conf['pdfB']=FAKEPDF()

    conf['aux']=aux.AUX()
    conf['datasets']={}
    conf['datasets']['dy']={}
    conf['datasets']['dy']['xlsx']={}
    conf['datasets']['dy']['xlsx'][10001]='dy/expdata/10001.xlsx'
    conf['datasets']['dy']['xlsx'][10002]='dy/expdata/10002.xlsx'
    conf['datasets']['dy']['norm']={}
    conf['datasets']['dy']['filters']=[]
    conf['datasets']['dy']['filters'].append("Q2>1") 
    conf['dy tabs']=READER().load_data_sets('dy')

    dymellin=DYMELLIN()
    #dymellin.gen_melltab()
    dymellin.test()





