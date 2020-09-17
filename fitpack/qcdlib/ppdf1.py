#!/bin/env python
from scipy.special import gamma
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from tools.config import conf
      
class PPDF: 

    def __init__(self, mellin=None):
  
        self.spl='pol'
        self.Q20=conf['Q20']
        self.mc2=conf['aux'].mc2
        self.mb2=conf['aux'].mb2
        if mellin==None:
            self.kernel=KERNELS(conf['mellin'],self.spl)
            if 'mode' in conf: mode=conf['mode']
            else: mode='truncated'
            self.dglap=DGLAP(conf['mellin'],conf['alphaS'],self.kernel,mode,conf['order'])
            self.mellin=conf['mellin']
        else:
            self.kernel=KERNELS(mellin,self.spl)
            if 'mode' in conf: mode=conf['mode']
            else: mode='truncated'
            self.dglap=DGLAP(mellin,conf['alphaS'],self.kernel,mode,conf['order'])
            self.mellin=mellin

        self.fub=1
        self.fdb=1

        if 'fub' in conf: self.fub=conf['fub']
        if 'fdb' in conf: self.fdb=conf['fdb']
        
        self.set_params()
        self.setup()
  
    def set_params(self):
        """
        f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        """
        self.params={}
        self.params['g1']   = np.array([1,-0.5,3,0,0])
        self.params['uv1']  = np.array([1,-0.5,3,0,0])
        self.params['dv1']  = np.array([1,-0.5,3,0,0])
        self.params['ub1']  = np.array([1, 0.5,3,0,0])
        self.params['sea1'] = np.array([1,-0.5,3,0,0])
        self.params['sea2'] = np.array([1,-0.5,3,0,0])
        self.params['db1']  = np.array([1, 0.5,3,0,0])
        self.params['ub1']  = np.array([0,-0.5,3,0,0])
        self.params['s1']   = np.array([1,-0.5,3,0,0])
        self.params['sb1']  = np.array([1,-0.5,3,0,0])
  
    def set_sumrules(self):

        gA=1.269
        g8=0.586

        self.params['uv1'][0]=1    
        self.params['dv1'][0]=1    
  
        ub=self.get_moments('ub1',1)+self.get_moments('sea1',1)
        db=self.get_moments('db1',1)+self.get_moments('sea1',1)
        uv=self.get_moments('uv1',1)
        dv=self.get_moments('dv1',1)
        s =self.get_moments('s1',1) +self.get_moments('sea2',1)
        sb=self.get_moments('sb1',1)+self.get_moments('sea2',1)
        sp=s+sb  

        ub*=self.fub
        db*=self.fdb

        self.params['uv1'][0]=(+gA + g8 -4*ub+2*sp)/2.0/uv
        self.params['dv1'][0]=(-gA + g8 -4*db+2*sp)/2.0/dv
  
        gA= (self.params['uv1'][0]+2*self.params['ub1'][0])\
           -(self.params['dv1'][0]+2*self.params['db1'][0])
        print 'gA=',gA

    def setup(self):
    
        if 'su2+su3' in conf:
            if conf['su2+su3']==True: self.set_sumrules()
  
        moms={}

        ub=self.get_moments('ub1')+self.get_moments('sea1')
        db=self.get_moments('db1')+self.get_moments('sea1')
        uv=self.get_moments('uv1')
        dv=self.get_moments('dv1')
        s =self.get_moments('s1') +self.get_moments('sea2')
        sb=self.get_moments('sb1')+self.get_moments('sea2')
        sp=s+sb  

        ub*=self.fub
        db*=self.fdb
        
        moms['g']=self.get_moments('g1')

        moms['up']=uv+2*ub
        moms['dp']=dv+2*db
        moms['sp']=sp
        moms['um']=uv
        moms['dm']=dv
        moms['sm']=s-sb
  
        self.moms0=moms
        self.get_BC(moms)
  
        #--we will store all Q2 values that has been precalc
        self.storage={}
  
    def beta(self,a,b):
        return gamma(a)*gamma(b)/gamma(a+b)
  
    def get_moments(self,flav,N=None):
        """
        if N==None: then parametrization is to be use to compute moments along mellin contour
        else the Nth moment is returned
        """
        if N==None: N=self.mellin.N
        M,a,b,c,d=self.params[flav]
        norm=self.beta(1+a,b+1)+c*self.beta(1+a+0.5,b+1)+d*self.beta(1+a+1.0,b+1)
        mom=self.beta(N+a,b+1)+c*self.beta(N+a+0.5,b+1)+d*self.beta(N+a+1.0,b+1)
        return M*mom/norm 
  
    def _get_BC(self,g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm):
        N=self.mellin.N
        
        # flav composition
        vm,vp={},{}
        vm[35]= bm + cm + dm + sm - 5*tm + um
        vm[24]= -4*bm + cm + dm + sm + um 
        vm[15]= -3*cm + dm + sm + um
        vm[8] = dm - 2*sp + 2*(-sm + sp) + um
        vm[3] = -dm + um
        vm[0] = np.zeros(N.size,dtype=complex)
        vp[0] = np.zeros(N.size,dtype=complex)
        vp[3] = -dp + up
        vp[8] = dp - 2*sp + up 
        vp[15]= -3*cp + dp + sp + up 
        vp[24]= -4*bp + cp + dp + sp + up 
        vp[35]= bp + cp + dp + sp - 5*tp + up 
        qs    = bp + cp + dp + sp + tp + up 
        qv    = bm + cm + dm + sm + tm + um 
        q     = np.zeros((2,N.size),dtype=complex)
        q[0]=np.copy(qs)
        q[1]=np.copy(g)
  
        BC={}
        BC['vm']=vm 
        BC['vp']=vp 
        BC['qv']=qv
        BC['q'] =q
        return BC
  
    def get_state(self):
        return (self.BC3,self.BC4,self.BC5)
    
    def set_state(self,state):
        self.BC3, self.BC4, self.BC5 = state[:]
        self.storage = {}
    
    def get_BC(self,moms):
  
        N=self.mellin.N
        zero=np.zeros(N.size,dtype=complex)
  
        ###############################################
        # BC for Nf=3
        g   = moms['g']
        up  = moms['up']
        um  = moms['um']
        dp  = moms['dp']
        dm  = moms['dm']
        sp  = moms['sp']
        sm  = moms['sm']
        cp  = zero
        cm  = zero
        bp  = zero
        bm  = zero
        self.BC3=self._get_BC(g,up,um,dp,dm,sp,sm,zero,zero,zero,zero,zero,zero)
  
        ###############################################
        # BC for Nf=4
        BC4=self.dglap.evolve(self.BC3,self.Q20,self.mc2,3)
        g =BC4['g']
        up=BC4['up']
        dp=BC4['dp']
        sp=BC4['sp']
        cp=BC4['cp']
        bp=BC4['bp']
        tp=BC4['tp']
        um=BC4['um']
        dm=BC4['dm']
        sm=BC4['sm']
        cm=BC4['cm']
        bm=BC4['bm']
        tm=BC4['tm']
        self.BC4=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm)
  
        ###############################################
        # BC for Nf=5
        BC5=self.dglap.evolve(self.BC4,self.mc2,self.mb2,4)
        g =BC5['g']
        up=BC5['up']
        dp=BC5['dp']
        sp=BC5['sp']
        cp=BC5['cp']
        bp=BC5['bp']
        tp=BC5['tp']
        um=BC5['um']
        dm=BC5['dm']
        sm=BC5['sm']
        cm=BC5['cm']
        bm=BC5['bm']
        tm=BC5['tm']
        self.BC5=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm)
  
    def evolve(self,Q2):
  
        if Q2 not in self.storage:
            if self.mb2<Q2: 
                self.storage[Q2]=self.dglap.evolve(self.BC5,self.mb2,Q2,5)
            elif self.mc2<=Q2 and Q2<=self.mb2: 
                self.storage[Q2]=self.dglap.evolve(self.BC4,self.mc2,Q2,4)
            elif Q2<self.mc2: 
                self.storage[Q2] =self.dglap.evolve(self.BC3,self.Q20,Q2,3)
    
    def get_xF(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return x*self.mellin.invert(x,self.storage[Q2][flav])
  
    def get_xF0(self,x,flav,evolve=True):
        if flav=='um': mom=self.moms0['um']
        if flav=='dm': mom=self.moms0['dm']
        return x*self.mellin.invert(x,mom)
  
if __name__=='__main__':

    conf['order']='NLO'
    conf['Q20'] = 1.0
    conf['aux']=AUX()
    conf['mellin']=MELLIN(npts=4)
    conf['alphaS']=ALPHAS()
    conf['su2+su3']=True
    ppdf=PPDF()
    #print ppdf.get_xF(0.5,10.,'u')
    dga = lambda x:  ppdf.get_xF(x,2.0,'u')/x+ppdf.get_xF(x,2.0,'ub')/x\
                    -ppdf.get_xF(x,2.0,'d')/x-ppdf.get_xF(x,2.0,'db')/x
    from scipy.integrate import quad 
    gA=quad(dga, 0, 1)
    print gA

