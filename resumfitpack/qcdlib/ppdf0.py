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
        
        self.set_params()
        self.setup()

    def set_params(self):
        """
        f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        """
        self.params={}
        self.params['g1'] =np.array([1,-0.5,3,0,0])
        self.params['up1']=np.array([1, 0.5,3,0,0])
        self.params['dp1']=np.array([1, 0.5,3,0,0])
        self.params['sp1']=np.array([1,-0.5,3,0,0])
        self.params['ub1']=np.array([1,-0.5,3,0,0])
        self.params['db1']=np.array([1,-0.5,3,0,0])
        self.params['sb1']=np.array([0,-0.5,3,0,0])

    def set_sumrules(self):
        gA=1.269
        g8=0.586

        self.params['up1'][0]=1    
        self.params['up1'][0]=(gA + g8)/2. + self.get_moments('sp1',1)

        self.params['dp1'][0]=1    
        self.params['dp1'][0]=(g8 - gA)/2. + self.get_moments('sp1',1)

    def setup(self):
    
        if conf['su2+su3']: self.set_sumrules()

        moms={}

        moms['g']=self.get_moments('g1')
        moms['up']=self.get_moments('up1')
        moms['dp']=self.get_moments('dp1')
        moms['sp']=self.get_moments('sp1')

        ub=db=sb=self.get_moments('sp1')/2

        moms['um']=moms['up']-2*ub
        moms['dm']=moms['dp']-2*db
        moms['sm']=moms['sp']-2*sb

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

if  __name__=='__main__':

    conf['order']='NLO'
    conf['Q20'] = 1.0
    conf['aux']=AUX()
    conf['mellin']=MELLIN(npts=4)
    conf['alphaS']=ALPHAS()
    conf['su2+su3'] = True
    ppdf=PPDF()
    print ppdf.get_xF(0.5,10.,'u')
    print ppdf.storage[10.0]['g']
    print len(ppdf.storage[10.0]['g'])

