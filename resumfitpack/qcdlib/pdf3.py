#!/bin/env python
import sys
import numpy as np
from scipy.special   import gamma
from qcdlib.alphaS   import ALPHAS
from qcdlib.aux      import AUX
from qcdlib.dglap    import DGLAP
from qcdlib.kernels  import KERNELS
from qcdlib.mellin   import MELLIN
from tools.config    import conf
from scipy.integrate import quad
      
class PDF: 
  
    def __init__(self,mellin=None):
  
        self.spl='upol'
        self.Q20=conf['Q20']
        self.mc2=conf['aux'].mc2
        self.mb2=conf['aux'].mb2
  
        if mellin==None:
            self.kernel=KERNELS(conf['mellin'],self.spl)
            self.dglap=DGLAP(conf['mellin'],conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
            self.mellin=conf['mellin']
        else:
            self.kernel=KERNELS(mellin,self.spl)
            self.dglap=DGLAP(mellin,conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
            self.mellin=mellin
  
        if 'su2+su3' not in conf: conf['su2+su3']=False
        self.set_params()
        self.setup()
  
    def set_params(self):
        """
        f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        """
        self.params={}
  
        self.params['g+1']    =np.array([1,-0.5,3,0,0])
        self.params['uv+1']   =np.array([1, 0.5,3,0,0])
        self.params['dv+1']   =np.array([1, 0.5,4,0,0])
        self.params['sea1+1'] =np.array([1,-0.5,6,0,0])
        self.params['sea2+1'] =np.array([1,-0.5,6,0,0])
        self.params['db+1']   =np.array([1,-0.5,6,0,0])
        self.params['ub+1']   =np.array([1,-0.5,6,0,0])
        self.params['s+1']    =np.array([1,-0.5,6,0,0])
        self.params['sb+1']   =np.array([1,-0.5,6,0,0])
  
        self.params['g-1']    =np.array([1,-0.5,3,0,0])
        self.params['uv-1']   =np.array([1, 0.5,3,0,0])
        self.params['dv-1']   =np.array([1, 0.5,4,0,0])
        self.params['sea1-1'] =np.array([1,-0.5,6,0,0])
        self.params['sea2-1'] =np.array([1,-0.5,6,0,0])
        self.params['db-1']   =np.array([1,-0.5,6,0,0])
        self.params['ub-1']   =np.array([1,-0.5,6,0,0])
        self.params['s-1']    =np.array([1,-0.5,6,0,0])
        self.params['sb-1']   =np.array([1,-0.5,6,0,0])
  
    def set_sumrules1(self):
        m=lambda _,n: self.get_moments(_,n) 
  
        #--valence
        self.params['uv+1'][0]=1    
        self.params['uv+1'][0]=(2-m('uv-1',1) )/m('uv+1',1)
  
        self.params['dv+1'][0]=1    
        self.params['dv+1'][0]=(1-m('dv-1',1) )/m('dv+1',1)
  
        #--strange
        self.params['s+1'][0]=1    
        self.params['s+1'][0]=(m('sb+1',1)+m('sb-1',1)-m('s-1',1))/m('s+1',1)
    
        #--msr
        sea1=m('sea1+1',2)+m('sea1-1',2)
        sea2=m('sea2+1',2)+m('sea2-1',2)
        uv=m('uv+1',2)+m('uv-1',2)
        dv=m('dv+1',2)+m('dv-1',2)
        ub=m('ub+1',2)+m('ub-1',2)+sea1
        db=m('db+1',2)+m('db-1',2)+sea1
        s =m('s+1' ,2)+m('s-1' ,2)+sea2
        sb=m('sb+1',2)+m('sb-1',2)+sea2
  
  
        up=uv+2*ub
        dp=dv+2*db
        sp=s+sb
        self.params['g+1'][0]=1    
        self.params['g+1'][0]=(1-up-dp-sp-m('g-1',2))/m('g+1',2)
        g=m('g+1',2)+m('g-1',2)
        msr=g+up+dp+sp
  
        #--share
        self.sr={}
        self.sr['msr']      = msr
        self.sr['uv(1)']    = m('uv+1',1)+m('uv-1',1)
        self.sr['dv(1)']    = m('dv+1',1)+m('dv-1',1)
        self.sr['s-sb(1)']  = m('s+1' ,1)+m('s-1' ,1) - m('sb+1' ,1) -m('sb-1' ,1)
        self.sr['s-sb(2)']  = m('s+1' ,2)+m('s-1' ,2) - m('sb+1' ,2) -m('sb-1' ,2)
        self.sr['db-ub(1)'] = m('db+1',1)+m('db-1' ,1)- m('ub+1' ,1) -m('ub-1' ,1)
        self.sr['db-ub(2)'] = m('db+1',2)+m('db-1' ,2)- m('ub+1' ,2) -m('ub-1' ,2)
  
    def set_sumrules2(self):
        m=lambda _,n: self.get_moments(_,n) 
  
        sea1_p =m('sea1+1',1)
        sea1_m =m('sea1-1',1)
        sea2_p =m('sea2+1',1)
        sea2_m =m('sea2-1',1)
  
        ub_p=m('ub+1',1) + sea1_p
        ub_m=m('ub-1',1) + sea1_m
        db_p=m('db+1',1) + sea1_p
        db_m=m('db-1',1) + sea1_m
        s_p =m('s+1',1)  + sea2_p
        s_m =m('s-1',1)  + sea2_m
        sb_p=m('sb+1',1) + sea2_p
        sb_m=m('sb-1',1) + sea2_m
  
  
        #--impose zero net strangens in the nucleon      
        self.params['s+1'][0]=1    
        self.params['s+1'][0]=(sb_p+sb_m-s_m)/m('s+1',1)
        s_p =m('s+1',1) + sea2_p
  
  
        #--buld delta s+
        Ds  = s_p  - s_m
        Dsb = sb_p - sb_m
        Dsp = Ds + Dsb
    
        #--impose valence sum rule and g3, g8
        g3=1.269
        g8=0.586
  
        self.params['uv+1'][0]=1    
        self.params['uv-1'][0]=1    
        self.params['uv+1'][0]=(+Dsp/2-ub_p+ub_m+(g3+g8)/4+1)/m('uv+1',1)
        self.params['uv-1'][0]=(-Dsp/2+ub_p-ub_m-(g3+g8)/4+1)/m('uv-1',1)
  
        self.params['dv+1'][0]=1    
        self.params['dv-1'][0]=1    
        self.params['dv+1'][0]=(+2*Dsp-4*db_p+4*db_m-g3+g8+2)/m('dv+1',1)/4
        self.params['dv-1'][0]=(-2*Dsp+4*db_p-4*db_m+g3-g8+2)/m('dv-1',1)/4
  
  
    
        #--constrains from second moments
        sea1=m('sea1+1',2)+m('sea1-1',2)
        sea2=m('sea2+1',2)+m('sea2-1',2)
        uv=m('uv+1',2)+m('uv-1',2)
        dv=m('dv+1',2)+m('dv-1',2)
        ub=m('ub+1',2)+m('ub-1',2)+sea1
        db=m('db+1',2)+m('db-1',2)+sea1
        s =m('s+1' ,2)+m('s-1' ,2)+sea2
        sb=m('sb+1',2)+m('sb-1',2)+sea2
  
        up=uv+2*ub
        dp=dv+2*db
        sp=s+sb
        self.params['g+1'][0]=1    
        self.params['g+1'][0]=(1-up-dp-sp-m('g-1',2))/m('g+1',2)
        g=m('g+1',2)+m('g-1',2)
        msr=g+up+dp+sp
  
        #--share
        self.sr={}
        self.sr['msr']      = msr
        self.sr['uv(1)']    = m('uv+1',1)+m('uv-1',1)
        self.sr['dv(1)']    = m('dv+1',1)+m('dv-1',1)
        self.sr['s-sb(1)']  = m('s+1' ,1)+m('s-1' ,1) - m('sb+1' ,1) -m('sb-1' ,1)
        self.sr['s-sb(2)']  = m('s+1' ,2)+m('s-1' ,2) - m('sb+1' ,2) -m('sb-1' ,2)
        self.sr['db-ub(1)'] = m('db+1',1)+m('db-1' ,1)- m('ub+1' ,1) -m('ub-1' ,1)
        self.sr['db-ub(2)'] = m('db+1',2)+m('db-1' ,2)- m('ub+1' ,2) -m('ub-1' ,2)
  
    def set_moms(self):
  
        m=lambda _: self.get_moments(_) 
  
        sea1=m('sea1+1')+m('sea1-1')
        sea2=m('sea2+1')+m('sea2-1')
        uv=m('uv+1')+m('uv-1')
        dv=m('dv+1')+m('dv-1')
        ub=m('ub+1')+m('ub-1')+sea1
        db=m('db+1')+m('db-1')+sea1
        s =m('s+1' )+m('s-1' )+sea2
        sb=m('sb+1')+m('sb-1')+sea2
        g =m('g+1') +m('g-1')
  
  
        moms={}
        moms['g']  = g
        moms['up'] = uv+2*ub
        moms['dp'] = dv+2*db
        moms['sp'] = s+sb
        moms['um'] = uv 
        moms['dm'] = dv 
        moms['sm'] = s-sb 
        self.moms0=moms
        self.get_BC(moms)
  
    def setup(self):
  
        if conf['su2+su3']==True: self.set_sumrules2()
        else: self.set_sumrules1()
        self.set_moms()
  
        #--store moments of a given Q2 that has been already calculated
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
        norm=self.beta(2+a,b+1)+c*self.beta(2+a+0.5,b+1)+d*self.beta(2+a+1.0,b+1)
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
                self.storage[Q2]=self.dglap.evolve(self.BC3,self.Q20,Q2,3)
    
    def get_xF(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return x*conf['mellin'].invert(x,self.storage[Q2][flav])
  
    def get_xF0(self,x,flav,evolve=True):
        if   flav=='um': mom=self.moms0['um']
        elif flav=='dm': mom=self.moms0['dm']
        elif flav=='sm': mom=self.moms0['sm']
        return x*conf['mellin'].invert(x,mom)
  
if __name__=='__main__':

    from scipy.integrate import quad
    conf['order']='NLO'
    conf['Q20'] = 1.27**2
    conf['dglap mode']='truncated'
    conf['aux']=AUX()
    conf['mellin']=MELLIN(npts=16)
    conf['alphaS']=ALPHAS()
    conf['su2+su3']=True
    pdf=PDF()
    print pdf.get_xF(0.5,10.,'u')




