#!/bin/env python
import sys
from scipy.special import gamma
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from tools.config import conf
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

        self.set_params()
        self.setup()

    def set_params(self):
        """
        f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        """
        params={}

        ## first shapes
        params['g1']    =np.array([1,-0.5,3,0,0])
        params['uv1']   =np.array([1, 0.5,3,0,0])
        params['dv1']   =np.array([1, 0.5,4,0,0])
        params['sea1']  =np.array([1,-1.19,4,0,0])
        params['sea2']  =np.array([1,-1.19,4,0,0])
        params['db1']   =np.array([1,-0.5,6,0,0])
        params['ub1']   =np.array([1,-0.5,6,0,0])
        params['s1']    =np.array([ 0.02117505,-0.7834729,5.94912366,0.,0.])
        params['sb1']   =np.array([ 0.03185304,-0.4761831,10.,0.,0.])

        ## second shapes
        ##                            n,    a,   b,   c,   d
        params['g2']    = np.array([0.0, -0.5, 3.0, 0.0, 0.0])
        params['uv2']   = np.array([0.0,  0.5, 3.0, 0.0, 0.0])
        params['dv2']   = np.array([0.0,  0.5, 4.0, 0.0, 0.0])
        params['db2']   = np.array([0.0, -0.5, 6.0, 0.0, 0.0])
        params['ub2']   = np.array([0.0, -0.5, 6.0, 0.0, 0.0])
        params['s2']    = np.array([0.0, -0.7834729, 5.94912366, 0.0, 0.0])
        params['sb2']   = np.array([0.0, -0.4761831, 10.0      , 0.0, 0.0])

        ## third shapes
        ##                            n,    a,   b,   c,   d
        params['g3']    = np.array([0.0, -0.5, 3.0, 0.0, 0.0])
        params['uv3']   = np.array([0.0,  0.5, 3.0, 0.0, 0.0])
        params['dv3']   = np.array([0.0,  0.5, 4.0, 0.0, 0.0])
        params['db3']   = np.array([0.0, -0.5, 6.0, 0.0, 0.0])
        params['ub3']   = np.array([0.0, -0.5, 6.0, 0.0, 0.0])
        params['s3']    = np.array([0.0, -0.7834729, 5.94912366, 0.0, 0.0])
        params['sb3']   = np.array([0.0, -0.4761831, 10.0      , 0.0, 0.0])

        self.params=params

    def set_sumrules(self):

        #--valence
        self.params['uv1'][0]=1
        # self.params['uv1'][0]=2/self.get_moments('uv1',1)
        self.params['uv1'][0]=(2-self.get_moments('uv2', 1)-self.get_moments('uv3', 1))/self.get_moments('uv1',1)

        self.params['dv1'][0]=1
        # self.params['dv1'][0]=1/self.get_moments('dv1',1)
        self.params['dv1'][0]=(1-self.get_moments('dv2', 1)-self.get_moments('dv3', 1))/self.get_moments('dv1',1)

        #--strange
        self.params['s1'][0]=1
        # self.params['s1'][0]=self.get_moments('sb1',1)/self.get_moments('s1',1)
        self.params['s1'][0]=(self.get_moments('sb1',1)+self.get_moments('sb2',1)+self.get_moments('sb3',1)-\
                              self.get_moments('s2',1)-self.get_moments('s3',1))/self.get_moments('s1',1)

        #--msr
        sea1=self.get_moments('sea1',2)
        sea2=self.get_moments('sea2',2)
        # up=self.get_moments('uv1',2)+2*(sea1+self.get_moments('ub1',2))
        up=self.get_moments('uv1',2)+self.get_moments('uv2',2)+self.get_moments('uv3',2)+\
           2*(sea1+self.get_moments('ub1',2)+self.get_moments('ub2',2)+self.get_moments('ub3',2))
        # dp=self.get_moments('dv1',2)+2*(sea1+self.get_moments('db1',2))
        dp=self.get_moments('dv1',2)+self.get_moments('dv2',2)+self.get_moments('dv3',2)+\
           2*(sea1+self.get_moments('db1',2)+self.get_moments('db2',2)+self.get_moments('db3',2))
        # sp=(sea2+self.get_moments('s1',2))+(sea2+self.get_moments('sb1',2))
        sp=(sea2+self.get_moments('s1',2)+self.get_moments('s2',2)+self.get_moments('s3',2))+\
           (sea2+self.get_moments('sb1',2)+self.get_moments('sb2',2)+self.get_moments('sb3',2))
        self.params['g1'][0]=1
        # self.params['g1'][0]=(1-up-dp-sp)/self.get_moments('g1',2)
        # g=self.get_moments('g1',2)
        self.params['g1'][0]=(1-up-dp-sp-self.get_moments('g2', 2)-self.get_moments('g3', 2))/self.get_moments('g1',2)
        g=self.get_moments('g1',2)+self.get_moments('g2', 2)+self.get_moments('g3', 2)
        msr=g+up+dp+sp

        #--share
        self.sr={}
        self.sr['msr']      = msr
        self.sr['uv(1)']    = self.get_moments('uv1',1)
        self.sr['dv(1)']    = self.get_moments('dv1',1)
        self.sr['s-sb(1)']  = self.get_moments('s1',1)-self.get_moments('sb1',1)
        self.sr['s-sb(2)']  = self.get_moments('s1',2)-self.get_moments('sb1',2)
        self.sr['db-ub(1)'] = self.get_moments('db1',1)-self.get_moments('ub1',1)
        self.sr['db-ub(2)'] = self.get_moments('db1',2)-self.get_moments('ub1',2)

        #for _ in sorted(self.sr): print _, self.sr[_]

    def set_moms(self):

        sea1=self.get_moments('sea1')
        sea2=self.get_moments('sea2')

        moms={}
        # moms['g']  = self.get_moments('g1')
        # moms['up'] = self.get_moments('uv1')+2*(sea1+self.get_moments('ub1'))
        # moms['dp'] = self.get_moments('dv1')+2*(sea1+self.get_moments('db1'))
        # moms['sp'] = 2*sea2+self.get_moments('s1')+self.get_moments('sb1')
        # moms['um'] = self.get_moments('uv1')
        # moms['dm'] = self.get_moments('dv1')
        # moms['sm'] = self.get_moments('s1')-self.get_moments('sb1')
        moms['g']  = self.get_moments('g1')+self.get_moments('g2')+self.get_moments('g3')
        moms['up'] = self.get_moments('uv1')+self.get_moments('uv2')+self.get_moments('uv3')+\
                     2*(sea1+self.get_moments('ub1')+self.get_moments('ub2')+self.get_moments('ub3'))
        moms['dp'] = self.get_moments('dv1')+self.get_moments('dv2')+self.get_moments('dv3')+\
                     2*(sea1+self.get_moments('db1')+self.get_moments('db2')+self.get_moments('db3'))
        moms['sp'] = 2*sea2+self.get_moments('s1')+self.get_moments('s2')+self.get_moments('s3')+\
                     self.get_moments('sb1')+self.get_moments('sb2')+self.get_moments('sb3')
        moms['um'] = self.get_moments('uv1')+self.get_moments('uv2')+self.get_moments('uv3')
        moms['dm'] = self.get_moments('dv1')+self.get_moments('dv2')+self.get_moments('dv3')
        moms['sm'] = self.get_moments('s1')+self.get_moments('s2')+self.get_moments('s3')-\
                     self.get_moments('sb1')-self.get_moments('sb2')-self.get_moments('sb3')
        self.moms0=moms
        self.get_BC(moms)

    def setup(self):
        self.set_sumrules()
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

    def get_qF_real(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return conf['mell(qpdf)-real'].invert(x,self.storage[Q2][flav])

    def get_qF_imag(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return conf['mell(qpdf)-imag'].invert(x,self.storage[Q2][flav])

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
    pdf=PDF()
    # print pdf.get_xF(0.5,10.,'u')
    q2 = 1.27 ** 2.0
    flavor = 's'
    xfs = []
    xs = np.linspace(0.001, 1.0, 20)
    for x in xs:
        xfs.append(pdf.get_xF(x, q2, flavor))
    print xfs

    #f=lambda x: pdf.get_xF(x,conf['Q20'],'s')/x-pdf.get_xF(x,conf['Q20'],'sb')/x
    ##f=lambda x: pdf.get_xF(x,conf['Q20'],'sb')/x
    ##f=lambda x: pdf.get_xF(x,conf['Q20'],'sm')/x
    #print quad(f,0,1)[0]

    #f=lambda x: pdf.get_xF0(x,'sm')/x
    #print quad(f,0,1)[0]
