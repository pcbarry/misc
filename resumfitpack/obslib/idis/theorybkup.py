import sys,os
import numpy as np
from scipy.integrate import quad
from tools.storage import STORAGE
from obslib.idis.aux import AUX 
from tools.config import conf
from tools.bar import BAR
from tools.tools import checkdir, load

class STFUNCS(AUX):
  
    def __init__(self):
  
        self.helium_calc=False
        self.mellin  =conf['mellin']
        if 'dmellin' in conf: self.dmellin =conf['dmellin']
  
        if 't4F2' in conf: self.t4F2=conf['t4F2']
        if 't4FL' in conf: self.t4FL=conf['t4FL']
        if 't4F3' in conf: self.t4F3=conf['t4F3']
        
        if 't4W2' in conf: self.t4W2=conf['t4W2']
        if 't4WL' in conf: self.t4WL=conf['t4WL']
        if 't4W3' in conf: self.t4W3=conf['t4W3']

        if 'F2off' in conf: self.F2off=conf['F2off']
        if 'FLoff' in conf: self.FLoff=conf['FLoff']
        if 'F3off' in conf: self.F3off=conf['F3off']
  
        if 'dsmf' in conf: self.dsmf=conf['dsmf']
        if 'hsmf' in conf: self.hsmf=conf['hsmf']
      
        #--need ht parametrization to determine how power corrections are calculated.  0 for additive (proton/neutron), 1 for additive (proton/deuteron),
        #--2 for multiplicative (proton/deuteron). 
        self.ht=2
        if 'ht' in conf and conf['ht']==True:  self.ht = conf['ht parametrization']
         
         
  
        self.setup()

        if 'path2idistab' in conf:
  
            self.path2distab=conf['path2idistab']
            checkdir(conf['path2idistab'])
            if 'offshell'  not in conf: print('Loading DIS tables...')
            elif 'offshell' in conf and conf['offshell']==False: print('Loading DIS tables...')
            else: print('Loading DIS onshell and offshell tables...')
            self.setup_tables()


        #else:

            #print('ERR: path2idistab not found in conf')
            #sys.exit()

        if conf['order']=='LO': self.order=0
        if conf['order']=='NLO': self.order=1
  
    #--unpolarized TMCs and nuclear smearing
  
    def get_TMC_GP(self,N,x,Q2,proj):
        """
        - we follow Eq.(4) of Brady et al., PRD84, 074008 (2011)
        - extra factor of x included for F2 and FL projections
        """
        M2=conf['aux'].M2
        mu2=M2/Q2
        xx=x*x
        rho=(1+4*mu2*xx)**0.5
        xi=2*x/(1+rho)              #kinematical factors
  
        if proj=='22':
  
            C1= (1+rho)**2/(4*rho**3)
            C2= 3*x*(rho**2-1)/(2*rho**4)
            C3= (rho**2-1)/(2*x*rho)                      
            h2=1/N*(-1+xi**(-N))                          
            g2=(N+xi*(1-N)-xi**(1-N))/N/(1-N)             
            return C1*xi**(-N+1) + C2*(h2 + C3*g2)        
  
        elif proj=='LL':
  
            C1= (1+rho)**2/(4*rho)
            return C1*xi**(-N+1) 
  
        elif proj=='L2':
  
            C2= x*(rho**2-1)/(rho**2)
            C3= (rho**2-1)/(2*x*rho)
            h2=1/N*(-1+xi**(-N))
            g2=(N+xi*(1-N)-xi**(1-N))/N/(1-N)
            return C2*(h2 + C3*g2)
  
        elif proj=='33':
  
            C1= (1+rho)/(2*rho**2)
            C2= (rho**2-1)/(2*rho**3)
            h3=1/(-N)*(1-xi**(-N))
            return C1*xi**(-N) + C2*h3 
  
    def get_TMC_AOT(self,N,x,Q2,proj):
        """
        - we follow Eq.(17) of "What does kinematical target mass sensitivity in DIS reveal about hadron structure?"
        - extra factor of x included for F2 and FL projections
        """
        M2=conf['aux'].M2
        mu2=M2/Q2
        xx=x*x
        rho=(1+4*mu2*xx)**0.5
        xi=2*x/(1+rho)
  
        if proj=='22':
  
            C1= (1+rho)/(2*rho**2)
            return C1*xi**(-N+1)    
  
        elif proj=='LL':
            C1= 1 
            return C1*xi**(-N+1) 
  
        elif proj=='L2':
  
            C1= (rho-1)/2
            return C1*xi**(-N+1) 
  
        elif proj=='33':
            C1= 1/rho
            return C1*xi**(-N)
  
    def get_smf(self,x,Q2,target,twist,tmc,offshell,conj,fXX,proj):
  
        eps_d   = -0.00222
        eps_He3 = -0.007718
        eps_H3  = -0.008482
        mN   = (0.93827231 + 0.93956563)/2
        mD   = 2*mN + eps_d
        mHe3 = (2*0.93827231 + 0.93956563) + eps_He3
        mH3  = (0.93827231 + 2*0.93956563) + eps_H3
  
        if target=='deuteron': 
            SMF=self.dsmf
            ymax = mD/mN
        elif target=='helium': 
            SMF=self.hsmf
            ymax = mHe3/mN
        elif target=='tritium': 
            #SMF=self.hsmf
            #ymax = mH3/mN
            print("ERROR: no tritium yet")
            sys.exit()

        #--need xi for offshell with TMC
        M2=conf['aux'].M2
        mu2=M2/Q2
        xx=x*x
        rho=(1+4*mu2*xx)**0.5
        xi=2*x/(1+rho)             

        if offshell==False:
 
            f=lambda y:SMF.get_fXX(fXX,'onshell',x,Q2,y)
  
            if twist==2 and proj!=None:    # add TMC
                if tmc=='GP':  func=lambda N,x: self.get_TMC_GP(N,x,Q2,proj)
                if tmc=='AOT': func=lambda N,x: self.get_TMC_AOT(N,x,Q2,proj)
              
            elif twist==2 and proj==None:  # no TMC
                if 'f2' in fXX or 'fL' in fXX: func=lambda N,x: x**(-N+1)
                else:                          func=lambda N,x: x**(-N)
  
            elif twist==4 and self.ht==0:
                func=lambda N,x: x**(-N)


        if offshell==True:
  
            f=lambda y:SMF.get_fXX(fXX,'offshell',x,Q2,y)
 
            if twist==2 and proj!=None:    # add TMC
                if tmc=='GP':
                    if 'f2' in fXX or 'fL' in fXX: func=lambda N,M,x: self.get_TMC_GP(N,x,Q2,proj)*xi**(-M+1)
                    else:                          func=lambda N,M,x: self.get_TMC_GP(N,x,Q2,proj)*xi**(-M)
                
                if tmc=='AOT':
                    if 'f2' in fXX or 'fL' in fXX: func=lambda N,M,x: self.get_TMC_AOT(N,x,Q2,proj)*xi**(-M+1)
                    else:                          func=lambda N,M,x: self.get_TMC_AOT(N,x,Q2,proj)*xi**(-M)
              
            elif twist==2 and proj==None:  # no TMC
                if 'f2' in fXX or 'fL' in fXX:  func=lambda N,M,x: x**(-N+1-M+1)
                else:                           func=lambda N,M,x: x**(-N-M)
  
            elif twist==4 and self.ht==0:
                func=lambda N,M,x: x**(-N-M)


        if offshell==False: 
            N=self.mellin.N
            fN=np.zeros(N.size,dtype=complex)
            #bar=BAR('generating grid target:%s twist:%s fXX:%s proj:%s'%(target,twist,fXX,proj),N.size)
            for i in range(N.size):
                intRE=lambda y: f(y)*func(N[i],x/y).real
                intIM=lambda y: f(y)*func(N[i],x/y).imag
                fN[i]=complex(quad(intRE,x,ymax,full_output=1)[0],quad(intIM,x,ymax,full_output=1)[0])
                #bar.next()
            #bar.finish()
            return fN
 
        if offshell==True: 
            N=self.dmellin.N
            if conj: M=np.conj(self.dmellin.M)
            else:    M=self.dmellin.M
            fNM=np.zeros((N.size,M.size),dtype=complex)
            #bar=BAR('generating offshell grid target:%s twist:%s fXX:%s proj:%s'%(target,twist,fXX,proj),N.size*M.size)
            for i in range(N.size):
                for j in range(M.size):
                    intRE=lambda y: f(y)*func(N[i],M[j],x/y).real
                    intIM=lambda y: f(y)*func(N[i],M[j],x/y).imag
                    fNM[i][j]=complex(quad(intRE,x,ymax,full_output=1)[0],quad(intIM,x,ymax,full_output=1)[0])
                    #bar.next()
            #bar.finish()
            return fNM
 
    #--twist level unpolarized structure functions in Mellin space
  
    def get_T2CFX(self,stf,nucleon,Q2,Nf=None,evolve=True):   #Twist-2 SFs without TMCs multiplied by PDFs in Mellin space.
        """
        CF(2,L,3) = F2/x,FL/x,F3  
        """ 
        if evolve: conf['pdf'].evolve(Q2)
  
        g =conf['pdf'].storage[Q2]['g']     # gluon[N] were N mellin contour
        if Nf==None: Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
        a=conf['alphaS'].get_a(Q2)   # a=alphaS/4pi

        if stf=='F2':
            CQ = self.C2Q[0] + a*self.order*self.C2Q[1]  # hard kernels 
            CG = self.C2G[0] + a*self.order*self.C2G[1]
            q=np.copy(conf['pdf'].storage[Q2]['qp'])    # (q + qb)[N], q[1] up, q[2] down, q[3] strange, q[4] charm, q[5] bottom, q[0]=q[6]=0
            aX='ap'

 
        elif stf=='FL':
            CQ = a*self.order*self.CLQ[1]
            CG = a*self.order*self.CLG[1]
            q=np.copy(conf['pdf'].storage[Q2]['qp'])    # (q + qb)[N]
            aX='ap'
  
        elif stf=='F3':
            CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
            CG = 0
            q=np.copy(conf['pdf'].storage[Q2]['qm'])    # (q - qb)[N]
            aX='am'
  
        if nucleon=="neutron":
            qup=np.copy(q[1])
            qdn=np.copy(q[2])
            q[1]=qdn
            q[2]=qup
  
        FX  = np.zeros(self.mellin.N.size,dtype=complex) # initialize FX[N] array
        for i in range(1,Nf+1):
            aXval = self.get_aX(aX,i,Q2)
            FX+=aXval*(CQ*q[i] + 2*CG*g)
  
        return FX
  
    def get_T4FX(self,stf,nucleon,x,Q2):
        if   stf=='F2': return self.t4F2.H[nucleon]/Q2
        elif stf=='FL': return self.t4FL.H[nucleon]/Q2
        elif stf=='F3': return self.t4F3.H[nucleon]/Q2
 
    def get_T4WX(self,stf,nucleon,x,Q2):
        if   stf=='W2': return self.t4W2.H[nucleon]/Q2
        elif stf=='WL': return self.t4WL.H[nucleon]/Q2
        elif stf=='W3': return self.t4W3.H[nucleon]/Q2

    def get_FXoff(self,stf,nucleon):
        if   stf=='F2': return self.F2off.off[nucleon]
        elif stf=='FL': return self.FLoff.off[nucleon]
        elif stf=='F3': return self.F3off.off[nucleon] 
  
    def get_T2CWX(self,stf,nucleon,Q2,sign,evolve=True):  
        """
        CF(2,L,3) = W2/x,WL/x,W3  
        """ 
        if evolve: conf['pdf'].evolve(Q2)
        g =conf['pdf'].storage[Q2]['g']     # gluon[N] were N mellin contour
        Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
        a=conf['alphaS'].get_a(Q2)   # a=alphaS/4pi
  
        if stf=='W2':
            CQ = self.C2Q[0] + a*self.order*self.C2Q[1]  # hard kernels 
            CG = self.C2G[0] + a*self.order*self.C2G[1]
  
        elif stf=='WL':
            CQ = a*self.order*self.CLQ[1]
            CG = a*self.order*self.CLG[1]
  
        elif stf=='W3':
            CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
            CG = 0
  
        u=conf['pdf'].storage[Q2]['u']
        d=conf['pdf'].storage[Q2]['d']
        s=conf['pdf'].storage[Q2]['s']
        c=conf['pdf'].storage[Q2]['c']
        b=conf['pdf'].storage[Q2]['b']
  
        ub=conf['pdf'].storage[Q2]['ub']
        db=conf['pdf'].storage[Q2]['db']
        sb=conf['pdf'].storage[Q2]['sb']
        cb=conf['pdf'].storage[Q2]['cb']
        bb=conf['pdf'].storage[Q2]['bb']
  
        U =(CQ*u  + CG*g)
        D =(CQ*d  + CG*g) + (CQ*s  + CG*g) 
        UB=(CQ*ub + CG*g)
        DB=(CQ*db + CG*g) + (CQ*sb + CG*g)
  
        if Nf>3:
            U += CQ*c  + CG*g
            UB+= CQ*cb + CG*g
  
        # cannot produce a top
        #if Nf>4:
        #  D+=  CQ*b  + CG*g
        #  DB+= CQ*bb + CG*g 
  
        if sign==+1:
            if   stf=='W2': return D+UB
            elif stf=='WL': return D+UB
            elif stf=='W3': return D-UB
        elif sign==-1:
            if   stf=='W2': return U+DB
            elif stf=='WL': return U+DB
            elif stf=='W3': return U-DB
  
    #--unpolarized functions
  
    def get_FXN(self,x,Q2,stf='F2',twist=2,nucleon='proton',Nf=None,tmc=False,precalc=False,evolve=True):
 
 
        if twist==2:
  
            if tmc==False:    #Since no TMCs, this simply inverts the Mellin moment of the SF without TMCs.
                if   stf=='F2': FX= x*self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)  #Extra factor of x.
                elif stf=='FL': FX= x*self.get_T2CFX('FL',nucleon,Q2,Nf,evolve)
                elif stf=='F3': FX=   self.get_T2CFX('F3',nucleon,Q2,Nf,evolve)
                FX*=x**(-self.mellin.N) 
  
            elif tmc=='GP':
                if   stf=='F2': FX= self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'22')
                elif stf=='FL': FX= self.get_T2CFX('FL',nucleon,Q2,Nf,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'LL')\
                                   +self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'L2')
                elif stf=='F3': FX= self.get_T2CFX('F3',nucleon,Q2,Nf,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'33')
  
            elif tmc=='AOT':  
              if   stf=='F2': FX= self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'22')
              elif stf=='FL': FX= self.get_T2CFX('FL',nucleon,Q2,Nf,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'LL')\
                                 +self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'L2')
              elif stf=='F3': FX= self.get_T2CFX('F3',nucleon,Q2,Nf,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'33')
  
        if twist==4: 
   
          if self.ht==0: FX = x**(-self.mellin.N)*self.get_T4FX(stf,nucleon,x,Q2)
          elif self.ht==1:
            if nucleon=='proton':  FX = x**(-self.mellin.N)*self.get_T4FX(stf,nucleon,x,Q2)
            if nucleon=='neutron': return 0
          elif self.ht==2:
            if nucleon=='proton':  FX = x**(-self.mellin.N)*self.get_FXN(x,Q2,stf=stf,twist=2,nucleon=nucleon,tmc=tmc)*self.get_T4FX(stf,nucleon,x,Q2)
            if nucleon=='neutron': return 0
 
        return self.mellin.invert(1,FX) # here x=1 since the factor x**(-N) is included in the moments
  
    def get_WXN(self,x,Q2,stf='W2',sign=+1,twist=2,nucleon='proton',tmc=False,precalc=False,evolve=True):
  
        if twist==2:

            if tmc==False:
                if   stf=='W2': WX= x*self.get_T2CWX('W2',nucleon,Q2,sign,evolve)
                elif stf=='WL': WX= x*self.get_T2CWX('WL',nucleon,Q2,sign,evolve)
                elif stf=='W3': WX=   self.get_T2CWX('W3',nucleon,Q2,sign,evolve)
                WX*=x**(-self.mellin.N)
            elif tmc=='GP':
                if   stf=='W2': WX= self.get_T2CWX('W2',nucleon,Q2,sign,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'22')
                elif stf=='WL': WX= self.get_T2CWX('WL',nucleon,Q2,sign,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'LL')\
                                   +self.get_T2CWX('W2',nucleon,Q2,sign,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'L2')
                elif stf=='W3': WX= self.get_T2CWX('W3',nucleon,Q2,sign,evolve)*self.get_TMC_GP(self.mellin.N,x,Q2,'33')
            elif tmc=='AOT':
                if   stf=='W2': WX= self.get_T2CWX('W2',nucleon,Q2,sign,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'22')
                elif stf=='WL': WX= self.get_T2CWX('WL',nucleon,Q2,sign,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'LL')\
                                   +self.get_T2CWX('W2',nucleon,Q2,sign,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'L2')
                elif stf=='W3': WX= self.get_T2CWX('W3',nucleon,Q2,sign,evolve)*self.get_TMC_AOT(self.mellin.N,x,Q2,'33')

 
        if twist==4:  

          if self.ht==0: WX = x**(-self.mellin.N)*self.get_T4WX(stf,nucleon,x,Q2)
          elif self.ht==1:
            if nucleon=='proton':  WX = x**(-self.mellin.N)*self.get_T4WX(stf,nucleon,x,Q2)
            if nucleon=='neutron': return 0
          elif self.ht==2:
            if nucleon=='proton':  WX = x**(-self.mellin.N)*self.get_WXN(x,Q2,stf=stf,twist=2,nucleon=nucleon,tmc=tmc)*self.get_T4WX(stf,nucleon,x,Q2)
            if nucleon=='neutron': return 0

        return self.mellin.invert(1,WX) # here x=1 since the factor x**(-N) is included in the moments
  
    def get_F2D(self,x,Q2,twist=2,Nf=None,tmc=False,nuc=False,offshell=False,evolve=True): 
    
        if offshell==False:

            if twist==2:
                T2F2 =  self.get_T2CFX('F2','proton',Q2,Nf,evolve)+self.get_T2CFX('F2','neutron',Q2,Nf,evolve)
                if    tmc==False and nuc==False: f22 = 1
                elif  tmc==False and nuc==True : f22 = self.tables[self.get_tname('2',fXX='f22',proj=None)].retrieve(x,Q2)
                elif  tmc!=False and nuc==True : f22 = self.tables[self.get_tname('2',fXX='f22',proj='22')].retrieve(x,Q2)
            
                F2D= f22*T2F2 
  
            if twist==4:
                #--no such thing as tmc for twist 4
                if self.ht==0:

                    T4F2 = self.get_T4FX('F2','proton',Q2) + self.get_T4FX('F2','neutron',Q2)

                    if nuc==False : f22==1
                    else          : f22=self.tables[self.get_tname('4',fXX='f22',proj=None)].retrieve(x,Q2)

                    F2D= f22*T4F2
                
                elif self.ht==1: F2D =  x**(-self.mellin.N)*self.get_T4FX('F2','deuteron',x,Q2)
                elif self.ht==2: F2D =  x**(-self.mellin.N)*self.get_F2D(x,Q2,twist=2,tmc=tmc,nuc=nuc,offshell=False)*self.get_T4FX('F2','deuteron',x,Q2)
 
            return self.mellin.invert(1,F2D)

        if offshell==True:

            if nuc==False:

                print ('ERR: There cannot be offshell corrections without nuc. smearing')
                return

            if twist==4: return 0

            F2Doff  = complex(0,0)
            F2Doffc = complex(0,0)
 
            for nucleon in ['proton','neutron']:

                TXF2 = self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)

                #--offshell corrections
                df2   = self.get_FXoff('F2',nucleon) 
                df2c  = self.get_FXoff('F2',nucleon+'_conj') 

                #--deuteron offshell smearing grids
                if tmc==False: 
                   f22off  =  self.tables[self.get_tname('2',fXX='f22',proj=None,offshell=True,conj=False)].retrieve(x,Q2)
                   f22offc  = self.tables[self.get_tname('2',fXX='f22',proj=None,offshell=True,conj=True)].retrieve(x,Q2)
                if tmc!=False:             
                   f22off   = self.tables[self.get_tname('2',fXX='f22',proj='22',offshell=True,conj=False)].retrieve(x,Q2)
                   f22offc  = self.tables[self.get_tname('2',fXX='f22',proj='22',offshell=True,conj=True)].retrieve(x,Q2)

                F2Doff  += x**(-self.dmellin.M)*df2*f22off*TXF2
                F2Doffc += x**(-np.conj(self.dmellin.M))*df2c*f22offc*TXF2  #--extra powers of x**-M are needed? (may be incorrect)
            
            return self.dmellin.invert(1,F2Doff,F2Doffc)  #--extra powers of x**-N not needed? (may be incorrect)

    def get_FLD(self,x,Q2,twist=2,Nf=None,tmc=False,nuc=False,evolve=True,offshell=False):
 
        if offshell==False:

              if twist==2: 

                  T2F2 = self.get_T2CFX('F2','proton',Q2,Nf,evolve)+self.get_T2CFX('F2','neutron',Q2,Nf,evolve)
                  T2FL = self.get_T2CFX('FL','proton',Q2,Nf,evolve)+self.get_T2CFX('FL','neutron',Q2,Nf,evolve)

                  if    tmc==False and nuc==False:
                      fL2p22 = 0
                      fLLpLL = 1
                      fLLpL2 = 0
                  elif  tmc==False and nuc==True :
                      fL2p22  = self.tables[self.get_tname('2',fXX='fL2',proj=None)].retrieve(x,Q2)
                      fLLpLL  = self.tables[self.get_tname('2',fXX='fLL',proj=None)].retrieve(x,Q2)
                      fLLpL2  = 0
                  elif  tmc!=False  and nuc==True : 
                      fL2p22=self.tables[self.get_tname('2',fXX='fL2',proj='22')].retrieve(x,Q2)
                      fLLpLL=self.tables[self.get_tname('2',fXX='fLL',proj='LL')].retrieve(x,Q2)
                      fLLpL2=self.tables[self.get_tname('2',fXX='fLL',proj='L2')].retrieve(x,Q2)

                  FLD= fL2p22*T2F2 + fLLpLL*T2FL + fLLpL2*T2F2
  
              if twist==4:

                  #--no such thing as tmc for twist 4

                  if self.ht==0:
                      T4F2 = self.get_T4FX('F2','proton',x,Q2) + self.get_T4FX('F2','neutron',x,Q2)
                      T4FL = self.get_T4FX('FL','proton',x,Q2) + self.get_T4FX('FL','neutron',x,Q2)

                      if nuc==False:  
                          fL2p22  = 0
                          fLLpLL  = 1
                          fLLpL2  = 0
                      else:  
                          fL2p22  = self.tables[self.get_tname('4',fXX='fL2',proj=None)].retrieve(x,Q2)
                          fLLpLL  = self.tables[self.get_tname('4',fXX='fLL',proj=None)].retrieve(x,Q2)
                          fLLpL2  = 0

                      FLD= fL2p22*T2F2 + fLLpLL*T2FL + fLLpL2*T2F2

                  elif self.ht==1: FLD =  x**(-self.mellin.N)*self.get_T4FX('FL','deuteron',x,Q2)
                  elif self.ht==2: FLD =  x**(-self.mellin.N)*self.get_FLD(x,Q2,twist=2,tmc=tmc,nuc=nuc,offshell=False)*self.get_T4FX('FL','deuteron',x,Q2)

              return self.mellin.invert(1,FLD)

        if offshell==True:

            if nuc==False:

                print ('ERR: There cannot be offshell corrections without nuc. smearing')
                return

            if twist==4: return 0

            FLDoff  = complex(0,0)
            FLDoffc = complex(0,0)
 
            for nucleon in ['proton','neutron']:

                TXF2 = self.get_T2CFX('F2',nucleon,Q2,Nf,evolve)
                TXFL = self.get_T2CFX('FL',nucleon,Q2,Nf,evolve)

                #--offshell corrections
                df2  = self.get_FXoff('F2',nucleon)
                dfL  = self.get_FXoff('FL',nucleon)
                df2c  = self.get_FXoff('F2',nucleon+'_conj')
                dfLc  = self.get_FXoff('FL',nucleon+'_conj')

                #--deuteron offshell smearing tables
                if tmc==False:
                    fL2p22off   = self.tables[self.get_tname('2',fXX='fL2',proj=None,offshell=True,conj=False)].retrieve(x,Q2)
                    fLLpLLoff   = self.tables[self.get_tname('2',fXX='fLL',proj=None,offshell=True,conj=False)].retrieve(x,Q2)
                    fLLpL2off   = 0
                    fL2p22offc  = self.tables[self.get_tname('2',fXX='fL2',proj=None,offshell=True,conj=True)].retrieve(x,Q2)
                    fLLpLLoffc  = self.tables[self.get_tname('2',fXX='fLL',proj=None,offshell=True,conj=True)].retrieve(x,Q2)
                    fLLpL2offc  = 0
                if tmc!=False:
                    fL2p22off   = self.tables[self.get_tname('2',fXX='fL2',proj='22',offshell=True,conj=False)].retrieve(x,Q2)
                    fLLpLLoff   = self.tables[self.get_tname('2',fXX='fLL',proj='LL',offshell=True,conj=False)].retrieve(x,Q2)
                    fLLpL2off   = self.tables[self.get_tname('2',fXX='fLL',proj='L2',offshell=True,conj=False)].retrieve(x,Q2)
                    fL2p22offc  = self.tables[self.get_tname('2',fXX='fL2',proj='22',offshell=True,conj=True)].retrieve(x,Q2)
                    fLLpLLoffc  = self.tables[self.get_tname('2',fXX='fLL',proj='LL',offshell=True,conj=True)].retrieve(x,Q2)
                    fLLpL2offc  = self.tables[self.get_tname('2',fXX='fLL',proj='L2',offshell=True,conj=True)].retrieve(x,Q2)

                FLDoff   += x**(-self.dmellin.M)*df2*fL2p22off*TXF2 + dfL*fLLpLLoff*TXFL + df2*fLLpL2*TXF2
                FLDoffc  += x**(-np.conj(self.dmellin.M))*df2c*fL2p22offc*TXF2 + dfLc*fLLpLLoffc*TXFL + df2c*fLLpL2c*TXF2

            return self.dmellin.invert(1,FLDoff,FLDoffc)

    def get_F3D(self,x,Q2,twist=2,Nf=None,tmc=False,nuc=False,offshell=False,evolve=True): 
    
        if offshell==False:

            if twist==2:
                T2F3 =  self.get_T2CFX('F3','proton',Q2,Nf,evolve)+self.get_T2CFX('F3','neutron',Q2,Nf,evolve)
                if    tmc==False and nuc==False: f33 = 1
                elif  tmc==False and nuc==True : f33 = self.tables[self.get_tname('2',fXX='f33',proj=None)].retrieve(x,Q2)
                elif  tmc!=False and nuc==True : f33 = self.tables[self.get_tname('2',fXX='f33',proj='33')].retrieve(x,Q2)
            
                F3D= f33*T2F3 
  
            if twist==4:

                #--no such thing as tmc for twist 4

                if self.ht==0:
                    T4F3 = self.get_T4FX('F3','proton',Q2) + self.get_T4FX('F3','neutron',Q2)

                    if nuc==False : f33==1
                    else          : f33=self.tables[self.get_tname('4',fXX='f33',proj=None)].retrieve(x,Q2)

                    F3D= f33*T4F3

                elif self.ht==1: F3D =  x**(-self.mellin.N)*self.get_T4FX('F3','deuteron',x,Q2)
                elif self.ht==2: F3D =  x**(-self.mellin.N)*self.get_F3D(x,Q2,twist=2,tmc=tmc,nuc=nuc,offshell=False)*self.get_T4FX('F3','deuteron',x,Q2)
  
            return self.mellin.invert(1,F3D)


        if offshell==True:

            if nuc==False:

                print ('ERR: There cannot be offshell corrections without nuc. smearing')
                return

            if twist==4: return 0

            F3Doff  = complex(0,0)
            F3Doffc = complex(0,0)
 
            for nucleon in ['proton','neutron']:

                TXF3 = self.get_T2CFX('F3',nucleon,Q2,Nf,evolve)

                #--offshell corrections
                df3   = self.get_FXoff('F3',nucleon) 
                df3c  = self.get_FXoff('F3',nucleon+'_conj') 

                #--deuteron offshell smearing grids
                if tmc==False:            
                    f33off   = self.tables[self.get_tname('2',fXX='f33',proj=None,offshell=True,conj=False)].retrieve(x,Q2)
                    f33offc  = self.tables[self.get_tname('2',fXX='f33',proj=None,offshell=True,conj=True)].retrieve(x,Q2)
                if tmc!=False:           
                    f33off   = self.tables[self.get_tname('2',fXX='f33',proj='33',offshell=True,conj=False)].retrieve(x,Q2)
                    f33offc  = self.tables[self.get_tname('2',fXX='f33',proj='33',offshell=True,conj=True)].retrieve(x,Q2)

                F3Doff  += x**(-self.dmellin.M)*df3*f33off*TXF3
                F3Doffc += x**(-np.conj(self.dmellin.M))*df3c*f33offc*TXF3

            return self.dmellin.invert(1,F3Doff,F3Doffc)

    def _dev_get_F2H(self,x,Q2,twist=2,Nf=None,tmc=True,evolve=True):    #F2 for Hydrogen
  
        if twist==2: 
            T2F2p = self.get_T2CFX('F2','proton',Q2,Nf,evolve)
            T2F2n = self.get_T2CFX('F2','neutron',Q2,Nf,evolve)
            if tmc=='GP':
                f22p=self.tables[self.get_tname('2',fXX='f22p',proj='22')].retrieve(x,Q2)
                f22n=self.tables[self.get_tname('2',fXX='f22n',proj='22')].retrieve(x,Q2)
            elif tmc==False:
                f22p=self.tables[self.get_tname('2',fXX='f22p',proj=None)].retrieve(x,Q2)
                f22n=self.tables[self.get_tname('2',fXX='f22n',proj=None)].retrieve(x,Q2)
            #Need tmc=='AOT'
            F2H = 2*f22p*T2F2p + f22n*T2F2n
          
  
        if twist==4:  
            T4F2p = self.get_T4FX('F2','proton',Q2)
            T4F2n = self.get_T4FX('F2','neutron',Q2)
            f22p=self.tables[self.get_tname('4',fXX='f22p',proj=None)].retrieve(x,Q2)
            f22n=self.tables[self.get_tname('4',fXX='f22n',proj=None)].retrieve(x,Q2)
            F2H = 2*f22p*T4F2p + f22n*T4F2n
  
        return self.mellin.invert(1,F2H)
  
    def _dev_get_FLH(self,x,Q2,twist=2,Nf=None,tmc=True,evolve=True):
  
        if twist==2: 

            T2F2p = self.get_T2CFX('F2','proton',Q2,Nf,evolve)
            T2F2n = self.get_T2CFX('F2','neutron',Q2,Nf,evolve)
            T2FLp = self.get_T2CFX('FL','proton',Q2,Nf,evolve)
            T2FLn = self.get_T2CFX('FL','neutron',Q2,Nf,evolve)
            if tmc=='GP': 
                fL2p22p=self.tables[self.get_tname('2',fXX='fL2p',proj='22')].retrieve(x,Q2)
                fLLpLLp=self.tables[self.get_tname('2',fXX='fLLp',proj='LL')].retrieve(x,Q2)
                fLLpL2p=self.tables[self.get_tname('2',fXX='fLLp',proj='L2')].retrieve(x,Q2)
                fL2p22n=self.tables[self.get_tname('2',fXX='fL2n',proj='22')].retrieve(x,Q2)
                fLLpLLn=self.tables[self.get_tname('2',fXX='fLLn',proj='LL')].retrieve(x,Q2)
                fLLpL2n=self.tables[self.get_tname('2',fXX='fLLn',proj='L2')].retrieve(x,Q2)
                FLH = 2*(fL2p22p*T2F2p+fLLpLLp*T2FLp+fLLpL2p*T2F2p)\
                        + (fL2p22n*T2F2n+fLLpLLn*T2FLn+fLLpL2n*T2F2n)
            elif tmc==False: 
                fL2p=self.tables[self.get_tname('2',fXX='fL2p',proj=None)].retrieve(x,Q2)
                fLLp=self.tables[self.get_tname('2',fXX='fLLp',proj=None)].retrieve(x,Q2)
                fL2n=self.tables[self.get_tname('2',fXX='fL2n',proj=None)].retrieve(x,Q2)
                fLLn=self.tables[self.get_tname('2',fXX='fLLn',proj=None)].retrieve(x,Q2)
                FLH = 2*(fL2p*T2F2p+fLLp*T2FLp) + (fL2n*T2F2n+fLLn*T2FLn)
  
        if twist==4:  

            T4F2p = self.get_T4FX('F2','proton',Q2)
            T4F2n = self.get_T4FX('F2','neutron',Q2)
            T4FLp = self.get_T4FX('FL','proton',Q2)
            T4FLn = self.get_T4FX('FL','neutron',Q2)
            fL2p=self.tables[self.get_tname('4',fXX='fL2p',proj=None)].retrieve(x,Q2)
            fL2n=self.tables[self.get_tname('4',fXX='fL2n',proj=None)].retrieve(x,Q2)
            fLLp=self.tables[self.get_tname('4',fXX='fLLp',proj=None)].retrieve(x,Q2)
            fLLn=self.tables[self.get_tname('4',fXX='fLLn',proj=None)].retrieve(x,Q2)
            FLH = 2*(fL2p*T4F2p+fLLp*T4FLp) + (fL2n*T4F2n+fLLn*T4FLn)
  
        return self.mellin.invert(1,FLH)
  
    def _dev_get_F3H(self,x,Q2,twist=2,tmc=True,evolve=True):
  
        if twist==2: 
            T2F3p = self.get_T2CFX('F3','proton',Q2,evolve)
            T2F3n = self.get_T2CFX('F3','neutron',Q2,evolve)
            if tmc=='GP':
                f3p33p=self.tables[self.get_tname('2',fXX='f3p',proj='33')].retrieve(x,Q2)
                f3p33n=self.tables[self.get_tname('2',fXX='f3n',proj='33')].retrieve(x,Q2)
                F3H = 2*f3p33p*T2F3p + f3p33n*T2F3n
            elif tmc==False:
                f3p=self.tables[self.get_tname('2',fXX='f3p',proj=None)].retrieve(x,Q2)
                f3n=self.tables[self.get_tname('2',fXX='f3n',proj=None)].retrieve(x,Q2)
                T2F3H+= 2*f3p*T2F3p + f3n*T2F3n
  
        if twist==4:  
            T4F3p = self.get_T4FX('F3','proton',Q2)
            T4F3n = self.get_T4FX('F3','neutron',Q2)
            f3p=self.tables[self.get_tname('4',fXX='f3p',proj=None)].retrieve(x,Q2)
            f3n=self.tables[self.get_tname('4',fXX='f3n',proj=None)].retrieve(x,Q2)
            F3H = 2*f3p*T4F3p + f3n*T4F3n
  
        return self.mellin.invert(1,F3H)
  
    #--routines to setup/fill/save mellin tables
  
    def get_tname(self,twist,fXX,proj,offshell=False,conj=False):
      if offshell==False:
        if proj!=None: return 'T%s-%s-%s'%(twist,fXX,proj)
        else:   return 'T%s-%s'%(twist,fXX)
      if offshell==True:
        if conj==False:
          if proj!=None: return 'T%s-%s-%s-offshell'%(twist,fXX,proj)
          else:   return 'T%s-%s-offshell'%(twist,fXX)
        if conj==True:
          if proj!=None: return 'T%s-%s-%s-offshell-conj'%(twist,fXX,proj)
          else:   return 'T%s-%s-offshell-conj'%(twist,fXX)
 
    def load_table(self,target,twist,fXX,proj,offshell=False,conj=False):
        name=self.get_tname(twist,fXX,proj,offshell,conj)
        if offshell==False:
          print('loading idis table: %s/%s/%s'%(self.path2distab,target,name))
          self.tables[name]= STORAGE('%s/%s'%(self.path2distab,target),name)

        if offshell==True:
          if conj==False:
            print('loading idis offshell table: %s/%s/%s'%(self.path2distab,target,name))
            self.tables[name]= STORAGE('%s/%s'%(self.path2distab,target),name)
          if conj==True:
            print('loading idis offshell table: %s/%s/%s'%(self.path2distab,target,name))
            self.tables[name]= STORAGE('%s/%s'%(self.path2distab,target),name)
  
    def setup_tables(self):

        tmc = False
        nuc = False
        offshell = False
        if 'tmc'  in conf: tmc = conf['tmc']
        if 'nuc'  in conf: nuc = conf['nuc']
        if 'offshell' in conf: offshell = conf['offshell']
        self.tables={}

        # for F2D
        if tmc==False and nuc==True: 

            self.load_table('deuteron',2,'f22',None) 

            if offshell==True:
                self.load_table('deuteron',2,'f22',None,offshell=True,conj=False) 
                self.load_table('deuteron',2,'f22',None,offshell=True,conj=True) 

        if tmc!=False and nuc==True: 

            self.load_table('deuteron',2,'f22','22') 

            if offshell==True:
                #self.load_table('deuteron',2,'f22','22',offshell=True,conj=False) 
                self.load_table('deuteron',2,'f22','22',offshell=True,conj=True)

 
        # for FLD
        
        if tmc==False and nuc==True: 
            return
            self.load_table('deuteron',2,'fL2',None)
            self.load_table('deuteron',2,'fLL',None)

            #if offshell==True:

               #self.load_table('deuteron',2,'fL2',None,offshell=True,conj=False)
               #self.load_table('deuteron',2,'fLL',None,offshell=True,conj=False)
               #self.load_table('deuteron',2,'fL2',None,offshell=True,conj=True)
               #self.load_table('deuteron',2,'fLL',None,offshell=True,conj=True)

        if tmc!=False and nuc==True: 
            return
            self.load_table('deuteron',2,'fL2','22')
            self.load_table('deuteron',2,'fLL','LL')
            self.load_table('deuteron',2,'fLL','L2')

  
            #if offshell==True:  

                #self.load_table('deuteron',2,'fL2','22',offshell=True,conj=False)
                #self.load_table('deuteron',2,'fLL','LL',offshell=True,conj=False)
                #self.load_table('deuteron',2,'fLL','L2',offshell=True,conj=False)
                #self.load_table('deuteron',2,'fL2','22',offshell=True,conj=True)
                #self.load_table('deuteron',2,'fLL','LL',offshell=True,conj=True)
                #self.load_table('deuteron',2,'fLL','L2',offshell=True,conj=True)

        """
        # for F3D

        if tmc==False and nuc==True:

            self.load_table('deuteron',2,'f33',None)

            if offshell==True:

                self.load_table('deuteron',2,'f33',None,offshell=True,conj=False)
                self.load_table('deuteron',2,'f33',None,offshell=True,conj=True)


        if tmc!=False and nuc==True and offshell==False:

            self.load_table('deuteron',2,'f33','33')

            if offshell==True:

                self.load_table('deuteron',2,'f33','33',offshell=True,conj=False)
                self.load_table('deuteron',2,'f33','33',offshell=True,conj=True)
        """
        ## for F2H
        #self.load_table('helium',2,'f22p','22')
        #self.load_table('helium',2,'f22n','22')
        #self.load_table('helium',2,'f22p',None)
        #self.load_table('helium',2,'f22n',None)
        #self.load_table('helium',4,'f22p',None)
        #self.load_table('helium',4,'f22n',None)
  
        ## for FLH
        #self.load_table('helium',2,'fL2p','22')
        #self.load_table('helium',2,'fLLp','LL')
        #self.load_table('helium',2,'fLLp','L2')
        #self.load_table('helium',2,'fL2n','22')
        #self.load_table('helium',2,'fLLn','LL')
        #self.load_table('helium',2,'fLLn','L2')
        #self.load_table('helium',2,'fL2p',None)
        #self.load_table('helium',2,'fLLp',None)
        #self.load_table('helium',2,'fLLp',None)
        #self.load_table('helium',2,'fL2n',None)
        #self.load_table('helium',2,'fLLn',None)
        #self.load_table('helium',2,'fLLn',None)
        #self.load_table('helium',4,'fL2p',None)
        #self.load_table('helium',4,'fLLp',None)
        #self.load_table('helium',4,'fL2n',None)
        #self.load_table('helium',4,'fLLn',None)
  
        ## for F3H
        #self.load_table('helium',2,'f3p','33')
        #self.load_table('helium',2,'f3n','33')
        #self.load_table('helium',2,'f3p',None)
        #self.load_table('helium',2,'f3n',None)
        #self.load_table('helium',4,'f3p',None)
        #self.load_table('helium',4,'f3n',None)
 
    def fill_table(self,x,Q2,target,twist,tmc,offshell,conj,fXX,proj):
        name=self.get_tname(twist,fXX,proj,offshell,conj)
        func=lambda x,Q2:self.get_smf(x,Q2,target,twist,tmc,offshell,conj,fXX,proj)
        self.tables[name].register(func,x,Q2,name)
        self.modified_tables[name]=True
  
    def fill_F2D(self,x,Q2,tmc=False,offshell=True,conj=False):

        #--t2:
        if tmc==False:
            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'f22',None)

        else:
            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'f22','22')

        #--t4:  remove any tmc effect

        #if self.ht==0: self.fill_table(x,Q2,'deuteron',4,False,offshell,conj,'f22',None)
  
    def fill_FLD(self,x,Q2,tmc=False,offshell=False,conj=False):
        return
        #--t2
        if tmc==False:

            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'fL2',None)
            #self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'fLL',None)

        else:

            #self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'fL2','22')
            #self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'fLL','LL')
            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'fLL','L2')

        #--t4:  remove any tmc effect

        #if self.ht==0: self.fill_table(x,Q2,'deuteron',4,False,offshell,'fL2',None)
        #if self.ht==0: self.fill_table(x,Q2,'deuteron',4,False,offshell,'fLL',None)

    def fill_F3D(self,x,Q2,tmc=False,offshell=False):

        #--t2:
        if tmc==False:

            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'f33',None)

        else:

            self.fill_table(x,Q2,'deuteron',2,tmc,offshell,conj,'f33','33')

        #--t4:  remove any tmc effect

        #if self.ht==0: self.fill_table(x,Q2,'deuteron',4,False,offshell,'f33',None)
  
    def _dev_fill_F2H(self,x,Q2):
        #self.fill_table(x,Q2,'helium',2,'f22p','22')
        #self.fill_table(x,Q2,'helium',2,'f22n','22')
        self.fill_table(x,Q2,'helium',2,'f22p',None)
        self.fill_table(x,Q2,'helium',2,'f22n',None)
        #self.fill_table(x,Q2,'helium',4,'f22p',None)
        #self.fill_table(x,Q2,'helium',4,'f22n',None)
  
    def _dev_fill_FLH(self,x,Q2):
        #self.fill_table(x,Q2,'helium',2,'fL2p','22')
        #self.fill_table(x,Q2,'helium',2,'fLLp','LL')
        #self.fill_table(x,Q2,'helium',2,'fLLp','L2')
        #self.fill_table(x,Q2,'helium',2,'fL2n','22')
        #self.fill_table(x,Q2,'helium',2,'fLLn','LL')
        #self.fill_table(x,Q2,'helium',2,'fLLn','L2')
        self.fill_table(x,Q2,'helium',2,'fL2p',None)
        self.fill_table(x,Q2,'helium',2,'fLLp',None)
        self.fill_table(x,Q2,'helium',2,'fLLp',None)
        self.fill_table(x,Q2,'helium',2,'fL2n',None)
        self.fill_table(x,Q2,'helium',2,'fLLn',None)
        self.fill_table(x,Q2,'helium',2,'fLLn',None)
        #self.fill_table(x,Q2,'helium',4,'fL2p',None)
        #self.fill_table(x,Q2,'helium',4,'fLLp',None)
        #self.fill_table(x,Q2,'helium',4,'fL2n',None)
        #self.fill_table(x,Q2,'helium',4,'fLLn',None)
  
    def _dev_fill_F3H(self,x,Q2):
        #self.fill_table(x,Q2,'helium',2,'f3p','33')
        #self.fill_table(x,Q2,'helium',2,'f3n','33')
        self.fill_table(x,Q2,'helium',2,'f3p',None)
        self.fill_table(x,Q2,'helium',2,'f3n',None)
        #self.fill_table(x,Q2,'helium',4,'f3p',None)
        #self.fill_table(x,Q2,'helium',4,'f3n',None)
  
    def gen_tables(self,tabname,tmc,offshell=False,conj=False,ncores=2):
        self.modified_tables={}
        tab=conf[tabname]
        K=sorted(tab.keys())
        for k in K:
            npts=len(tab[k]['value'])
            for i in range(npts):
                #print('generating tables for %d i=%d of %d'%(k,i,npts))
                x=tab[k]['X'][i]
                Q2=tab[k]['Q2'][i]
                target=tab[k]['target'][i].strip()
                obs=tab[k]['obs'][i]
  
                proton=False
                neutron=False
                deuteron=False
                helium=False
  
                if target=='p' or target=='proton'   : proton=True 
                if target=='n' or target=='neutron'  : neutron=True 
                if target=='d' or target=='deuteron' : deuteron=True 
                if target=='h' or target=='helium'   : helium=True 
                if target=='d/p': deuteron=True 
                if target=='n/d': deuteron=True 
  
                if deuteron: 
                    if obs=='F2':self.fill_F2D(x,Q2,tmc,offshell,conj=conj)
                    if obs=='sig_r': 
                        self.fill_F2D(x,Q2,tmc=tmc,offshell=offshell,conj=conj)
                        self.fill_FLD(x,Q2,tmc=tmc,offshell=offshell,conj=conj)
                        self.fill_F3D(x,Q2,tmc=tmc,offshell=offshell,conj=conj)
                    if obs=='F2d/F2p':self.fill_F2D(x,Q2,tmc,offshell,conj=conj)
                    if obs=='F2n/F2d':self.fill_F2D(x,Q2,tmc,offshell,conj=conj)
                    if obs=='multiplicity':
                        self.fill_F2D(x,Q2,tmc,offshell,conj=conj)
                        self.fill_FLD(x,Q2,tmc,offshell,conj=conj)
                    if any([obs==_obs for _obs in ['A1','A2','Apa','Ape','Atpe']]):
                        self.fill_F2D(x,Q2,tmc,offshell,conj=conj)
                        self.fill_FLD(x,Q2,tmc,offshell,conj=conj)
  
                #if helium: 
                #    if obs=='F2': self.fill_F2H(x,Q2)
                #    if obs=='FL': self.fill_FLH(x,Q2)
                #    if obs=='F3': self.fill_F3H(x,Q2)
                #    if obs=='sig_r': 
                #        self.fill_F2H(x,Q2)
                #        self.fill_FLH(x,Q2)
                #        self.fill_F3H(x,Q2)
  
                #    if any([obs==_obs for _obs in ['A1','A2','Apa','Ape','Atpe']]):
                #        self.fill_F2H(x,Q2)
                #        self.fill_FLH(x,Q2)

      
        for _ in self.tables.keys():
            if _ in self.modified_tables:
                npts=len(self.tables[_].requests)
                print('processing table %s entries=%d'%(_,npts))

                self.tables[_].process_request(ncores)
                self.tables[_].Save()
  
  

