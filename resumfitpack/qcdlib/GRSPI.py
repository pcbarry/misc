#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.integrate import fixed_quad,quad
from scipy.special import gamma as Gamma

class GRSPI(object):

  def __init__(self,evolution=True):
    self.evolution=evolution
    self.s=lambda Q2: np.log(np.log(Q2/0.299**2)/np.log(0.40/0.299**2))

  def get_valence(self,x,Q2):
    s=self.s(Q2)
    N=1.500+0.525*s-0.050*s**2
    a=0.560-0.034*s
    A=-0.357-0.458*s
    B=0.427+0.220*s
    D=0.475+0.550*s
    xv=N*x**a*(1+A*x**0.5+B*x)*(1-x)**D
    return 0.5*xv/x

  def get_glu(self,x,Q2):
    s=self.s(Q2)
    alpha=0.793
    a=1.418-0.215*s**0.5
    A=5.392+0.553*s-0.385*s**2
    C=11.548-4.316*s+0.382*s**2
    E=0.104+1.980*s
    beta=1.722
    b=0.0
    B=-11.928+1.844*s
    D=1.347+1.135*s
    Ep=2.375-0.188*s
    xg=(x**a*(A+B*x**0.5+C*x)*(np.log(1.0/x))**b+s**alpha*np.exp(-E+(Ep*s**beta*np.log(1.0/x))**0.5))*(1-x)**D
    return xg/x

  def get_light_sea(self,x,Q2):
    s=self.s(Q2)
    alpha=1.118
    a=0.111-0.326*s**0.5
    A=1.035-0.295*s
    C=4.111-1.575*s
    E=5.035+0.997*s
    beta=0.457
    b=-0.978-0.488*s**0.5
    B=-3.008+1.165*s
    D=6.192+0.705*s
    Ep=1.486+1.288*s
    xls=4*(x**a*(A+B*x**0.5+C*x)*(np.log(1/x))**b+s**alpha*np.exp(-E+(Ep*s**beta*np.log(1/x))**0.5))*(1-x)**D
    return xls/x

  def get_strange_sea(self,x,Q2):
    s=self.s(Q2)
    alpha=0.908
    a=-0.567-0.466*s
    B=4.403
    E=3.796+1.618*s
    beta=0.812
    A=-2.348+1.433*s
    D=2.061
    Ep=0.309+0.355*s
    xss=2*s**alpha/(np.log(1/x))**a*(1+A*x**0.5+B*x)*(1-x)**D*np.exp(-E+(Ep*s**beta*np.log(1/x))**0.5)
    return xss/x

  def get_pdfs(self,x,Q2):
    # for pi minus
    if self.evolution==True:
      val=self.get_valence(x,Q2)
      light_sea=self.get_light_sea(x,Q2)/4.0
      strange_sea=self.get_strange_sea(x,Q2)/2.0
      glu=self.get_glu(x,Q2)
    elif self.evolution==False:
      val=0.5*1.391*x**(-0.447)*(1-x)**0.426
      light_sea=0.417*x**(-0.793)*(1-2.466*x**0.5+3.855*x)*(1-x)**4.454
      strange_sea=0.0
      glu=5.90*x**0.270*(1-2.074*x**0.5+1.824*x)*(1-x)**1.290
    pdf=np.zeros(11)
    pdf[0]=glu
    pdf[1]=light_sea
    pdf[2]=val+light_sea
    pdf[3]=val+light_sea
    pdf[4]=light_sea
    pdf[5]=strange_sea
    pdf[6]=strange_sea
    return pdf

  def get_pdfsN(self,N,Q2):
    # for pi minus
    if self.evolution==True:
      print 'no analytic form for Q2!=Q20'
      sys.exit()
    else:
      val=1.391*Gamma(N-0.447)*Gamma(1.426)/Gamma(N+0.979)
      light_sea=0.417*Gamma(5.454)*(Gamma(N-0.793)/Gamma(N+4.661)\
        -2.466*Gamma(N-0.293)/Gamma(N+5.161)+3.855*Gamma(N+0.207)/Gamma(N+5.661))
      strange_sea=0.0
      glu=5.90*Gamma(2.290)*(Gamma(N+0.270)/Gamma(N+2.56)\
        -2.074*Gamma(N+0.770)/Gamma(N+3.06)+1.824*Gamma(N+1.270)/Gamma(N+3.56))
      pdf=np.zeros(11,dtype=complex)
      pdf[0]=glu
      pdf[1]=light_sea
      pdf[2]=val+light_sea
      pdf[3]=val+light_sea
      pdf[4]=light_sea
      pdf[5]=strange_sea
      pdf[6]=strange_sea
      return pdf


if __name__=="__main__":

  grs=GRSPI(evolution=False)
  print grs.get_pdfs(0.5,10.0)

  pdfs=lambda x: grs.get_pdfs(x,10.0)

  g=lambda x: pdfs(x)[0]
  u=lambda x: pdfs(x)[1]
  ub=lambda x: pdfs(x)[2]
  d=lambda x: pdfs(x)[3]
  db=lambda x: pdfs(x)[4]
  s=lambda x: pdfs(x)[5]
  sb=lambda x: pdfs(x)[6]
  c=lambda x: pdfs(x)[7]
  cb=lambda x: pdfs(x)[8]
  b=lambda x: pdfs(x)[9]
  bb=lambda x: pdfs(x)[10]

  print quad(lambda x: ub(x)-u(x),0,1)
  print quad(lambda x: x*(g(x)+u(x)+ub(x)+d(x)+db(x)+s(x)+sb(x)+c(x)+cb(x)+b(x)+bb(x)),0,1)
