import numpy as np
import math as mth

import scipy.special.lambertw as lambertw
from scipy.special import gammainc as gamma

import scipy
from pynverse import inversefunc
from mpmath import gammainc as gamma2
from scipy.special import gamma as gamma
from scipy.special import gammaincc as gamma3
import scipy.integrate as integrate
#from mpmath import fp as float2


def rbounsch(m,ru0,dr0,x):
    try:
        #r=2*m*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)*(1/(2*m)+dr0*x/(4*m)*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m),k=0))/lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0),k=0)),k=0))
        C1=-dr0*mth.exp(ru0/(2*m))*m*(2*m-ru0)*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m)))/lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))
        C2=2*lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))/(dr0*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))))
        r=2*m*(1+lambertw(C1*(x+C2)/(4*mth.exp(1)*m**2)))
    except:
        r=np.nan
    return r.real


def rbounschinv(m,ru0,dr0,r):
    try:
        C1=-dr0*mth.exp(ru0/(2*m))*m*(2*m-ru0)*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m)))/lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))
        C2=2*lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))/(dr0*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))))
    
        #r=r*(32*m**3.0/r*mth.exp(-r/(2*m)))**(1.0)
    
        x=(2*m*mth.exp(r/(2*m))*(r-2*m)-C1*C2)/C1
    except:
        x=np.nan
    return x.real


def drbounsch(m,ru0,dr0,x):
    C1=-dr0*mth.exp(ru0/(2*m))*m*(2*m-ru0)*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m)))/lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))
    #C2=-2*mth.exp(ru0/(2*m))*m*(2*m-ru0)/C1
    C2=2*lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))/(dr0*(1+lambertw(-mth.exp(-1+ru0/(2*m))*(2*m-ru0)/(2*m))))
    drdx=2*m*lambertw(C1*(x+C2)/(4*mth.exp(1)*m**2))/((x+C2)*(1+lambertw(C1*(x+C2)/(4*mth.exp(1)*m**2))))*(32*m**(3.0)/ru0*mth.exp(-ru0/(2.0*m)))**(-1.0)
    return drdx.real


def dr(m,dru,r):
    drv=-1/4/(dru)*(1-2*m/r)*32*m**3/r*mth.exp(-r/(2*m))
    return drv

def dr2(m,dru,r): 
    drv=-1/4/(dru)*(1-2*m/r)
    return drv

def rsch(m,ru0,dr0,u,v):
    #dr02=dr0*(32*m**(3.0)/ru0*mth.exp(-ru0/(2.0*m)))**(1.0)
    r=rbounsch(m,rbounsch(m,ru0,dr(m,dr0,ru0),u),dr(m,drbounsch(m,ru0,dr(m,dr0,ru0),u),rbounsch(m,ru0,dr(m,dr0,ru0),u)),v)
                                                                       
    return r
                                                                       
     
def esig(m,r):
    sig=32*m**3/r*mth.exp(-r/(2*m))**(1.0)
    
    return sig
    
def rsch_to_file(m,ru0,dr0,u,v):
    dr02=dr0*(32*m**(3.0)/ru0*mth.exp(-ru0/(2.0*m)))**(1.0)
    r=rbounsch(m,rbounsch(m,ru0,dr(m,dr0,ru0),u),dr(m,drbounsch(m,ru0,dr(m,dr0,ru0),u),rbounsch(m,ru0,dr(m,dr0,ru0),u)),v)
    #r=rbounsch(m,rbounsch(m,ru0,-.375,u),dr(m,drbounsch(m,ru0,dr02,u),rbounsch(m,ru0,-.375,u)),v)
                                                                         
    with open('rsch_output.txt','a') as f:
        f.write(r)

######################################################################################################################################################################################

###Reissner Nordstrom Functions###

def esigrn(m,Q,r):
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
    
    sig=rplus*rminus/(kplus**2*r**2)*mth.exp(-2*kplus*r)*(r/rminus-1)**(1+kplus/kminus)
    
    return sig.real
    
def drrn(m,Q,dru,r): 
    drv=-1/4/(dru)*(1-2*m/r+Q**2/r**2)*esigrn(m,Q,r)
    return drv

def dr2rn(m,Q,dru,r): 
    drv=-1/4/(dru)*(1-2*m/r+Q**2/r**2)
    return drv

#def Gamma(s,x):
    #def real_func(t,s):
        #return scipy.real(t**(s-1)*mth.exp(-t))
    #def imag_func(t,s):
        #return scipy.imag(t**(s-1)*mth.exp(-t))
    #real_integral = integrate.quad(real_func, x, np.inf,args=s)
    #imag_integral = integrate.quad(imag_func, x, np.inf,args=s)
    
    
    #print(real_integral[1]+1j*imag_integral[1])
    #return real_integral[0]+1j*imag_integral[0]
    

#def Gamma2(x,s):
    
    
    #try:
        #gval=(gamma3(s+1,x)*gamma(s+1)-(x)**(s)*mth.exp(-x))/s
    #except RuntimeWarning:
        #gval=np.nan
    #return gval

def Gamma(s,x):
    
    if type(x) is np.ndarray:
        x=x[0]
    #print(type(x))
    #print(s,x)
    gval=float(gamma2(s,x).real)+1j*float(gamma2(s,x).imag)
    #gval=gamma2(s,x)
    #print(gval)
    #print(gamma2(s,x))
    return gval

def drbounrn(r,u,m,Q,ru0,dr0):
    
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
        
    return 1/r**2*dr0*mth.exp(2*kplus*(ru0-r))*(r-rminus)**(1+kplus/kminus)*(ru0-rminus)**(-1-kplus/kminus)*ru0**2

    
def rbounrn(m,Q,ru0,dr0,u):
    
    return scipy.integrate.odeint(drbounrn,ru0,u,args=(m,Q,ru0,dr0))
  
def rbounrninv(r,m,Q,ru0,dr0):
   
    #m=1.0
    #Q=.95
    #ru0=5.0
    #dr0=.4
    
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
    

    #try:
    C1=mth.exp(2*ru0*kplus)*ru0**2*(ru0-rminus)**(-1-(kplus)/kminus)*dr0
    
    C2=-1/(kplus**2*C1)*2**(-2+kplus/kminus)*mth.exp(2*rminus*kplus)*(ru0-rminus)**(-kplus/kminus)*((rminus-ru0)*kplus)**(kplus/kminus)*(-4*rminus*kplus*Gamma(1-kplus/kminus,-2*(ru0-rminus)*kplus)+Gamma(2-kplus/kminus,-2*(ru0-rminus)*kplus)+4*(rminus*kplus)**2*Gamma(-kplus/kminus,-2*(ru0-rminus)*kplus))
        
    x=-C2-1/(kplus**2*C1)*2**(-2+kplus/kminus)*mth.exp(2*rminus*kplus)*(r-rminus)**(-kplus/kminus)*((rminus-r)*kplus)**(kplus/kminus)*(-4*rminus*kplus*Gamma(1-kplus/kminus,-2*(r-rminus)*kplus)+Gamma(2-kplus/kminus,-2*(r-rminus)*kplus)+4*(rminus*kplus)**2*Gamma(-kplus/kminus,-2*(r-rminus)*kplus))
    
    #except:
        #x=np.nan
    #return float(gamma2(-kplus/kminus,2*(ru0-rminus)*kplus))

    return x.real
    
def rbounrninv2(r,m,Q,ru0,dr0,x):   

    return rbounrninv(r,m,Q,ru0,dr0)-x
    
def rrn(m,Q,ru0,dr0,u,v,i):
        
    #rplus=m+(m**2-Q**2)**(0.5)
    #rminus=m-(m**2-Q**2)**(0.5)
    #kplus=(rplus-rminus)/(2*(rplus)**2)
    #kminus=abs((rminus-rplus)/(2*(rminus)**2))
      
    try:
        ru=rbounrn(m,Q,ru0,drrn(m,Q,dr0,ru0),u)[i]
        ru=float(ru)
        #print(ru)
        drdu=drbounrn(ru,u,m,Q,ru0,dr2rn(m,Q,dr0,ru0))
        drdu=float(drdu)
        #print(drdu)
        rv=rbounrn(m,Q,ru,drrn(m,Q,drdu,ru),v)
        
  
        
        
        
        
        #ru= inversefunc(rbounrninv,args=(m,Q,ru0,drrn(m,Q,dr0,ru0)),y_values=u,domain=[0.68,5.1],open_domain=[False,False],image=[0,13.0],accuracy=3)
    
        #drdu= esigrn(m,Q,ru0)**(0)*drrn(m,Q,dr0,ru0)*mth.exp(2*kplus*(ru-ru0))*(ru-rminus)**(1+kplus/kminus)*(ru0-rminus)**(-kplus/kminus)*ru0**(2)/(ru**2*(ru0-rminus))
    
       
        #print(ru,drdu)
        
        #rv= inversefunc(rbounrninv,args=(m,Q,ru,drrn(m,Q,drdu,ru)*esigrn(m,Q,ru0)**(0)),y_values=v,domain=[.68,20],open_domain=[False,False],image=[0,1000],accuracy=3)
        
        
        
        
        ####
        
        #rv= inversefunc(rbounrninv,args=(m,Q,ru0,dr0*esigrn(m,Q,ru0)),y_values=v,domain=[0.68,20.0],open_domain=[False,False],image=[0,13.0],accuracy=1)
    
        #drdv=dr0*esigrn(m,Q,ru0)**(1)*mth.exp(2*kplus*(rv-ru0))*(rv-rminus)**(1+kplus/kminus)*(ru0-rminus)**(-kplus/kminus)*ru0**(2)/(rv**2*(ru0-rminus))
    
        #ru= inversefunc(rbounrninv,args=(m,Q,rv,dr2rn(m,Q,drdv,rv)),y_values=u,domain=[.68,20.0],open_domain=[False,False],image=[0,1000],accuracy=1)
    
    except ValueError:
        rv=np.nan
    
    #print(rv,drdu)
    return rv

#def rrn2(m,Q,ru0,dr0,u,v):
    
    #rplus=m+(m**2-Q**2)**(0.5)
    #rminus=m-(m**2-Q**2)**(0.5)
    #kplus=(rplus-rminus)/(2*(rplus)**2)
    #kminus=abs((rminus-rplus)/(2*(rminus)**2))
    
    #ru = scipy.optimize.fsolve(rbounrninv2,ru0,args=(m,Q,ru0-drrn(m,Q,dr0,ru0)*u,drrn(m,Q,dr0,ru0),u))
                               
    #drdu= esigrn(m,Q,ru0)**(0)*dr2rn(m,Q,dr0,ru0)*mth.exp(2*kplus*(ru-ru0))*(ru-rminus)**(1+kplus/kminus)*(ru0-rminus)**(-kplus/kminus)*ru0**(2)/(ru**2*(ru0-rminus))
                               
    #rv= scipy.optimize.fsolve(rbounrninv2,ru0,args=(m,Q,ru,dr2rn(m,Q,drdu,ru)*esigrn(m,Q,ru0)**(0),v))
                               
    #return ru
    
##############################################################

###Reissner Nordstrom Interior Functions###    

def esigrni(m,Q,r):
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
    
    sig=rplus*rminus/(kminus**2*r**2)*mth.exp(-2*kminus*r)*(1-r/rplus)**(1+kminus/kplus)
    
    return sig.real
    
def drrni(m,Q,dru,r): 
    drv=-1/4/(dru)*(1-2*m/r+Q**2/r**2)*esigrni(m,Q,r)
    return drv

def drbounrni(r,u,m,Q,ru0,dr0):
    
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
        
    return 1/r**2*dr0*mth.exp(2*kminus*(r-ru0))*(r-rplus)**(-1-kminus/kplus)*(ru0-rminus)**(1+kminus/kplus)*ru0**2

    
def rbounrni(m,Q,ru0,dr0,u):
    
    return scipy.integrate.odeint(drbounrni,ru0,u,args=(m,Q,ru0,dr0))
                              
def rbounrninvi(r,m,Q,ru0,dr0):
   
    #m=1.0
    #Q=.95
    #ru0=5.0
    #dr0=.4
    
    rplus=m+(m**2-Q**2)**(0.5)
    rminus=m-(m**2-Q**2)**(0.5)
    kplus=(rplus-rminus)/(2*(rplus)**2)
    kminus=abs((rminus-rplus)/(2*(rminus)**2))
    

    #try:
    C1=mth.exp(2*ru0*kplus)*ru0**2*(ru0-rminus)**(-1-(kplus)/kminus)*dr0
    
    C2=-1/(C1)*2**(-4-kminus/kplus)*(kminus)**(-4-kminus/kplus)*(4*rplus*kminus*Gamma(3+kminus/kplus,2*(ru0-rplus)*kminus)+Gamma(4+kminus/kplus,2*(ru0-rplus)*kminus)+4*(rplus*kminus)**2*Gamma(2+kminus/kplus,2*(ru0-rplus)*kminus))
        
    x=-C2-1/(C1)*2**(-4-kminus/kplus)*(kminus)**(-4-kminus/kplus)*(4*rplus*kminus*Gamma(3+kminus/kplus,2*(r-rplus)*kminus)+Gamma(4+kminus/kplus,2*(r-rplus)*kminus)+4*(rplus*kminus)**2*Gamma(2+kminus/kplus,2*(r-rplus)*kminus))

    return x.real
     
def rrni(m,Q,ru0,dr0,u,v,i):
        
    #rplus=m+(m**2-Q**2)**(0.5)
    #rminus=m-(m**2-Q**2)**(0.5)
    #kplus=(rplus-rminus)/(2*(rplus)**2)
    #kminus=abs((rminus-rplus)/(2*(rminus)**2))
      
    try:
        ru=rbounrni(m,Q,ru0,drrni(m,Q,dr0,ru0),u)[i]
        ru=float(ru.real)
        #print(ru)
        drdu=drbounrni(ru,u,m,Q,ru0,dr2rn(m,Q,dr0,ru0))
        drdu=float(drdu.real)
        #print(drdu)
        rv=rbounrni(m,Q,ru,drrni(m,Q,drdu,ru),v)
        
    
    except ValueError:
        rv=np.nan
    
    #print(rv,drdu)
    return rv