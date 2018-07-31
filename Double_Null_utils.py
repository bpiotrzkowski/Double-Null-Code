# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:37:15 2018

@author: Brandon Piotrzkowski
"""

#################################
###Boundary Condition Function###

#################################
import numpy as np
import math as mth
#from mpmath import *
from decimal import *

def boundary(scal,Es,bdytype,Nu,Nv,ru0,dr0v,du0,vmax,M0,Q,scalarfield):
    
    rnpu=np.zeros((Nu*scal*Es))
    rnpv=np.zeros((Nv*scal*Es))
    
    signpu=np.zeros((Nu*scal*Es))
    signpv=np.zeros((Nv*scal*Es))
    
    phinpu=np.zeros((Nu*scal*Es))
    phinpv=np.zeros((Nv*scal*Es))
    
    drnpu=np.zeros((Nu*scal*Es))
    drnpv=np.zeros((Nv*scal*Es))
    
    dsignpu=np.zeros((Nu*scal*Es))
    dsignpv=np.zeros((Nv*scal*Es))
    
    dphinpu=np.zeros((Nu*scal*Es))
    dphinpv=np.zeros((Nv*scal*Es))
    
    scalf=float(scal)

    rnpu[0]=ru0
    rnpv[0]=ru0 
    drnpv[0]=dr0v
    
    dt=du0/(scalf*Es)
    
    if scalarfield==True:
        #A=.115
        #A=10.
        A=0.01
      
        v1=1.0
        v2=7.0
        v1n=int(v1*(Nv*scal*Es)/vmax)
        v2n=int(v2*(Nv*scal*Es)/vmax)
        
        
        for i in range(0,Nu*scal*Es):
            #dsignpu[i]=0.0
            phinpu[i]=0.0
            dphinpu[i]=0.0
        for i in range(v1n,v2n):
                #dsignpv[i]=0.0
            v=i/(Nv*scal*Es)*vmax
            dphinpv[i]=192*A*(v-v1)**2.0*(v-v2)**2.0*(-2*v+v1+v2)/(v1-v2)**6.0
            phinpv[i]=A*64*(v-v1)**3.0*(v2-v)**3.0/(v2-v1)**6.0

    if bdytype=="stan" or bdytype=="max" or bdytype=="hor":
        sigu0=0.0
        sigv0=0.0
        dsignpu[0]=0.0
        
        drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v)*(Q**2.0/ru0**2.0-2*M0/ru0+1.0)
    
        dsignpu[0]=dsignpu[0]
        
        
        for j in range(0,Nv*scal*Es-1):
            dsignpv[j]=0.0
            if rnpv[j]+dt*drnpv[j]>0.0:
                rnpv[j+1]=rnpv[j]+dt*drnpv[j]
                drnpv[j+1]=drnpv[j]+dt*(drnpv[j]*dsignpv[j]-rnpv[j]*dphinpv[j]**2.0)
                #print(dphinpv[j])
                #signpv[i]=sigv0
            else:
                break
        
        for i in range(0,Nu*scal*Es-1):
            dsignpu[i]=0.0
            if rnpu[i]+dt*drnpu[i]>0.0:
                rnpu[i+1]=rnpu[i]+dt*drnpu[i]
                drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnpu[i]*dphinpu[i]**2.0)
                #signpu[i]=sigu0
            else:
                break
                
            
    #elif bdytype=="max" or bdytype=="hor":
        #sigu0=0.0
        #sigv0=0.0
        
        #signpu[0]=sigu0
        #drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v)*(Q**2.0/ru0**2.0-2*M0/ru0+1.0)
        #dsignpv[0]=0.0
        #dsignpu[0]=0.0
        
        
        #for i in range(0,Nv*scal*Es-1):
            #dsignpv[i+1]=dsignpv[i]
            #signp[0][i+1]=dsignpv[i]*dt+signp[0][i]
            #rnp[0][i+1]=rnp[0][i]+dt*drnpv[i]
            #drnpv[i+1]=drnpv[i]+dt*(drnpv[i]*dsignpv[i]-rnp[0][i]*dphinpv[i]**2.0)
            
        #for i in range(0,Nu*scal*Es-1):
            #dsignpu[i+1]=dsignpu[i]
            #signp[i+1][0]=dsignpu[i+1]*dt+signp[i][0]
            #rnp[i+1][0]=rnp[i][0]+dt*drnpu[i]
            #drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnp[i][0]*dphinpu[i]**2.0)
            
    
    elif bdytype=="edd":
        
        signpu[0]=mth.log(1.0-2.0*M0/ru0+Q**2.0/ru0**2.0)
        signpv[0]=signpu[0]
        
     
        
        drnpu[0]=-mth.exp(signpu[0])/(4.0*dr0v)*(Q**2.0/ru0**2.0-2*M0/ru0+1.0)
        dsignpu[0]=2.0*(M0*ru0-Q**2.0)/(ru0*(Q**2.0+ru0*(-2*M0+ru0)))*drnpu[0]
        
        
        dsignpv[0]=0.0
        
        for i in range(0,Nu*scal*Es-1):
            rnpu[i+1]=rnpu[i]+dt*drnpu[i]
            dsignpu[i+1]=2.0*(M0*rnpu[i+1]-Q**2.0)/(rnpu[i+1]*(Q**2.0+rnpu[i+1]*(-2*M0+rnpu[i+1])))*drnpu[i]
            
            signpu[i+1]=signpu[i]+dt*dsignpu[i]
            drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnpu[i]*dphinpu[i]**2.0)
        for j in range(0,Nv*scal*Es-1):
            
            rnpv[j+1]=rnpv[j]+dt*drnpv[j]
            signpv[j+1]=signpv[j]+dt*dsignpv[j]
            drnpv[j+1]=drnpv[j]+dt*(drnpv[j]*dsignpv[j]-rnpv[j]*dphinpv[j]**2.0)
            
        
    
        
           
    rnpu=rnpu[::scal]
    rnpv=rnpv[::scal]
    signpu=signpu[::scal]
    signpv=signpv[::scal]
    phinpu=phinpu[::scal]
    phinpv=phinpv[::scal]
    #dphinpu=dphinpu[::scal]  
    #dphinpv=dphinpv[::scal]
    #dsignpu=dsignpu[::scal]
    #dsignpv=dsignpv[::scal]
    #drnpu=drnpu[::scal]
    #drnpv=drnpv[::scal]
    
    
    
    drnpu=None
    drnpv=None
    dsignpu=None
    dsignpv=None
    dphinpu=None
    dphinpv=None
    return (rnpu, rnpv, signpu, signpv, phinpu, phinpv) #dphinpu)#dsignpu, dphinpu)

###Alternate Boundary Condition Function (v-only)###

#################################

def boundaryv(scal,bdytype,Nv,ru0,dr0v,dv0,vmax,M0,Q,Lambda,scalarfield,datatype):
    
    
    if datatype==object:
        rnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        signpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        phinpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        drnpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dsignpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dphinpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        scald=Decimal(scal)
        dt=dv0/(scald)
    else:
        rnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        signpv=np.zeros((Nv*scal),dtype=datatype)
        phinpv=np.zeros((Nv*scal),dtype=datatype)
        drnpv=np.zeros((Nv*scal),dtype=datatype)
        dsignpv=np.zeros((Nv*scal),dtype=datatype)
        dphinpv=np.zeros((Nv*scal),dtype=datatype)
        scalf=float(scal)
        dt=dv0/(scalf)
    

    
    rnpv[0]=ru0 
    drnpv[0]=dr0v
    
   
    
    if scalarfield==True:
        A=.115
        #A=.2
      
        v1=1.0
        v2=7.0
        v1n=int(v1*(Nv*scal)/vmax)
        v2n=int(v2*(Nv*scal)/vmax)
        
       
        for i in range(v1n,v2n):
                #dsignpv[i]=0.0
            v=i/(Nv*scal)*vmax
            dphinpv[i]=192*A*(v-v1)**2.0*(v-v2)**2.0*(-2*v+v1+v2)/(v1-v2)**6.0
            phinpv[i]=A*64*(v-v1)**3.0*(v2-v)**3.0/(v2-v1)**6.0

    if bdytype=="stan" or bdytype=="max" or bdytype=="hor":
    
        sigv0=0.0
        if datatype==object:
            sigv0=Decimal(0)
        
        #dsignpv[0]=0.0
        
        
        for j in range(0,Nv*scal-1):
            dsignpv[j]=0.0
            if rnpv[j]+dt*drnpv[j]>0.0:
                rnpv[j+1]=rnpv[j]+dt*drnpv[j]
                drnpv[j+1]=drnpv[j]
                #print(dphinpv[j])
                signpv[j]=sigv0
            else:
                break
        
        print("Using Standard Coordinates Along V") 
            
    #######
    elif bdytype=="edd":
        
        signpv[0]=mth.log(1.0-2.0*M0/ru0+Q**2.0/ru0**2.0)
        
        dsignpv[0]=0.0
        drnpv[0]=dr0v
        
      
        for j in range(0,Nv*scal-1):
            dsignpv[j+1]=dsignpv[j]
            rnpv[j+1]=rnpv[j]+dt*drnpv[j]
            signpv[j+1]=signpv[j]+dt*dsignpv[j]
            drnpv[j+1]=drnpv[j]+dt*(drnpv[j]*dsignpv[j]-rnpv[j]*dphinpv[j]**2.0)
        print("Using Eddington Coordinates Along U")  
        
    ########  
    elif bdytype=="fulledd":
        
        signpv[0]=np.log(1.0-2.0*M0/ru0+Q**2.0/ru0**2.0-Lambda*ru0**2/3)
        
        dsignpv[0]=2*(3*Q**2-3*M0*ru0+ru0**4*Lambda)/(ru0*(-3*Q**2+ru0*(6*M0-3*ru0+ru0**3*Lambda)))*dr0v
        drnpv[0]=dr0v
        
      
        for j in range(0,Nv*scal-1):
            rnpv[j+1]=rnpv[j]+dt*drnpv[j]
            #signpv[j+1]=signpv[j]+dt*dsignpv[j]
            signpv[j+1]=np.log(1.0-2.0*M0/rnpv[j+1]+Q**2.0/rnpv[j+1]**2.0)
            drnpv[j+1]=drnpv[j]+dt*(drnpv[j]*dsignpv[j]-rnpv[j]*dphinpv[j]**2.0)  
            dsignpv[j+1]=2*(3*Q**2-3*M0*rnpv[j+1]+rnpv[j+1]**4*Lambda)/(rnpv[j+1]*(-3*Q**2+rnpv[j+1]*(6*M0-3*rnpv[j+1]+rnpv[j+1]**3*Lambda)))*drnpv[j+1]
            #print(dsignpv[j+1])
            
            
        print("Using Full Eddington Coordinates")
        
           
    
    rnpv=rnpv[::scal]
    
    signpv=signpv[::scal]
    
    phinpv=phinpv[::scal]
    #dphinpu=dphinpu[::scal]  
    #dphinpv=dphinpv[::scal]
    #dsignpu=dsignpu[::scal]
    #dsignpv=dsignpv[::scal]
    #drnpu=drnpu[::scal]
    #drnpv=drnpv[::scal]
    
    
    
    drnpu=None
    drnpv=None
    dsignpu=None
    dsignpv=None
    dphinpu=None
    dphinpv=None
    return (rnpv, signpv, phinpv) #dphinpu)#dsignpu, dphinpu)    
    
#############################################
###Defining Propagation Algorithm/Function###
def x4giver(un,vn,Es,k,du00,rnpf,phinpf,signpf,Q):
    
    
    Es=float(Es)
    du0=du00/(Es)
    dv0=du0
    #dv0
    
    
    r1=rnpf[un][vn][k]
    r2=rnpf[un][vn+1][k]
    r3=rnpf[un+1][vn][k]
    phi1=phinpf[un][vn][k]
    phi2=phinpf[un][vn+1][k]
    phi3=phinpf[un+1][vn][k]
    sig1=signpf[un][vn][k]
    sig2=signpf[un][vn+1][k]
    sig3=signpf[un+1][vn][k]
    
    
    r0=(r2+r3)/2.0
    
    sig0=(sig2+sig3)/2.0
    
    try:
        r4i=(r3+r2-r1)-(r3-r1)*(r2-r1)/r0-du0*dv0*mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
        r4=(r3+r2-r1)-(r3-r1+r4i-r2)*(r2-r1+r4i-r3)/(4*r0)-du0*dv0*mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
    
        dru0=(r3-r1+r4-r2)/(2.0*du0)
        drv0=(r2-r1+r4-r3)/(2.0*dv0)
        
        r0=(r1+r2+r3+r4)/4.0

        phi4i=(phi3+phi2-phi1)-du0*dv0/r0*(dru0*(phi2-phi1)/dv0+drv0*(phi3-phi1)/du0)
        phi4=(phi3+phi2-phi1)-du0*dv0/r0*(dru0*(phi2-phi1+phi4i-phi3)/(2*dv0)+drv0*(phi3-phi1+phi4i-phi2)/(2*du0))
    
        dphiu0=(phi3-phi1+phi4-phi2)/(2.0*du0)
        dphiv0=(phi2-phi1+phi4-phi3)/(2.0*dv0)

        sig4=(sig3+sig2-sig1)+du0*dv0*(2.0*dru0*drv0/(r0)**2.0+mth.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2*dphiu0*dphiv0)
        sig4=(sig3+sig2-sig1)+du0*dv0*(mth.exp(sig0)*(m/r0**3.0-3*Q**2/(2*r0**4.0))-2*dphiu0*dphiv0)
        
        sig0=(sig1+sig2+sig3+sig4)/4.0
    
        m=(1.0+4.0*mth.exp(-sig0)*dru0*drv0)*r0/2.0+(Q**2.0)/(2.0*r0)
        
        dsigu0=(sig3-sig1+sig4-sig2)/(2.0*du0)
        #dsigv0=(sig2-sig1+sig4-sig3)/(2.0*dv0)
    
        #druv=-(r3-r1+r4-r2)*(r2-r1+r4-r3)/(4*r0*du0*dv0)-mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
        #dphiuv=-1/r0*(dru0*(phi2-phi1+phi4-phi3)/(2*dv0)+drv0*(phi3-phi1+phi4-phi2)/(2*du0))
        #dsiguv=2.0*dru0*drv0/(r0)**2.0+mth.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2*dphiu0*dphiv0
    except OverflowError:
        
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
        
    if r4<.001:
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
        
    answer=[r4,phi4,sig4,m,dru0,dsigu0]
    
    return answer

#############################################
###Defining Alternate Propagation Algorithm/Function###
def x4giveralt(un,vn,du0,dv0,rnpf,phinpf,signpf,massnpf,M0,Q,Lambda,datatype):
    
    
    #Es=float(Es)
    #du0=du00#/(Es)
    #dv0=dv00#/(Es)
    #dv0
    
    
    r1=rnpf[un][vn]
    r2=rnpf[un][vn+1]
    r3=rnpf[un+1][vn]
    phi1=phinpf[un][vn]
    phi2=phinpf[un][vn+1]
    phi3=phinpf[un+1][vn]
    sig1=signpf[un][vn]
    sig2=signpf[un][vn+1]
    sig3=signpf[un+1][vn]
    m1=massnpf[un][vn]
    m2=massnpf[un][vn+1]
    m3=massnpf[un+1][vn]
   
    
    
    r0=(r2+r3)/2.0
    dru0=(r3-r1)/du0
    drv0=(r2-r1)/dv0
    dphiu0=(phi3-phi1)/du0
    dphiv0=(phi2-phi1)/dv0
    sig0=(sig2+sig3)/2.0
    m0=(m2+m3)/2.0
    
    try:
        #r4i=(r3+r2-r1)-(r3-r1)*(r2-r1)/(r0)-du0*dv0*np.exp(sig0)/(4*r0)*(1.0-(Q**2.0)/(r0**2.0)-Lambda*r0**2.0)
        r4i=(r3+r2-r1)+du0*dv0*np.exp(sig0)*(Q**2.0-M0*r0)/(2.0*r0**3.0)
        
        r0=(r1+r2+r3+r4i)/4.0
        dru0=(r3-r1+r4i-r2)/(2.0*du0)
        drv0=(r2-r1+r4i-r3)/(2.0*dv0)
        
        phi4i=(phi3+phi2-phi1)-du0*dv0/r0*(dru0*(phi2-phi1)/dv0+drv0*(phi3-phi1)/du0)
        dphiu0=(phi3-phi1+phi4i-phi2)/(2.0*du0)
        dphiv0=(phi2-phi1+phi4i-phi3)/(2.0*dv0)
        
        #sig4i=(sig3+sig2-sig1)+2.0*(r3-r1+r4i-r2)*(r2-r1+r4i-r3)/(4.0*(r0)**2.0)+du0*dv0*(np.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2.0*dphiu0*dphiv0)
        sig4i=(sig3+sig2-sig1)+du0*dv0*np.exp(sig0)*(2.0*M0*r0-3.0*Q**2.0)/(2.0*r0**4.0)
        
        sig0=(sig1+sig2+sig3+sig4i)/4.0
           
        #r4=(r3+r2-r1)-(r3-r1+r4i-r2)*(r2-r1+r4i-r3)/(4.0*r0)-du0*dv0*np.exp(sig0)/(4.0*r0)*(1.0-(Q**2.0)/(r0**2.0)-Lambda*r0**2.0)
        r4=(r3+r2-r1)+du0*dv0*np.exp(sig0)*(Q**2.0-M0*r0)/(2.0*r0**3.0)
        
        
        r0=(r1+r2+r3+r4)/4.0
        dru0=(r3-r1+r4-r2)/(2.0*du0)
        drv0=(r2-r1+r4-r3)/(2.0*dv0)
        
        
        phi4=(phi3+phi2-phi1)-du0*dv0/r0*(dru0*(phi2-phi1+phi4i-phi3)/(2*dv0)+drv0*(phi3-phi1+phi4i-phi2)/(2*du0))
    
        dphiu0=(phi3-phi1+phi4-phi2)/(2.0*du0)
        dphiv0=(phi2-phi1+phi4-phi3)/(2.0*dv0)

        #sig4=(sig3+sig2-sig1)+2.0*(r3-r1+r4-r2)*(r2-r1+r4-r3)/(4.0*(r0)**2.0)+du0*dv0*(np.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2.0*dphiu0*dphiv0)
        sig4=(sig3+sig2-sig1)+du0*dv0*np.exp(sig0)*(2.0*M0*r0-3.0*Q**2.0)/(2.0*r0**4.0)
            
           
        sig0=(sig1+sig2+sig3+sig4)/4.0
        
         
        
        #m4=(m3+m2-m1)+du0*dv0*(2*r0**3*np.exp(-sig0)*(dphiu0*dphiv0)**2.0-r0*(1-2*m0/r0+Q**2/r0**2)*dphiu0*dphiv0)
        
        #dru0=(r3-r1)/du0-(r3-r1)*(r2-r1)/(2*r0*du0)-dv0*np.exp(sig0)/(8*r0)*(1-Q**2/r0**2)
        #drv0=(r2-r1)/dv0-(r3-r1)*(r2-r1)/(2*r0*dv0)-du0*np.exp(sig0)/(8*r0)*(1-Q**2/r0**2)
        
        m4=(1.0+4.0*mth.exp(-sig0)*dru0*drv0)*r0/2.0+(Q**2.0)/(2.0*r0)
        
        #dru0=(r3-r1)/du0
        #drv0=(r2-r1)/dv0
        
        #m2=(1.0+4.0*np.exp(-sig1)*dru0*drv0)*r1/2.0+(Q**2.0)/(2.0*r1)
        
        #m3=-(r4-r3-r2+r1)/(dv0*du0)*np.exp(-sig0)*2*r0**2+Q**2/r0
        #m4=(sig4-sig3-sig2+sig1)/(dv0*du0)*np.exp(-sig0)*r0**3+3/2*Q**2/r0
        
        
        
        dsigu0=(sig3-sig1+sig4-sig2)/(2.0*du0)
        dsigv0=(sig2-sig1+sig4-sig3)/(2.0*dv0)
        
        exterm1=dru0
        exterm2=drv0
        exterm3=dsigu0
        
        #exterm1=(du0*dv0*np.exp(sig0)*(Q**2.0-M0*r0)/(2.0*r0**3.0)+1.0)-1.0
        #exterm2=(du0*dv0*np.exp(sig0)*(2.0*M0*r0-3.0*Q**2.0)/(2.0*r0**4.0)+1.0)-1.0
        #exterm3=(du0*dv0/r0*(dru0*(phi2-phi1+phi4i-phi3)/(2*dv0)+drv0*(phi3-phi1+phi4i-phi2)/(2*du0))+1.0)-1.0
        
    except OverflowError or RuntimeWarning:
        
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
        
    if r4<.001:
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
    
    #mf=M0+min(abs(m-M0),abs(m2-M0),abs(m3-M0),abs(m4-M0))
    
    #answer=[r4,phi4,sig4,m,m2,m3,m4,dru0,dsigu0]
    #answer=np.array([r4,phi4,sig4,m4,dru0,dsigu0],dtype=datatype)
    answer=np.array([r4,phi4,sig4,m4,dru0,dsigu0,exterm1,exterm2,exterm3],dtype=datatype)
    
    return answer

#############################################
###Defining Alternate Propagation Algorithm/Function###
def x4giveraltd(un,vn,du0,dv0,rnpf,phinpf,signpf,massnpf,M0,Q,Lambda,datatype):
    
    getcontext().prec=40
    
   
    
    
    r1=rnpf[un][vn]
    r2=rnpf[un][vn+1]
    r3=rnpf[un+1][vn]
    phi1=phinpf[un][vn]
    phi2=phinpf[un][vn+1]
    phi3=phinpf[un+1][vn]
    sig1=signpf[un][vn]
    sig2=signpf[un][vn+1]
    sig3=signpf[un+1][vn]
    m1=massnpf[un][vn]
    m2=massnpf[un][vn+1]
    m3=massnpf[un+1][vn]
   
    
    
    r0=(r2+r3)/Decimal(2)
    dru0=(r3-r1)/du0
    drv0=(r2-r1)/dv0
    #print(type(phi3),type(phi1),type(du0))
    dphiu0=(phi3-phi1)/du0
    dphiv0=(phi2-phi1)/dv0
    sig0=(sig2+sig3)/Decimal(2.0)
    m0=(m2+m3)/Decimal(2.0)
    
    try:
        #r4i=(r3+r2-r1)-(r3-r1)*(r2-r1)/(r0)-du0*dv0*(sig0).exp()/(Decimal(4)*r0)*(Decimal(1)-(Q**Decimal(2.0))/(r0**Decimal(2.0))-Lambda*r0**Decimal(2))
        
        r4i=(r3+r2-r1)+du0*dv0*(sig0).exp()*(Q**Decimal(2)-M0*r0)/(Decimal(2)*r0**Decimal(3))
        
        r0=(r1+r2+r3+r4i)/Decimal(4.0)
        dru0=(r3-r1+r4i-r2)/(Decimal(2.0)*du0)
        drv0=(r2-r1+r4i-r3)/(Decimal(2.0)*dv0)
        
        phi4i=(phi3+phi2-phi1)-1/r0*((r3-r1)*(phi2-phi1)+(r2-r1)*(phi3-phi1))
        
        dphiu0=(phi3-phi1+phi4i-phi2)/(Decimal(2.0)*du0)
        dphiv0=(phi2-phi1+phi4i-phi3)/(Decimal(2.0)*dv0)
        
        #sig4i=(sig3+sig2-sig1)+du0*dv0*(Decimal(2.0)*dru0*drv0/(r0)**Decimal(2.0)+(sig0).exp()/(Decimal(2.0)*r0**Decimal(2.0))*(Decimal(1.0)-Decimal(2.0)*(Q**Decimal(2.0))/(r0**Decimal(2.0)))-Decimal(2)*dphiu0*dphiv0)
        sig4i=(sig3+sig2-sig1)+du0*dv0*(sig0).exp()*(Decimal(2)*M0*r0-Decimal(3)*Q**Decimal(2))/(Decimal(2)*r0**Decimal(4))
        
        sig0=(sig1+sig2+sig3+sig4i)/Decimal(4.0)
           
        #r4=(r3+r2-r1)-(r3-r1+r4i-r2)*(r2-r1+r4i-r3)/(Decimal(4)*r0)-du0*dv0*(sig0).exp()/(Decimal(4)*r0)*(Decimal(1)-(Q**Decimal(2.0))/(r0**Decimal(2.0))-Lambda*r0**Decimal(2))
        r4=(r3+r2-r1)+du0*dv0*(sig0).exp()*(Q**Decimal(2)-M0*r0)/(Decimal(2)*r0**Decimal(3))
        
        
        r0=(r1+r2+r3+r4)/Decimal(4.0)
        dru0=(r3-r1+r4-r2)/(Decimal(2.0)*du0)
        drv0=(r2-r1+r4-r3)/(Decimal(2.0)*dv0)
        
        
        phi4=(phi3+phi2-phi1)-1/r0*((r3-r1+r4-r2)*(phi2-phi1+phi4i-phi3)/Decimal(2)+(r2-r1+r4-r3)*(phi3-phi1+phi4i-phi2)/Decimal(2))
    
        dphiu0=(phi3-phi1+phi4-phi2)/(Decimal(2.0)*du0)
        dphiv0=(phi2-phi1+phi4-phi3)/(Decimal(2.0)*dv0)

        #sig4=(sig3+sig2-sig1)+du0*dv0*((Decimal(2.0)*dru0*drv0/(r0)**Decimal(2.0)+(sig0).exp()/(Decimal(2.0)*r0**Decimal(2.0))*(Decimal(1.0)-Decimal(2.0)*(Q**Decimal(2.0))/(r0**Decimal(2.0)))-Decimal(2)*dphiu0*dphiv0))
        sig4=(sig3+sig2-sig1)+du0*dv0*(sig0).exp()*(Decimal(2)*M0*r0-Decimal(3)*Q**Decimal(2))/(Decimal(2)*r0**Decimal(4))
            
           
        sig0=(sig1+sig2+sig3+sig4)/Decimal(4.0)
        
        #print(phi4i,phi4)
        
        #m4=(m3+m2-m1)+du0*dv0*(Decimal(2)*r0**Decimal(3)*(-sig0).exp()*(dphiu0*dphiu0)**Decimal(2.0)-r0*(Decimal(1)-Decimal(2)*m0/r0+Q**Decimal(2)/r0**Decimal(2))*dphiu0*dphiu0)
        
        m4=(Decimal(1.0)+Decimal(4.0)*(-sig0).exp()*dru0*drv0)*r0/Decimal(2.0)+(Q**Decimal(2.0))/(Decimal(2.0)*r0)
        
        #######
        #dru0=(r3-r1)/du0
        #drv0=(r2-r1)/dv0
        
        #m2=(1.0+4.0*np.exp(-sig1)*dru0*drv0)*r1/2.0+(Q**2.0)/(2.0*r1)
        
        #m3=-(r4-r3-r2+r1)/(dv0*du0)*np.exp(-sig0)*2*r0**2+Q**2/r0
        #m4=(sig4-sig3-sig2+sig1)/(dv0*du0)*np.exp(-sig0)*r0**3+3/2*Q**2/r0
        
        dsigu0=(sig3-sig1+sig4-sig2)/(Decimal(2.0)*du0)
        
        #exterm1=dru0
        #exterm2=(sig0).exp()
        #exterm3=dsigu0
        
        exterm1=(du0*dv0*np.exp(sig0)*(Q**Decimal(2.0)-M0*r0)/(Decimal(2.0)*r0**Decimal(3.0))+Decimal(1.0))-Decimal(1.0)
        exterm2=(du0*dv0*(sig0).exp()*(Decimal(2)*M0*r0-Decimal(3)*Q**Decimal(2))/(Decimal(2)*r0**Decimal(4))+Decimal(1.0))-Decimal(1.0)
        exterm3=(du0*dv0/r0*(dru0*(phi2-phi1+phi4i-phi3)/(Decimal(2)*dv0)+drv0*(phi3-phi1+phi4i-phi2)/(Decimal(2)*du0))+Decimal(1.0))-Decimal(1.0)
        
    except OverflowError or RuntimeWarning:
        
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
        
    if r4<.001:
        r4i=np.nan
        r4=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi4i=np.nan
        phi4=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig4=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        #dsigv0=np.nan
    
        #druv=np.nan
        #dphiuv=np.nan
        #dsiguv=np.nan
    
    #mf=M0+min(abs(m-M0),abs(m2-M0),abs(m3-M0),abs(m4-M0))
    
    #answer=[r4,phi4,sig4,m,m2,m3,m4,dru0,dsigu0]
    answer=np.array([r4,phi4,sig4,m4,dru0,dsigu0,exterm1,exterm2,exterm3],dtype=datatype)
    
    return answer
