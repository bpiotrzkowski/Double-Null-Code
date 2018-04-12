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

def boundary(scal,Es,bdytype,Nu,Nv,ru0,dr0v,du0,vmax,M0,Q,scalarfield):
    
    rnp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    phinp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    signp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    #drnp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    #dphinp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    #dsignp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    
    drnpu=np.zeros((Nu*scal*Es))
    drnpv=np.zeros((Nv*scal*Es))
    
    dsignpu=np.zeros((Nu*scal*Es))
    dsignpv=np.zeros((Nv*scal*Es))
    
    dphinpu=np.zeros((Nu*scal*Es))
    dphinpv=np.zeros((Nv*scal*Es))
    
    scalf=float(scal)

    rnp[0][0]=ru0 
    drnpv[0]=dr0v
    
    dt=du0/(scalf*Es)
    
    if scalarfield==True:
        A=.115
        #v=1.0
        v1=1.0
        v2=7.0
        v1n=int(v1*(Nv*scal*Es)/vmax)
        v2n=int(v2*(Nv*scal*Es)/vmax)
        
        
        for i in range(0,Nu*scal*Es):
            #dsignpu[i]=0.0
            phinp[i][0]=0.0
            dphinpu[i]=0.0
        for i in range(v1n,v2n):
                #dsignpv[i]=0.0
            v=i/(Nv*scal*Es)*vmax
            dphinpv[i]=192*A*(v-v1)**2.0*(v-v2)**2.0*(-2*v+v1+v2)/(v1-v2)**6.0
            phinp[0][i]=A*64*(v-v1)**3.0*(v2-v)**3.0/(v2-v1)**6.0

    if bdytype=="stan":
        sigu0=0.0
        sigv0=0.0
        dsignpu[0]=0.0
        
        drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v)*(Q**2.0/ru0**2.0-2*M0/ru0+1.0)
    
        dsignpv[0]=dsignpu[0]
        
        
        for i in range(0,Nv*scal*Es-1):
            dsignpv[i]=0.0
            rnp[0][i+1]=rnp[0][i]+dt*drnpv[i]
            drnpv[i+1]=drnpv[i]+dt*(drnpv[i]*dsignpv[i]-rnp[0][i]*dphinpv[i]**2.0)
            signp[0][i]=sigv0
        
        for i in range(0,Nu*scal*Es-1):
            dsignpu[i]=0.0
            rnp[i+1][0]=rnp[i][0]+dt*drnpu[i]
            drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnp[i][0]*dphinpu[i]**2.0)
            signp[i][0]=sigu0
            
    elif bdytype=="max" or bdytype=="hor":
        sigu0=0.0
        sigv0=0.0
        
        signp[0][0]=sigu0
        drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v)*(Q**2.0/ru0**2.0-2*M0/ru0+1.0)
        dsignpv[0]=0.0
        dsignpu[0]=0.0
        
        
        for i in range(0,Nv*scal*Es-1):
            dsignpv[i+1]=dsignpv[i]
            signp[0][i+1]=dsignpv[i]*dt+signp[0][i]
            rnp[0][i+1]=rnp[0][i]+dt*drnpv[i]
            drnpv[i+1]=drnpv[i]+dt*(drnpv[i]*dsignpv[i]-rnp[0][i]*dphinpv[i]**2.0)
            
        for i in range(0,Nu*scal*Es-1):
            dsignpu[i+1]=dsignpu[i]
            signp[i+1][0]=dsignpu[i+1]*dt+signp[i][0]
            rnp[i+1][0]=rnp[i][0]+dt*drnpu[i]
            drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnp[i][0]*dphinpu[i]**2.0)
            
    
    elif bdytype=="edd":
        
        sigu0=mth.log(1.0-2.0*M0/ru0+Q**2.0/ru0**2.0)
        sigv0=sigu0
        dsignpu[0]=2.0*(M0*ru0-Q**2.0)/(ru0*(Q**2.0+ru0*(-2*M0+ru0)))
        
        drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v*ru0**2.0)*(Q**2.0-2*M0*ru0+ru0**2.0)
        dsignpv[0]=dsignpu[0]
        
        for i in range(0,Nu*scal*Es-1):
            dsignpu[i]=2.0*(M0*rnp[i][0]-Q**2.0)/(rnp[i][0]*(Q**2.0+rnp[i][0]*(-2*M0+rnp[i][0])))*drnpu[i]
            rnp[i+1][0]=rnp[i][0]+dt*drnpu[i]
            drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnp[i][0]*dphinpu[i]**2.0)
        for i in range(0,Nv*scal*Es-1):
            rnp[0][i+1]=rnp[0][i]+dt*drnpv[i]
            drnpv[i+1]=drnpv[i]+dt*(drnpv[i]*dsignpv[i]-rnp[0][i]*dphinpv[i]**2.0)
        
    
        
           
    rnp=rnp[::scal,::scal]
    signp=signp[::scal,::scal]
    phinp=phinp[::scal,::scal]
    dphinpu=dphinpu[::scal]  
    dphinpv=dphinpv[::scal]
    dsignpu=dsignpu[::scal]
    dsignpv=dsignpv[::scal]
    drnpu=drnpu[::scal]
    drnpv=drnpv[::scal]
    
    return (rnp, signp, phinp, dsignpu,dphinpu)
    
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

        m=(1.0+4.0*mth.exp(-sig0)*dru0*drv0)*r0/2.0+(Q**2.0)/(2.0*r0)

        sig4=(sig3+sig2-sig1)+du0*dv0*(2.0*dru0*drv0/(r0)**2.0+mth.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2*dphiu0*dphiv0)
        #sig4=(sig3+sig2-sig1)+du0*dv0*(mth.exp(sig0)*(m/r0**3.0-3*Q**2/(2*r0**4.0))-2*dphiu0*dphiv0)
    
        
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
