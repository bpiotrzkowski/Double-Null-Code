# -*- coding: utf-8 -*-
"""
Created on Wed Sep 2 21:21:15 2018

@author: Brandon Piotrzkowski
"""

#################################
###Boundary Condition Function###

#################################
import numpy as np
import math as mth
#from mpmath import *
from decimal import *
from scipy import interpolate
from scipy import optimize

def rc(M,Q,Lambda):
    def fr(r,M,Q,Lambda):
        return 1-2*M/r+Q**2.0/r**2.0-Lambda*r**2.0/3.0
    if Lambda>0 or Lambda<0:
        sol = optimize.root(fr,[50.0],args=(M,Q,Lambda), method='hybr')
        #rminus=sol.x[0]
        #rplus=sol.x[1]
        rcosm=sol.x[0]
    else:
        rcosm=np.inf
    return rcosm

def boundaryv(scal,bdytype,Nv,ru0,dr0v,dv0,vmax,M0,Q,Lambda,scalarfield,A,rcosmtol,datatype):
    
    
    if datatype==object:
        rnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        signpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        phinpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        drnpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dsignpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dphinpv=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        drnpu=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dsignpu=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        dphinpu=np.zeros((Nv*scal),dtype=datatype)*Decimal(0)
        scald=Decimal(scal)
        dt=dv0/(scald)
    else:
        rnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        signpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        phinpv=np.zeros((Nv*scal),dtype=datatype)#*np.nan
        drnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        dsignpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        dphinpv=np.zeros((Nv*scal),dtype=datatype)#*np.nan
        drnpu=np.zeros((Nv*scal),dtype=datatype)*np.nan
        dsignpu=np.zeros((Nv*scal),dtype=datatype)*np.nan
        dphinpu=np.zeros((Nv*scal),dtype=datatype)*np.nan
        massnpv=np.zeros((Nv*scal),dtype=datatype)*np.nan
        scalf=float(scal)
        dt=dv0/(scalf)
    

    
    rnpv[0]=ru0 
    drnpv[0]=dr0v
    
   
    
    if scalarfield==True:
        #A=.115
        #A=.2
      
        v1=1.0
        v2=4.0
        v1n=int(v1*(Nv*scal)/vmax)
        v2n=int(v2*(Nv*scal)/vmax)
        
       
        for j in range(v1n,v2n):
                #dsignpv[i]=0.0
            v=j/(Nv*scal)*vmax
            dphinpv[j]=192*A*(v-v1)**2.0*(v-v2)**2.0*(-2*v+v1+v2)/(v1-v2)**6.0
            phinpv[j]=A*64*(v-v1)**3.0*(v2-v)**3.0/(v2-v1)**6.0
   

    if bdytype=="stan" or bdytype=="max" or bdytype=="hor":
    
        sigv0=0.0
        if datatype==object:
            sigv0=Decimal(0)
        
        #dsignpv[0]=0.0
        
        dsignpv[:]=0.0
        drnpu[0]=-1/(4*drnpv[0])*(1-2*M0/rnpv[0]+(Q/rnpv[0])**2-Lambda*rnpv[0]**2/3)
        dphinpu[0]=0.0
        dphinpv[0]=0.0
        dsignpu[0]=0.0
        signpv[0]=0.0
        massnpv[0]=M0
        
        rcosm=rc(massnpv[0],Q,Lambda)
        for j in range(0,Nv*scal-1):
            
            #print(rcosm)
            if rnpv[j]+dt*drnpv[j]>0.0 and rnpv[j]+dt*drnpv[j]<rcosm-rcosmtol:
                ###Predictor###
                signpv[j+1]=sigv0
                rnpv[j+1]=rnpv[j]+dt*drnpv[j]
                drnpv[j+1]=drnpv[j]+dt*Coneq(drnpv[j],dsignpv[j],dphinpv[j],rnpv[j])
                drnpu[j+1]=drnpu[j]+dt*Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)
                dphinpu[j+1]=dphinpu[j]+dt*Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j]) 
                dsignpu[j+1]=dsignpu[j]+dt*Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q) 
                
                ###Corrector###           
                signpv[j+1]=sigv0
                rnpv[j+1]=rnpv[j]+1/2*dt*(drnpv[j]+drnpv[j+1])
                drnpv[j+1]=drnpv[j]+1/2*dt*(Coneq(drnpv[j],dsignpv[j],dphinpv[j],rnpv[j])+Coneq(drnpv[j+1],dsignpv[j+1],dphinpv[j+1],rnpv[j+1]))
                drnpu[j+1]=drnpu[j]+1/2*dt*(Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)+Rfunc(drnpv[j+1],drnpu[j+1],rnpv[j+1],signpv[j+1],Q,Lambda))
                dphinpu[j+1]=dphinpu[j]+1/2*dt*(Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j])+Phifunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1]))                     
                dsignpu[j+1]=dsignpu[j]+1/2*dt*(Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q)+Sigfunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1],signpv[j+1],Q))
                
                massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1])*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
                rcosm=rc(massnpv[j+1],Q,Lambda)
                #print(rcosm)
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
    dphinpu=dphinpu[::scal]  
    dphinpv=dphinpv[::scal]
    dsignpu=dsignpu[::scal]
    dsignpv=dsignpv[::scal]
    drnpu=drnpu[::scal]
    drnpv=drnpv[::scal]
    massnpv=massnpv[::scal]
    
    print(rnpv,Nv)
    
    
    
    #drnpu=None
    #drnpv=None
    #dsignpu=None
    #dsignpv=None
    #dphinpu=None
    #dphinpv=None
    return (rnpv, signpv, phinpv,drnpv,dsignpv,dphinpv,drnpu,dsignpu,dphinpu,massnpv) #dphinpu)#dsignpu, dphinpu)    



#########ODE functions###############
def Rfunc(drnpvf,drnpuf,rnpf,signpf,Q,Lambda):
    return -drnpvf*drnpuf/rnpf-np.exp(signpf)/(4.0*rnpf)*(1.0-np.power((Q/rnpf),2.0)-Lambda*np.power(rnpf,2.0))

def Sigfunc(drnpvf,drnpuf,dphinpuf,dphinpvf,rnpf,signpf,Q):
    return 2.0*drnpuf*drnpvf/np.power(rnpf,2.0)+np.exp(signpf)/(2.0*np.power(rnpf,2.0))*(1.0-2.0*np.power((Q/rnpf),2.0))-2.0*dphinpuf*dphinpvf

def Phifunc(drnpvf,drnpuf,dphinpuf,dphinpvf,rnpf):
    return -1/rnpf*(drnpuf*dphinpvf+drnpvf*dphinpuf)

def Coneq(drnpf,dsignpf,dphinpf,rnpf):
    return drnpf*dsignpf-rnpf*np.power(dphinpf,2.0)


def interp(i,array,urange,du):
    #array=array[0]
    
    #print(x,y)
    try:
        x=[urange[i-1],urange[i],urange[i+1],urange[i+2]]
        y=[array[i-1],array[i],array[i+1],array[i+2]]
        value=interpolate.interp1d(x,y,kind='cubic')(urange[i]+du/2)#+1/2*(urange[i+1]-urange[i]))
    except:
        print("Used Linear Interpolation")
        y=[array[i],array[i+1]]
        value=(y[0]+y[1])/2
    return value