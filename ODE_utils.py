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
#import scipy
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

def boundaryv(scal,ubdytype,bdytype,Nv,ru0,dr0v,dv0,vmax,M0,Q,Lambda,scalarfield,A,rcosmtol,scalarfieldtype,power):
    datatype=float
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
        
        #scalarfieldtype='pulse'#'power'
        if scalarfieldtype=='pulse':
            v1=1.0
            v2=5.0
            v1n=int(v1*(Nv*scal)/vmax)
            v2n=int(v2*(Nv*scal)/vmax)
        
       
            for j in range(v1n,v2n):
                    #dsignpv[i]=0.0
                v=j/(Nv*scal)*vmax
                dphinpv[j]=192*A*(v-v1)**2.0*(v-v2)**2.0*(-2*v+v1+v2)/(v1-v2)**6.0
                phinpv[j]=A*64*(v-v1)**3.0*(v2-v)**3.0/(v2-v1)**6.0
        
        if scalarfieldtype=='power':
            v1=1
            v2=vmax
            #power=-1
            v1n=int(v1*(Nv*scal)/vmax)
            v2n=int(v2*(Nv*scal)/vmax)


            for j in range(v1n,v2n):
                    #dsignpv[i]=0.0
                v=j/(Nv*scal)*vmax
                #try:
                dphinpv[j]=A*((v)**power*np.cosh(v-v1)**(-2.0)+power*(v)**(power-1)*np.tanh(v-v1))
                phinpv[j]=A*np.tanh(v-v1)*v**power
                #except:
                    #phinpv[j]=0.0
                    #phinpv[j]=0.0

    if bdytype=="stan" or bdytype=="max" or bdytype=="hor":
    
        sigv0=0.0
        if datatype==object:
            sigv0=Decimal(0)
        
        #dsignpv[0]=0.0
        
        dsignpv[:]=0.0
        drnpu[0]=-1/(4*drnpv[0])*(1-2*M0/rnpv[0]+(Q/rnpv[0])**2-Lambda*rnpv[0]**2/3)
        dphinpu[0]=0.0
        dphinpv[0]=0.0
        signpv[0]=0.0
        if ubdytype=="stan":
            dsignpu[0]=0.0
        elif ubdytype=="adapt":
            dsignpu[0]=0.0#-(0-np.log(dr0v)+3/2*np.log(2))
        
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
        
        signpv[0]=np.log(1.0-2.0*M0/ru0+Q**2.0/ru0**2.0-Lambda*ru0**2/3)
        
        dsignpv[0]=2*(3*Q**2-3*M0*ru0+ru0**4*Lambda)/(ru0*(-3*Q**2+ru0*(6*M0-3*ru0+ru0**3*Lambda)))*dr0v
        drnpv[0]=dr0v
        print('dr/dv is '+str(drnpv[0]))
        #sigv0=0.0
        if datatype==object:
            sigv0=Decimal(0)
        
        #dsignpv[0]=0.0
        
        #dsignpv[:]=0.0
        drnpu[0]=-1/(4*drnpv[0])*np.exp(signpv[0])*(1-2*M0/rnpv[0]+(Q/rnpv[0])**2-Lambda*rnpv[0]**2/3)
        dphinpu[0]=0.0
        #dphinpv[0]=0.0
        if ubdytype=="stan":
            dsignpu[0]=0.0
        elif ubdytype=="adapt":
            dsignpu[0]=-(-1/2*np.log(2)-np.log(dr0v)+3/2*np.log(2))
        elif ubdytype=="edd":
            dsignpu[0]=2*(3*Q**2-3*M0*ru0+Lambda*ru0**4)/(ru0*(-3*Q**2+ru0*(6*M0-3*ru0+Lambda*ru0**3)))*drnpu[0]
            print("Using Full Eddington Coordinates")
        #signpv[0]=0.0
        massnpv[0]=M0
        
        rcosm=rc(massnpv[0],Q,Lambda)
        for j in range(0,Nv*scal-1):
            
            #print(rcosm)
            if rnpv[j]+dt*drnpv[j]>0.0 and rnpv[j]+dt*drnpv[j]<rcosm-rcosmtol:
                ###Predictor###
                
                rnpv[j+1]=rnpv[j]+dt*drnpv[j]
                drnpv[j+1]=drnpv[j]+dt*Coneq(drnpv[j],dsignpv[j],dphinpv[j],rnpv[j])
                drnpu[j+1]=drnpu[j]+dt*Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)
                
                dphinpu[j+1]=dphinpu[j]+dt*Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j]) 
                dsignpu[j+1]=dsignpu[j]+dt*Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q) 
                ###
                signpv[j+1]=signpv[j]+dt*dsignpv[j]
                
                
                
                massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
                
               
                
                
                dsignpv[j+1]=2*(3*Q**2-3*massnpv[j+1]*rnpv[j+1]+rnpv[j+1]**4*Lambda)/(rnpv[j+1]*(-3*Q**2+rnpv[j+1]*(6*massnpv[j+1]-3*rnpv[j+1]+rnpv[j+1]**3*Lambda)))*drnpv[j+1]
                #signpv[j+1]=np.log(np.abs(1.0-2.0*massnpv[j+1]/rnpv[j+1]+Q**2.0/rnpv[j+1]**2.0-Lambda*rnpv[j+1]**2/3))
                
                #massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
                
                
                #dsignpv[j+1]=2*((massnpv[j+1]*rnpv[j+1]-Q**2-1/3*Lambda*rnpv[j+1]**4)/(rnpv[j+1]**3-2*massnpv[j+1]*rnpv[j+1]**2+Q**2*rnpv[j+1]-Lambda/3*rnpv[j+1]**5))*drnpv[j+1]
                
                
                ###Corrector###           
                #signpv[j+1]=sigv0
                rnpv[j+1]=rnpv[j]+1/2*dt*(drnpv[j]+drnpv[j+1])
                drnpv[j+1]=drnpv[j]+1/2*dt*(Coneq(drnpv[j],dsignpv[j],dphinpv[j],rnpv[j])+Coneq(drnpv[j+1],dsignpv[j+1],dphinpv[j+1],rnpv[j+1]))
                drnpu[j+1]=drnpu[j]+1/2*dt*(Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)+Rfunc(drnpv[j+1],drnpu[j+1],rnpv[j+1],signpv[j+1],Q,Lambda))
                dphinpu[j+1]=dphinpu[j]+1/2*dt*(Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j])+Phifunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1]))                     
                dsignpu[j+1]=dsignpu[j]+1/2*dt*(Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q)+Sigfunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1],signpv[j+1],Q))
                ###
                signpv[j+1]=signpv[j]+1/2*dt*(dsignpv[j]+dsignpv[j+1])
                
                
                
                massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
                
                
                
                
                dsignpv[j+1]=2*(3*Q**2-3*massnpv[j+1]*rnpv[j+1]+rnpv[j+1]**4*Lambda)/(rnpv[j+1]*(-3*Q**2+rnpv[j+1]*(6*massnpv[j+1]-3*rnpv[j+1]+rnpv[j+1]**3*Lambda)))*drnpv[j+1]
                #signpv[j+1]=np.log(1.0-2.0*massnpv[j+1]/rnpv[j+1]+Q**2.0/rnpv[j+1]**2.0-Lambda*rnpv[j+1]**2/3)
                
                #massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
               
                
                #dsignpv[j+1]=2*((massnpv[j+1]*rnpv[j+1]-Q**2-1/3*Lambda*rnpv[j+1]**4)/(rnpv[j+1]**3-2*massnpv[j+1]*rnpv[j+1]**2+Q**2*rnpv[j+1]-Lambda/3*rnpv[j+1]**5))*drnpv[j+1]
                
                rcosm=rc(massnpv[j+1],Q,Lambda)
                #print(rcosm)
            else:
                break
  
    ############################
    ############################
                
    elif bdytype=="EH":
        beta=.1
        #power=3
        #power=-power
        print(power)
        rplus=M0+(M0**2.0-Q**2.0)**(.5)
        v0=((rplus-ru0)/beta)**(-1/power)
        
        signpv=np.zeros((Nv*scal),dtype=datatype)
        dsignpv=np.zeros((Nv*scal),dtype=datatype)
        
        
        
        
        drnpv[0]=beta*power*(v0)**(-power-1) 
        #dsignpv[0]=-(power+1)/v0#kappaq(ru0,Q)*0.0
        rnpv[0]=ru0
        phinpv[0]=0.0
        
        dphinpv[0]=(beta*power*(power+1)*(v0)**(-power-2)/(rplus-beta*(v0)**(-power)))**(.5)
        print('dr/dv is '+str(drnpv[0]))
        #sigv0=0.0
        
        dsignpv[0]=kappaq(ru0,Q)#np.abs(Rfunc(drnpv[0],drnpu[0],rnpv[0],signpv[0],Q,Lambda)/drnpu[0])#-(power+1)/v0#kappaq(ru0,Q)*0.0
        signpv[0]=kappaq(ru0,Q)*v0#np.log(.5*np.abs(drnpu[0]))
        
        M=Mq(.90,Q)
        print(0.0,0.0,M)
        #signpv[0]=.1#np.log(.5*np.abs(drnpv[0]))
        drnpu[0]=-1/(4*drnpv[0])*np.exp(signpv[0])*(1-2*M/rnpv[0]+(Q/rnpv[0])**2-Lambda*rnpv[0]**2/3)
        
        
        
        print(0.0,0.0,drnpu[0],drnpv[0])
        
        dphinpu[0]=0.0
        #dphinpv[0]=0.0
        if ubdytype=="stan":
            dsignpu[0]=0.0
        elif ubdytype=="adapt":
            dsignpu[0]=-(-1/2*np.log(2)-np.log(dr0v)+3/2*np.log(2))
        elif ubdytype=="edd":
            dsignpu[0]=2*(3*Q**2-3*M0*ru0+Lambda*ru0**4)/(ru0*(-3*Q**2+ru0*(6*M0-3*ru0+Lambda*ru0**3)))*drnpu[0]
            print("Using Full Eddington Coordinates")
        #signpv[0]=0.0
        massnpv[0]=M0
        v=v0
        
        rcosm=rc(massnpv[0],Q,Lambda)
        for j in range(0,Nv*scal-1):
            
            #print(rcosm)
            if rnpv[j]+dt*drnpv[j]>0.0 and rnpv[j]+dt*drnpv[j]<rcosm-rcosmtol:
                ###Predictor###
                v=v+dt
                
                rnpv[j+1]=rnpv[j]+dt*drnpv[j]
                drnpv[j+1]=beta*power*(v)**(-power-1)
                drnpu[j+1]=drnpu[j]+dt*Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)
                
                dphinpv[j+1]=((beta*power*(power+1)*(v)**(-power-2)+dsignpv[j+1]*beta*power*(v)**(-power-1))/(rplus-beta*(v)**(-power)))**(.5)
                phinpv[j+1]=phinpv[j]+dt*dphinpv[j]
                
                dsignpv[j+1]=kappaq(rnpv[j+1],Q)
                #signpv[j+1]=kappaq(ru0,Q)*v0
                #dsignpv[j+1]=-(power+1)/v#kappaq(rnpv[j+1],Q)*0.0
                signpv[j+1]=signpv[j]+dt*dsignpv[j]
                
                dphinpu[j+1]=dphinpu[j]+dt*Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j]) 
                dsignpu[j+1]=dsignpu[j]+dt*Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q) 
                ###
                
                #dsignpv[j+1]=0.0
                #signpv[j+1]=0.0
       
                #massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0

                
                ###Corrector###           
                #signpv[j+1]=sigv0
                rnpv[j+1]=rnpv[j]+0.5*dt*(drnpv[j]+drnpv[j+1])
                #drnpv[j+1]=beta*power*(v)**(-power-1)
                drnpu[j+1]=drnpu[j]+0.5*dt*(Rfunc(drnpv[j],drnpu[j],rnpv[j],signpv[j],Q,Lambda)+Rfunc(drnpv[j+1],drnpu[j+1],rnpv[j+1],signpv[j+1],Q,Lambda))
                
                dsignpv[j+1]=kappaq(rnpv[j+1],Q)
                #signpv[j+1]=kappaq(ru0,Q)*v0
                #dsignpv[j+1]=-(power+1)/v#kappaq(rnpv[j+1],Q)*0.0
                signpv[j+1]=signpv[j]+0.5*dt*(dsignpv[j]+dsignpv[j+1])
                
                phinpv[j+1]=phinpv[j]+0.5*dt*(dphinpv[j]+dphinpv[j+1])
                
                dphinpu[j+1]=dphinpu[j]+0.5*dt*(Phifunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j])+Phifunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1])) 
                dsignpu[j+1]=dsignpu[j]+0.5*dt*(Sigfunc(drnpv[j],drnpu[j],dphinpu[j],dphinpv[j],rnpv[j],signpv[j],Q)+Sigfunc(drnpv[j+1],drnpu[j+1],dphinpu[j+1],dphinpv[j+1],rnpv[j+1],signpv[j+1],Q)) 
                massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
                
                
                #signpv[j+1]=np.log(1.0-2.0*massnpv[j+1]/rnpv[j+1]+Q**2.0/rnpv[j+1]**2.0-Lambda*rnpv[j+1]**2/3)
                
                #massnpv[j+1]=(1+4.0*drnpu[j+1]*drnpv[j+1]*np.exp(-signpv[j+1]))*rnpv[j+1]/2.0+Q**2.0/(2*rnpv[j+1])-Lambda*rnpv[j+1]**3.0/6.0
               
                
                #dsignpv[j+1]=2*((massnpv[j+1]*rnpv[j+1]-Q**2-1/3*Lambda*rnpv[j+1]**4)/(rnpv[j+1]**3-2*massnpv[j+1]*rnpv[j+1]**2+Q**2*rnpv[j+1]-Lambda/3*rnpv[j+1]**5))*drnpv[j+1]
                
                rcosm=rc(massnpv[j+1],Q,Lambda)
                #print(rcosm)
            else:
                break   
        
        
    #print(Nv)      
    
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
    
    #print(rnpv,Nv)
    
    
    
    #drnpu=None
    #drnpv=None
    #dsignpu=None
    #dsignpv=None
    #dphinpu=None
    #dphinpv=None
    return (rnpv, signpv, phinpv,drnpv,dsignpv,dphinpv,drnpu,dsignpu,dphinpu,massnpv) #dphinpu)#dsignpu, dphinpu)    



#########ODE functions###############
def trapint(array):
    return 0.5*(np.cumsum(np.roll(array,1)+array)-(np.roll(array,1)+array)[0])


def Rfunc(drnpvf,drnpuf,rnpf,signpf,Q,Lambda):
    return -drnpvf*drnpuf/rnpf-np.exp(signpf)/(4.0*rnpf)*(1.0-np.power((Q/rnpf),2.0)-Lambda*np.power(rnpf,2.0))

def Sigfunc(drnpvf,drnpuf,dphinpuf,dphinpvf,rnpf,signpf,Q):
    return 2.0*drnpuf*drnpvf/np.power(rnpf,2.0)+np.exp(signpf)/(2.0*np.power(rnpf,2.0))*(1.0-2.0*np.power((Q/rnpf),2.0))-2.0*dphinpuf*dphinpvf

def Phifunc(drnpvf,drnpuf,dphinpuf,dphinpvf,rnpf):
    return -1/rnpf*(drnpuf*dphinpvf+drnpvf*dphinpuf)

def Coneq(drnpf,dsignpf,dphinpf,rnpf):
    return drnpf*dsignpf-rnpf*np.power(dphinpf,2.0)

def Kretsch(r,drdv,drdu,dphidv,dphidu,m,Q,Lambda):
    return 16/r**6.0*((m-3*Q**2/(2*r)+Lambda/6*r**3)+r/2*(1-2*m/r+Q**2.0/r**2.0-Lambda/3*r**2.0)*(r*dphidu/drdu)*(r*dphidv/drdv))**2.0+16/r**6*(m-Q**2.0/(2*r)+Lambda*r**3/6)+16/r**6*(m-Q**2/r-r**3/3*Lambda)**2+4/r**4*(1-2*m/r+Q**2/r**2-Lambda/3*r**2)**2*(r*dphidu/drdu)**2*(r*dphidv/drdv)**2


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

#######
def drdvparam(x,y):
    
    def func(x,p,a):
        return p*np.log(x)+a     #A*x**p+b
    
    result,result2=optimize.curve_fit(func,x,y)
    
    return result

def massparam(x,y):
      
    def func(x,a,b,c):
        return a*x+b*np.log(x)+c#A*np.exp(a*x)*x**b
    
    result,result2=optimize.curve_fit(func,x,y)
    
    return result

def kappaq(rplus,Q):
    M=(Q**2+rplus**2)/(2*rplus)
    return (M**2-Q**2+M*(M**2-Q**2)**(.5))/(M+(M**2-Q**2)**(.5))**3

def Mq(rplus,Q):
    return (Q**2+rplus**2)/(2*rplus)
