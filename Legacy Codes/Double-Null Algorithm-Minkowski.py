# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:32:39 2018

@author: Brandon
"""

###Python Packages###
import numpy as np
import matplotlib.pyplot as plt
import math as mth
###Initial Values###

M0=0.0
Q=0.0
N=80
scal=1
umax=35
vmax=35
ru0=5.0

dr0v=.50
bdytype="max"
scalarfield=False
Elist=[1]

print
#####################

scalf=float(scal)
u0=0.0
v0=0.0
rv0=ru0
phiu0=0.0
phiv0=0.0


if M0==0.0:
    du00=1/N
    du0=du00
else:
    du00=M0/N
    du0=du00
dv00=du00
dv0=du0

Nu=int(umax/du0)
Nv=int(vmax/dv0)
print("Number of points for lowest iteration is "+str(Nu*Nv)+","+str(Nu)+"X"+str(Nv))
print("Number of points for highest iteration is "+str(Nu*Nv*max(Elist)**2)+","+str(Nu*max(Elist))+"X"+str(Nv*max(Elist)))

rplus=M0+(M0**2.0-Q**2.0)**.5
rminus=M0-(M0**2.0-Q**2.0)**.5
print("r+: "+str(rplus))
print("r-: "+str(rminus))
#####################
###Boundary Condition Function###
def boundary(scal,Es,bdytype):
    
    rnp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    phinp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    signp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    drnp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    dphinp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    dsignp=np.zeros((Nu*scal*Es,Nv*scal*Es))
    
    drnpu=np.zeros((Nu*scal*Es))
    drnpv=np.zeros((Nv*scal*Es))
    
    dsignpu=np.zeros((Nu*scal*Es))
    dsignpv=np.zeros((Nv*scal*Es))
    
    dphinpu=np.zeros((Nu*scal*Es))
    dphinpv=np.zeros((Nv*scal*Es))
    
    

    rnp[0][0]=ru0 
    drnpv[0]=dr0v
    
    dt=du0/(scalf*Es)
    
    if scalarfield==True:
        A=.115
        #v=1.0
        v1=0.0
        v2=1.0
        for i in range(0,Nu*scal*Es):
            #dsignpu[i]=0.0
            phinp[i][0]=0.0
            dphinpu[i]=0.0
        for i in range(0,Nv*scal*Es):
                #dsignpv[i]=0.0
            phinp[0][i]=0.0
            v=Nv*scal*Es
            dphinpv[i]=64*A*(v-v1)**2.0*(5*v-2*v1-3*v2)*(v-v2)/(v1-v2)**6.0
            phinp[0][i]=0.0
    else:
        for i in range(0,Nu*scal*Es):
            #dsignpu[i]=0.0
            dphinpu[i]=0.0
        for i in range(0,Nv*scal*Es):
                #dsignpv[i]=0.0
            dphinpv[i]=0.0

    if bdytype=="stan":
        sigu0= 0.0
        sigv0=0.0
        dsignpu[0]=0.0
        
        #drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v*ru0**2.0)*(Q**2.0-2*M0*ru0+ru0**2.0)
        dsignpv[0]=dsignpu[0]
        
        
        for i in range(0,Nv*scal*Es-1):
            dsignpv[i]=0.0
            rnp[0][i+1]=rnp[0][i]+dt*drnpv[i]
            drnpv[i+1]=drnpv[i]+dt*(drnpv[i]*dsignpv[i]-rnp[0][i]*dphinpv[i]**2.0)
            signp[0][i]=sigv0
        drnpu[0]=-mth.exp(sigu0)/(4.0*dr0v*rnp[0][-1]**2.0)*(Q**2.0-2*M0*rnp[0][-1]+rnp[0][-1]**2.0)
        for i in range(0,Nu*scal*Es-1):
            dsignpu[i]=0.0
            rnp[i+1][-1]=rnp[i][-1]+dt*drnpu[i]
            drnpu[i+1]=drnpu[i]+dt*(drnpu[i]*dsignpu[i]-rnp[i][0]*dphinpu[i]**2.0)
            signp[i][-1]=sigu0
            
    
    if bdytype=="edd":
        
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
    
    return (rnp, signp, phinp)
    


#####################
#####################
###Applying Boundary Conditions###
Emax=max([Elist])

rnpf=np.zeros((Nu*max(Elist),Nv*max(Elist),len(Elist)))
signpf=np.zeros((Nu*max(Elist),Nv*max(Elist),len(Elist)))
phinpf=np.zeros((Nu*max(Elist),Nv*max(Elist),len(Elist)))
temprnp=np.zeros((Nu*max(Elist),Nv*max(Elist)))
tempsignp=np.zeros((Nu*max(Elist),Nv*max(Elist)))
tempphinp=np.zeros((Nu*max(Elist),Nv*max(Elist)))
  


for k in range(0,len(Elist)):
    temprnp, tempsignp, tempphinp = boundary(scal,Elist[k],bdytype)
    
    for i in range(0,Nu*Elist[k]):
        for j in range(0,Nv*Elist[k]):
            rnpf[i][j][k]=temprnp[i][j]
            signpf[i][j][k]=tempsignp[i][j]
            phinpf[i][j][k]=tempphinp[i][j]
            





#####################
###Defining Propagation Algorithm###
def x4giver(un,vn,Es,k):
    Es=float(Es)
    du0=du00/(Es)
    dv0=du0
    #dv0
    
    r1=rnpf[un][vn][k]
    r2=rnpf[un][vn+1][k]
    r4=rnpf[un+1][vn+1][k]
    phi1=phinpf[un][vn][k]
    phi2=phinpf[un][vn+1][k]
    phi4=phinpf[un+1][vn+1][k]
    sig1=signpf[un][vn][k]
    sig2=signpf[un][vn+1][k]
    sig4=signpf[un+1][vn+1][k]
    
    
    r0=(r1+r4)/2.0
    
    sig0=(sig1+sig4)/2.0
    
    try:
        r3i=(-r2+r1+r4)+(r4-r2)*(r2-r1)/r0+du0*dv0*mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
        r3=(-r2+r1+r4)+(r3i-r1+r4-r2)*(r2-r1+r4-r3i)/(4*r0)+du0*dv0*mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
    
        dru0=(r3-r1+r4-r2)/(2.0*du0)
        drv0=(r2-r1+r4-r3)/(2.0*dv0)

        phi3i=(-phi2+phi1+phi4)+du0*dv0/r0*(dru0*(phi2-phi1)/dv0+drv0*(phi4-phi2)/du0)
        phi3=(-phi2+phi1+phi4)+du0*dv0/r0*(dru0*(phi2-phi1+phi4-phi3i)/(2*dv0)+drv0*(phi3i-phi1+phi4-phi2)/(2*du0))
    
        dphiu0=(phi3-phi1+phi4-phi2)/(2.0*du0)
        dphiv0=(phi2-phi1+phi4-phi3)/(2.0*dv0)

        #sig3=(-sig2+sig1+sig4)-du0*dv0*(2.0*dru0*drv0/(r0)**2.0+mth.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2*dphiu0*dphiv0)
        sig3=(sig4-sig2-sig1)-du0*dv0*(mth.exp(sig0)*(M0/r0**3.0-3*Q**2/(2*r0**4.0))-2*dphiu0*dphiv0)
    
        #m=(1.0+4.0*mth.exp(-sig0)*dru0*drv0)*r0/2.0+(Q**2.0)/(2.0*r0)
        m=M0
    
        dsigu0=(sig3-sig1+sig4-sig2)/(2.0*du0)
        dsigv0=(sig2-sig1+sig4-sig3)/(2.0*dv0)
    
        druv=-(r3-r1+r4-r2)*(r2-r1+r4-r3)/(4*r0*du0*dv0)-mth.exp(sig0)/(4*r0)*(1-(Q**2.0)/(r0**2.0))
        dphiuv=-1/r0*(dru0*(phi2-phi1+phi4-phi3)/(2*dv0)+drv0*(phi3-phi1+phi4-phi2)/(2*du0))
        dsiguv=2.0*dru0*drv0/(r0)**2.0+mth.exp(sig0)/(2.0*r0**2.0)*(1.0-2.0*(Q**2.0)/(r0**2.0))-2*dphiu0*dphiv0
    except OverflowError:
        
        r3i=np.nan
        r3=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi3i=np.nan
        phi3=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig3=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        dsigv0=np.nan
    
        druv=np.nan
        dphiuv=np.nan
        dsiguv=np.nan
        
    if r3<.01:
        r3i=np.nan
        r3=np.nan
    
        dru0=np.nan
        drv0=np.nan

        phi3i=np.nan
        phi3=np.nan
    
        dphiu0=np.nan
        dphiv0=np.nan

        sig3=np.nan
    
        m=np.nan
    
        dsigu0=np.nan
        dsigv0=np.nan
    
        druv=np.nan
        dphiuv=np.nan
        dsiguv=np.nan
        
    answer=[r3,phi3,sig3,m,drv0,dsigu0]
    
    return answer

drunp=np.zeros((Nu*max(Elist),Nv*max(Elist),len(Elist)))
dsigunp=np.zeros((Nu*max(Elist),Nv*max(Elist),len(Elist)))
drdvnplist=[]
drdvrlist=[]
dsigdvnplist=[]
dsigdvrlist=[]
rminuslist=[]
rpluslist=[]


###Applying Propagation Algorithm###
for k in range(0,len(Elist)):
    for j in range(0,Nv*Elist[k]-1):
        for i in range(0,Nu*Elist[k]-1):
        
            answer=x4giver(i,-j-2,Elist[k],k)
            rnpf[i+1][-j-2][k]=answer[0]
            phinpf[i+1][-j-2][k]=answer[1]
            signpf[i+1][-j-2][k]=answer[2]
            #massnp[i+1][j+1]=answer[3]
            drunp[i+1][-j-2][k]=answer[4]
            dsigunp[i+1][-j-2]=answer[5]
            #dsignp[i+1][j+1]=answer[6]


###Finding Horizons###           
for j in range(0,Nv*Elist[0]):
    for i in range(0,Nu*Elist[0]):
        if drunp[i][j][0]>0.0:
            drdvnplist.append(i/(Nu*Elist[0])*umax)
            drdvrlist.append(rnpf[i][j][0])
            rminuslist.append(rminus)
            break
        if i==Nu*Elist[0]-1:
            #drdvnplist.append(umax)
            #drdvrlist.append(rnpf[Nu*Elist[0]-1][j][0])
            drdvnplist.append(np.nan)
            drdvrlist.append(np.nan)
            rminuslist.append(rminus)
            #drdvnplist.append(nan)
            #drdvrlist.append(nan)
            break
        else:
            continue
            
for j in range(0,Nv*Elist[0]):
    for i in range(0,Nu*Elist[0]):
        if dsigunp[i][j][0]<0.0:
            dsigdvnplist.append(i/(Nu*Elist[0])*umax)
            dsigdvrlist.append(rnpf[i][j][0])
            rpluslist.append(rplus)
            break 
        if i==Nu*Elist[0]-1:
            #drdvnplist.append(umax)
            #drdvrlist.append(rnpf[Nu*Elist[0]-1][j][0])
            dsigdvnplist.append(np.nan)
            dsigdvrlist.append(np.nan)
            rpluslist.append(rplus)
            #drdvnplist.append(nan)
            #drdvrlist.append(nan)
            break
        else:
            continue       

        
urange=np.zeros((Nu*Elist[0]))
vrange=np.zeros((Nv*Elist[0]))
#vrange=np.empty([Nv*Elist[1]])


dt=du0/float(Elist[0])
for i in range(0,Nu*Elist[0]):
    urange[i]=i*dt
for i in range(0,Nv*Elist[0]):
    vrange[i]=i*dt
if len(Elist)>1:
    urange2=np.zeros((Nu*Elist[1]))
    dt=du0/float(Elist[1])
    for i in range(0,Nu*Elist[1]):
        urange2[i]=i*dt

    #for i in range(0,Nv*Elist[1]):
        #vrange[i]=i*dt




###Displaying Results###       
#print(rnpf)

rlist=[]
rlist2=[]
siglist=[]
siglist2=[]
philist=[]
diflist=[]
diflist2=[]
dudvlist=[]
uscale=[]
uvalue=u0

uloc=1/20

if len(Elist)>1:
    for i in range(0,Nu*Elist[0]):
        rlist.append(rnpf[i][int(uloc*Nv*Elist[0])][0])
        siglist.append(signpf[i][int(uloc*Nv*Elist[0])][0])
        philist.append(phinpf[i][int(uloc*Nv*Elist[0])][0])
    for i in range(0,Nu*Elist[1]):
        rlist2.append(rnpf[i][int(uloc*Nv*Elist[1])][1])
if len(Elist)==1:
    for i in range(0,Nu*Elist[0]):
        rlist.append(rnpf[i][int(uloc*Nv*Elist[0])][0])
        siglist.append(signpf[i][int(uloc*Nv*Elist[0])][0])
        philist.append(phinpf[i][int(uloc*Nv*Elist[0])][0]) 
       
for i in range(0,Nu*Elist[0]):    
    uvalue=du0/(1.0-2.0*M0/rlist[i]+Q**2.0/rlist[i]**2.0)+uvalue
    uscale.append(uvalue)
    
if len(Elist)==1:
    sample=plt.plot(urange, rlist)
    #sample=plt.plot(uscale, rlist)
    plt.xlabel('u')
    plt.ylabel('r(u,v)')
    plt.grid()
    plt.title('Radius vs u (v='+str(int(uloc*vmax))+')')
    #plt.savefig('ReisnerN r(u)-N=200.png',dpi=300)
plt.show()

if len(Elist)>1:
    
    #rlist2v2=rlist2[::int(Elist[1]/Elist[0])]
    #rlistdif=np.zeros((len(rlist)))
    
    #for i in range(0,len(rlist)):
        #rlistdif[i]=(rlist[i]-rlist2v2[i])/rlist2v2[i]
    #sample=plt.plot(urange,rlistdif)
    
    sample=plt.plot(urange, rlist)
    sample=plt.plot(urange2, rlist2)
    plt.xlabel('u')
    plt.ylabel('r(u,v)')
    plt.grid()
    plt.title('Radius vs u (v='+str(int(uloc*vmax))+')')
    #plt.savefig('1) Minkowski r(u,v)-N=200.png',dpi=300)
    plt.show()

sample2=plt.plot(urange,siglist)
#sample=plt.plot(uscale, siglist)
plt.xlabel('u')
plt.ylabel('sigma(u,v)')
plt.grid()
plt.title('Sigma vs u (v='+str(int(uloc*vmax))+')')
#plt.savefig('ReisnerN sigma(u,v)-N=200.png',dpi=300)
plt.show()

sample3=plt.plot(urange,philist)
plt.xlabel('u')
plt.ylabel('Scalar Field(u,v)')
plt.grid()
plt.title('Scalar Field vs u (v='+str(int(uloc*vmax))+')')
#plt.savefig('Scaling-Standard-to-Eddington.png',dpi=300)
plt.show()

sample3=plt.plot(vrange,drdvnplist)
plt.xlabel('v')
plt.ylabel('u')
plt.grid()
plt.title('Inner Event Horizon Position')
#plt.savefig('Event Horizon Position(u,v)-M=1-Q=.09-N=200.png',dpi=300)
plt.show()

sample3=plt.plot(vrange,drdvrlist)
sample3=plt.plot(vrange,rminuslist,'r--')
plt.xlabel('v')
plt.ylabel('r(v)')
plt.grid()
plt.title('Inner Event Horizon Position (r(u,v))')
#plt.savefig('Event Horizon Position r(u,v)-M=1-Q=.09-N=200.png',dpi=300)
plt.show()

#sample3=plt.plot(vrange,dsigdvnplist)
#plt.xlabel('v')
#plt.ylabel('u')
#plt.grid()
#plt.title('Outer Event Horizon Position')
#plt.savefig('Event Horizon Position(u,v)-M=1-Q=.09-N=200.png',dpi=300)
#plt.show()

#sample3=plt.plot(vrange,dsigdvrlist)
#sample3=plt.plot(vrange,rpluslist,'r--')
#plt.xlabel('v')
#plt.ylabel('r(v)')
#plt.grid()
#plt.title('Outer Event Horizon Position (r(u,v))')
#plt.savefig('Event Horizon Position r(u,v)-M=1-Q=.09-N=200.png',dpi=300)
#plt.show()

rnp=np.empty([Nu*Elist[0],Nv*Elist[0]])
signp=np.empty([Nu*Elist[0],Nv*Elist[0]])
phinp=np.empty([Nu*Elist[0],Nv*Elist[0]])

for i in range(0,Nu*Elist[0]): 
    for j in range(0,Nv*Elist[0]):
        rnp[i][j]=rnpf[i][j][0]
        signp[i][j]=signpf[i][j][0]
        phinp[i][j]=phinpf[i][j][0]

###Heatmaps for radius, sigma, scalar field###
levels=np.arange(10.0, np.nanmax(rnp), 5.0 )

if rminus==0.0:
    levels=np.insert(levels,0,[rplus,ru0])
else:
    levels=np.insert(levels,0,[rminus,rplus,ru0])

plot1=plt.imshow(rnp,cmap=plt.cm.cool,extent=(0,vmax,0,umax),aspect='auto',origin='lower')
plot4=plt.contour(rnp, levels,linewidths=0.5,colors='black',extent=(0,vmax,0,umax),aspect='auto',origin='lower')
plt.clabel(plot4, levels,  # label every second level
           inline=0, inline_spacing=0, 
           fmt='%1.0f',rightside_up=True,
           fontsize=8)
sample3=plt.plot(vrange,drdvnplist)
#sample3=plt.plot(vrange,dsigdvnplist)
plt.xlabel('v')
plt.ylabel('u')
plt.xlim()
plt.colorbar(plot1)
plt.title('Radius vs (u,v) Coordinates')
plt.savefig('Test.png',dpi=300)
plt.show()

plot3=plt.imshow(signp,cmap=plt.cm.afmhot,extent=(0,vmax,0,umax),aspect='auto',origin='lower')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar(plot3)
plt.title('Sigma vs (u,v) Coordinates')
plt.show()

plot2=plt.imshow(phinp,cmap=plt.cm.afmhot,extent=(0,vmax,0,umax),aspect='auto',origin='lower')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar(plot2)
plt.title('Scalar Field vs (u,v) Coordinates')
plt.show()



