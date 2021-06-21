# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 09:00:06 2019

@author: santiago valencia
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 20:32:08 2019

@author: santiago valencia
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

n=100
th=np.linspace(0,0.88*np.pi,n) #theta


P0=101325
T0=300

air1 = ct.Solution('gri30.xml')
YO2=0.16
VN2=(1-YO2)/YO2
YO2t=str(YO2)
YN2t=str(1-YO2)
comt='O2:'+YO2t+', N2:'+YN2t
air1.X = comt
pmol=air1.mean_molecular_weight

gasc = ct.Solution('gri30.xml')
ex=0.2
ngt=(17+12.5*((1+ex)*VN2+ex))
YgCO2=8/ngt
YgH2O=9/ngt
YgN2=(12.5*(1+ex)*VN2)/ngt
YgO2=(12.5*ex)/ngt
YgCO2t=str(YgCO2)
YgH2Ot=str(YgH2O)
YgN2t=str(YgN2)
YgO2t=str(YgO2)
comgt='CO2:'+YgCO2t+', H2O:'+YgH2Ot+', N2:'+YgN2t+', O2:'+YgO2t
gasc.X = comgt
pmolg=gasc.mean_molecular_weight



om=5100 #velocidad angular [rpm]
r=9 #relación de compresión
C1=130 #constante
C2=1.4 #constante
Ru=8.314472 #constante de gas kj/mol*K
R=Ru/pmol*1000 #constante del aire j/kg*K
Rg=Ru/pmolg*1000 #constante de la mezcla de gases
Nc=4 #Número de cilindros
Vm=1750 #Volumen desplazado [cc]
SB=1 #relación carera-diámetro
efin=0.8 #eficiencia de aletas
h0=150 #coeficiente de convección forzado
k=50 #coeficiente de conduccion por el cilindro 
Vd=Vm/Nc #volumen de cada cilindro total 
Vc=Vd/(r-1) #volumen clearence
B=(4.0/np.pi*Vd/SB)**(1./3.)/100 #Bore [m]
S=SB*B #Stroke [m]
re=10.4/2/100 #radio externo cilindro [m]
ri=B/2 #radio interno cilindro [m]
Vp=2*S*om*2*np.pi/60 #velocidad del pistón [m/s]
Ar=2 #relación de aeras con aletas
Ab=2*np.pi*(re)*S #area de la base [m2]
Ai=(2*np.pi*(ri)*S)/2 #área interna media [m2]
A0=Ab*Ar #area con aletas [m2]
l=14/100
Rc=l/S
Vmax=(Vd+Vc)/1e6
pci=44300000
efic=0.8
qin=pci*efic
Rad=B/2 #radio del cilindro
AFe=3.38+2.96*VN2
AF=(1+ex)*AFe
fvc=0.04
n1=1
efim=0.86

air1.TP=T0,P0
v0=air1.v
u0=air1.u

mm=Vd/v0/1e6
mF=1/(1+AF)*mm*(1-fvc)
mA=AF/(1+AF)*mm*(1-fvc)

gasc.TP=T0, P0
v0g=gasc.v
mg=Vc/1e6/v0g

mt=mm+mg
xi=mg/mt
Den=mm/Vmax
Deng=mg/Vmax
pmolt=1/(xi/pmolg+(1-xi)/pmol)



#Compresión

def Tsol(T,the):
    global Pi, Den, Deng, cp, cpg, cv, V1, va1, vg1

    air1.TD=T,Den
    Pa=air1.P
    cv=air1.cv
    cp=air1.cp
    va1=air1.v
    
    gasc.TD=T,Deng
    Pg=gasc.P
    cvg=gasc.cv
    cpg=gasc.cp
    vg1=gasc.v
    
    
    V1=Vmax*((1/r)+(1/2)*(1-1/r)*(1+np.cos(the)+Rc-np.sqrt(Rc**2-(np.sin(the))**2)))
    dcv=cvg-cv
    A=xi*dcv+cv
    Den=mm/V1
    Deng=mg/V1
    Dent=xi*Deng+(1-xi)*Den
    Pi=Den*pmolt/Dent/pmol*Pa+Deng*pmolt/Dent/pmolg*Pg 
   
    
    hg=C1*(((1+1/r)/2*V1)**(-0.06))*((T)**(-0.4))*(((Pi)/1e5)**0.8)*((Vp+C2)**0.8)
    UA=1/((1/(hg*Ai))+(np.log(re/ri)/(2*np.pi*k*S/100))+(1/(A0*efin*h0)))
    dvdth=Vmax/2*(1-1/r)*((np.sin(the)*np.cos(the))/(np.sqrt(Rc**2-(np.sin(the)**2)))-np.sin(the))
    Qdth=-(UA/(om*2*np.pi/60))*(T-T0)
    dTdth=(Qdth/mt/A)-(xi*Rg+(1-xi)*R)*T/A/V1*dvdth
    return dTdth

Pi=P0
Sol = odeint(Tsol,T0,th) 

Tfinal1=Sol[n-1,0]
Pfinal1=Pi


Vb=0.01*V1
mb=mt*Vb/V1
Tb=1/cpg*(cv*Tfinal1+(1-xi)*pci*mF/mm-Pi*Vb/mb) 
pcip=pci*mF/mm


th2=np.linspace(0.88*np.pi,1.2*np.pi,300)
tha=0.88*np.pi
thd=2*np.pi
thdd=2*np.pi
rad=0.0
ths=tha
#dT[:]=0.0
cf=0

#Combustión

def Tsol2(the,T):
    global mg, mm, Den, Deng, xi, Pi, tha, thd, rad,ths, cv, cp, cvg, cpg, UA, dT,thdd,Tcf,Pcf,Dcf,cf,vg2,mgg

    
    
    dth=the-tha
    tha=the
    
    if Den>0.0:
        air1.TD=T[0],Den
        Pa=air1.P
        cv=air1.cv
        cp=air1.cp
    else:
        Pa=0.0
        
    gasc.TD=T[2],Deng
    Pg=gasc.P
    cvg=gasc.cv
    cpg=gasc.cp
    vg2=gasc.v

    
    dTdth=np.zeros(3,float)
    v=Vmax*((1/r)+(1/2)*(1-1/r)*(1+np.cos(the)+Rc-np.sqrt(Rc**2-(np.sin(the))**2)))
    dvdth=Vmax/2*(1-1/r)*((np.sin(the)*np.cos(the))/(np.sqrt(Rc**2-(np.sin(the)**2)))-np.sin(the))
    
    Den=mm/v
    Deng=mg/v

    
    
    xi=mg/mt
    
    Dent=xi*Deng+(1-xi)*Den
    Pi=Den*pmolt/Dent/pmol*Pa+Deng*pmolt/Dent/pmolg*Pg
    
    if the<thd:
        Vf=30
        rad=Vf/(om*2*np.pi/60)*(the-(0.88*np.pi))
    else:
        Vf=0.
        rad=Rad
    if  rad>Rad-1e-8:
        thd=the
    
    hg=C1*(((1+1/r)/2*v)**(-0.06))*(((T[0]+T[2])/2)**(-0.4))*(((Pi)/1e5)**0.8)*((Vp+C2)**0.8)
    UA=1/((1/(hg*Ai))+(np.log(re/ri)/(2*np.pi*k*S/100))+(1/(A0*efin*h0)))
    Qdthu=-(UA)*(T[0]-T0)
    Qdthb=-(UA)*(T[2]-T0)
    Qu=(1-xi)*Qdthu
    Qb=xi*Qdthb
    
    dTdth[1]=(1/Rad**2)*(rad**2)*dvdth+(2*v*rad/Rad**2)*Vf/(om*2*np.pi/60)
    if mm>0.0:
        dmmdth=-mt*(1/v*dTdth[1]-T[1]/v**2*dvdth)
    else:
        dmmdth=0.0
    Mpb=-dmmdth*om*2*np.pi/60
    mg=mg-dmmdth*dth
    mm=mm+dmmdth*dth
    dVudth=(dmmdth+(v-T[1])/v**2*dvdth*mt)/(mt/v)
    if mm>1e-16:
        dTdth[0]=(Qu-Pi*dVudth*(om*2*np.pi/60)+Mpb*(cv-cp)*T[0])/(mm*cv*(om*2*np.pi/60))
    else:
        dTdth[0]=0
        mm=0.0
        Den=0.0
        
        if cf==0:
            thdd=the
            Tcf=T[2]
            Pcf=Pi
            Dcf=Deng
            mgg=mg
            cf=1
            
        
        
    dTdth[2]=(Qb-Pi*dTdth[1]*(om*2*np.pi/60)+Mpb*(cp*T[0]-cvg*T[2]+pcip))/(mg*cvg*(om*2*np.pi/60))
    #dT=dTdth
    
    
    return dTdth


Sol2=solve_ivp(Tsol2,(0.88*np.pi,1.16*np.pi),[Tfinal1,0,Tb],method='RK45')

Tfinal2=Tcf
Pfinal2=Pcf
mg=mgg

th3=np.linspace(thdd,2*np.pi-np.pi/9,n)

#Expansión

def Tsol3(T,the):
    global Pi, Deng, cpg, cvg, V, mg, vg3
    
    gasc.TD=T,Deng
    Pg=gasc.P
    cvg=gasc.cv
    cpg=gasc.cp
    vg3=gasc.v
   
    
    
    
    V=Vmax*((1/r)+(1/2)*(1-1/r)*(1+np.cos(the)+Rc-np.sqrt(Rc**2-(np.sin(the))**2)))
    Deng=mg/V
    Pi=Pg 
    hg=C1*(((1+1/r)/2*V1)**(-0.06))*((T)**(-0.4))*(((Pi)/1e5)**0.8)*((Vp+C2)**0.8)
    UA=1/((1/(hg*Ai))+(np.log(re/ri)/(2*np.pi*k*S/100))+(1/(A0*efin*h0)))
    dVdth=Vmax/2*(1-1/r)*((np.sin(the)*np.cos(the))/(np.sqrt(Rc**2-(np.sin(the)**2)))-np.sin(the))
    Qdth=-(UA/(om*2*np.pi/60))*(T-T0)
    
    dTdth=Qdth/mg/cvg-Rg*T/cvg/V*dVdth
    
    return dTdth

Sol3=odeint(Tsol3,Tcf,th3)


Tfinal3=Sol3[n-1,0]
Pfinal3=Pi
gasc.TP=Tfinal3,Pfinal3
u3=gasc.u
n=750
th4=np.linspace(2*np.pi-np.pi/9,3*np.pi,n)
tha=2*np.pi-np.pi/9
mg=mt
#Rechazo de gases
def Tsol4(T,the):
    global Pi, Deng, cpg, cvg, V, tha, u4,vg4,mg,mcr
    
    dth=the-tha
    tha=the
    
    
    gasc.TD=T,Deng
    Pg=gasc.P
    cvg=gasc.cv
    cpg=gasc.cp
    u4=gasc.u
    vg4=gasc.v
    
   
    
    V=Vmax*((1/r)+(1/2)*(1-1/r)*(1+np.cos(the)+Rc-np.sqrt(Rc**2-(np.sin(the))**2)))
    Deng=mg/V
    Pi=Pg
    hg=C1*(((1+1/r)/2*V1)**(-0.06))*((T)**(-0.4))*(((Pi)/1e5)**0.8)*((Vp+C2)**0.8)
    UA=1/((1/(hg*Ai))+(np.log(re/ri)/(2*np.pi*k*S/100))+(1/(A0*efin*h0)))
    dVdth=Vmax/2*(1-1/r)*((np.sin(the)*np.cos(the))/(np.sqrt(Rc**2-(np.sin(the)**2)))-np.sin(the))
    Qdth=-(UA/(om*2*np.pi/60))*(T[0]-T0)
    
    Dv=0.023
    Cf=0.4
    Av=np.pi/4*Dv**2
    K=cpg/cvg
    Co=np.sqrt(K*Rg*T[0])
    
    
    mcr=Deng*Cf*Av*Co*(2/(K+1))**((K+1)/2/(K+1))/(om*2*np.pi/60)
    mg=mg-mcr*dth
    
    dTdth=Qdth/mg/cvg-Rg*T[0]/V/cvg*dVdth-mcr/mg*Rg*T[0]/cvg
    
    return dTdth

Sol4=odeint(Tsol4,Tfinal3,th4)

Tfinal4=Sol4[n-1,0]
Pfinal4=Pi
th5=np.linspace(3*np.pi,4*np.pi,n)

#Admisión

tha=3*np.pi
def Tsol5(T,the):
    global Pi, Deng, Den, cv, cp, cpg, cvg, V,tha,mm,mg
    
    dth=the-tha
    tha=the
    

    air1.TP=(T+T0)/2,P0
    Pa=air1.P
    cv=air1.cv
    cp=air1.cp
    va5=air1.v
    
    gasc.TD=T,Deng
    Pg=gasc.P
    cvg=gasc.cv
    cpg=gasc.cp
    
    
    V=Vmax*((1/r)+(1/2)*(1-1/r)*(1+np.cos(the)+Rc-np.sqrt(Rc**2-(np.sin(the))**2)))
    Den=mm/V
    Deng=mg/V
    xi=mg/mt
    Dent=xi*Deng+(1-xi)*Den
    Pi=Den*pmolt/Dent/pmol*Pa+Deng*pmolt/Dent/pmolg*Pg
    
    hg=C1*(((1+1/r)/2*V1)**(-0.06))*((T)**(-0.4))*(((Pi)/1e5)**0.8)*((Vp+C2)**0.8)
    UA=1/((1/(hg*Ai))+(np.log(re/ri)/(2*np.pi*k*S/100))+(1/(A0*efin*h0)))
    dVdth=Vmax/2*(1-1/r)*((np.sin(the)*np.cos(the))/(np.sqrt(Rc**2-(np.sin(the)**2)))-np.sin(the))
    Qdth=-(UA/(om*2*np.pi/60))*(T-T0)
    
    Dv=0.035
    Cf=0.4
    Av=np.pi/4*Dv**2
    K=cp/cv
    Co=np.sqrt(K*R*T)
    
    mcr=Den*Cf*Av*Co*(2/(K+1))**((K+1)/2/(K+1))/(om*2*np.pi/60)
    mm=mm+mcr*dth
    
    mme=xi*mg+(1-xi)*mm
    dTdth=Qdth/mme/cv-R*T/V/cv*dVdth+mcr/mme/cv*(cp*T0-cv*T)
    
    return dTdth

Sol5=odeint(Tsol5,Tfinal4,th5)

Tfinal5=Sol5[:,0]
Pfinal5=Pi

#Flujos de energía
#qout=u3-u0
w01=R*((T0-Tfinal1)/(v0**n1-((va1+vg1)/2)**n1))*((((va1+vg1)/2)**n1-v0**n1)/n1-v0**n1*np.log(((va1+vg1)/2)/v0))+R*T0*np.log(((va1+vg1)/2)/v0)
w23=Rg*((Tfinal2-Tfinal3)/(vg2**n1-vg3**n1))*((vg3**n1-vg2**n1)/n1-vg2**n1*np.log(vg3/vg2))+Rg*Tfinal3*np.log(vg3/vg2)
wex=Pfinal4*(vg3-vg4)
wad=Pfinal5*(v0-vg4)
wnet=w23+w01+wad+wex
#Wnet=Qin-mm[0]*qout
efi=wnet/qin/mF*mt*100
Wnet=wnet*mt

#indicadores

imep=Wnet/((Vd+Vc)-Vc)*1000 #indicate mean efective pressure  [kPa]
Pi=(((Wnet*(om/60))/2)*Nc) #indicate power [w]
Up=2*S*(om/60) # mean piston speed [m/s]
Wb=efim*Wnet #break work [J]
Pb=efim*Pi #break power [w]
To=Pb/(2*np.pi*(om/60)) #torque [Nm]
Wf=Pi-Pb #friction power [W]
bmep=efim*imep #brake mean effective pressure [kPa]
BSP=Pb/((np.pi/4)*(B)**2)*Nc/1000 #break specific power [kW/m**2]
OPD=Pb/Vm         #output per displacement [kW/lt]
bsfc=mF*(om/60)/2*Nc/Pb*1000*3600*1000 #brake specific fuel consumption [gr/kW-h]
efiv=mA*v0/Vd*1e6*100 #Volumetric efficiency [%]

plt.plot(th,Sol[:,0],'b',label='Temperatura (K)')
plt.plot(Sol2.t[0:5],Sol2.y[2,0:5],'b',label='Temperatura (K)')
plt.plot(th3,Sol3[:,0],'b',label='Temperatura (K)')
plt.plot(th4,Sol4[:,0],'b',label='Temperatura (K)')
plt.plot(th5,Sol5[:,0],'b',label='Temperatura (K)')
plt.xlabel('Posición angular (Rad)')
plt.ylabel('Temperatura (K)')