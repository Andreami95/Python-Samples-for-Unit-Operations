# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:52:19 2021

@author: Davide Carli e Andrea Milazzo
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.integrate as integrate
plt.close('all')


valsat=np.array([ 7.36490e+01, -7.25820e+03, -7.30370e+00,  4.16530e-06,2.00000e+00,  2.73160e+02,  6.47096e+02])
valcpgACQUA=np.array([3.33630e+04, 2.67900e+04, 2.61050e+03, 8.89600e+03, 1.16900e+03,1.00000e+02, 2.27315e+03])
valcpgARIA=np.array([2.8958e+04, 9.3900e+03, 3.0120e+03, 7.5800e+03, 1.4840e+03,5.0000e+01, 1.5000e+03])
vallambda=np.array([5.2053e+07, 0.3199, -0.212, 0.25795])
valmul=np.array([-52.843, 3703.6, 5.866, -5.879e-29, 10])
valrol=np.array([-13.851, 0.64038, -0.00191, 1.8211e-6])
valmug=np.array([1.425e-6,0.5039,108.3])
valktl=np.array([-0.432, 0.0057255, -0.000008078, 1.861e-9])
valPsat=np.array([73.649,-7258.2,-7.3037,4.1653e-6,2])
Tcr=647.096

CP=lambda T,valcpg:valcpg[0]+valcpg[1]*((valcpg[2]/T)/(np.sinh(valcpg[2]/T)))**2 +valcpg[3]*((valcpg[4]/T)/(np.cosh(valcpg[4]/T)))**2
LAMBDA= lambda T: vallambda[0]*(1-T/Tcr)**(vallambda[1]+vallambda[2]*T/Tcr +vallambda[3]*(T/Tcr)**2)
RO= lambda T: valrol[0]+valrol[1]*T + valrol[2]*T**2 +valrol[3]* T**3
MUl= lambda T: np.exp(valmul[0]+ valmul[1]/T +valmul[2]*np.log(T) + valmul[3]*T**valmul[4])
ROg= lambda T:P*pmar/(8.314*T*1000)
MUg= lambda T: (valmug[0]*T**valmug[1])/(1+valmug[2]/T)
DIFFg= lambda T: 22.5e-6*(T/(273.15))**1.8 #rif a sentimento =)
KTER= lambda T:valktl[0]+valktl[1]*T+valktl[2]*T**2 +valktl[3]*T**3
PSAT= lambda T: np.exp(valPsat[0]+valPsat[1]/T + valPsat[2]*np.log(T) + valPsat[3]*T**valPsat[4])/(1e5)
pma=18
pmar=28.96

def Hg(T,xp):
    [Cpgaria,e]=integrate.quad(CP,273.15,T,args=(valcpgARIA))
    [Cpgacqua,e]=integrate.quad(CP,273.15,T,args=(valcpgACQUA))
    [LAM,e]=integrate.quad(LAMBDA,273.15,T)
    Cpgaria=Cpgaria/(1000*pmar*(T-273.15))
    Cpgacqua=Cpgacqua/(1000*pma*(T-273.15))
    LAM=LAM/(1000*pma*(T-273.15))
    return (Cpgaria+Cpgacqua*xp)*(T-273.15)+(xp)*LAM

def CPm(T0,T1,valcpg,PM):
    [Cpg,e]=integrate.quad(CP,T0,T1,args=(valcpg))
    return Cpg/(1000*PM*(T1-T0))

Tl2=50+273.15 #ingresso
Tl1=30+273.15 #uscita
Tg1=35+273.15 #ingresso
Hsin=0.0043 #a 15 0.0041 umidità assoluta gas kg H2O su kg aria secca
L=10000 #kg/h
Gminimo = 1600 # minimo con 15 C aria 50-25 h20 2145 # kg/H di aria
G=Gminimo*1.3
hr=24.8e-3
Br=55e-3
Sr=37e-3
alphar=45*(np.pi/180)
a=125 #m2/m3
eps=0.95
Pack_factor=33
Dc=0.6
Ac=Dc**2 *0.25*np.pi
P=1e5
viscoH2O= 8.926573449635064e-07
Tm=(Tl1+Tl2)/2
# CHECK WETTING RATE
Lv=L/(RO(Tm)*pma)       #m^3/h
wr=Lv/(Ac*a)            #m^3/h/m  
if wr<0.02:
    print('wr is too low')
    
FLG=L/G*(ROg(Tm)/(RO(Tm)*pma))**0.5 #guarda kister pag 490 grafico B
Cs=(G/(Ac*3600)/(ROg(Tm)**0.5*(RO(Tm)*pma-ROg(Tm))**0.5))
us_gas=G/ROg(Tm)/3600/Ac

CPgraf = us_gas*((Pack_factor)**(0.5)*(ROg(30)*pmar/(1000-1.16))**0.5*(viscoH2O*1e6)**0.05)*0.3048          #m/s VELOCITA SUPERFICIALE GAS


#CHECK F-FACTOR questi valori limite sono stati dati a lezione
F_Factor=us_gas*ROg(Tm)**0.5
if F_Factor<1.5:
    print('F-Factor is too low')
elif F_Factor>2.8:
    print('F-Factor is too high')
print('F-Factor =','%.2f'  %F_Factor, 'Pa^0.5')
#CHECK VELOCITA' SUPERFICIALE LIQUIDA questi valori limite sono stati dati a lezione
us_liq=(L/(RO(Tm)*pma))/Ac #m3/h/m2


ae=a*1.34*(((RO(Tm)*pma)/(72*1e-3))*9.81**(1/3)*((L/(RO(Tm)*pma*3600))/(4*Sr*Ac/(Br*hr)))**(4/3))**0.116
if ae<a:
    a=ae

curvaTH= lambda T: 1.202*((T-305.6)/13.27)**4 +6.106*((T-305.6)/13.27)**3 +22.79*((T-305.6)/13.27)**2 +75.79*((T-305.6)/13.27) + 112.2

Hg1=Hg(Tg1,Hsin) #entalpia gas in entrata
# plt.figure(1)
# plt.grid()
# plt.plot(Tg1-273,Hg1,'*')

Hg2=L * CPm(Tl1,Tl2,valcpgACQUA,pma)* (Tl2-Tl1)/(G) +Hg1 #entalpia gas in uscita
Hg2min=L * CPm(Tl1,Tl2,valcpgACQUA,pma)* (Tl2-Tl1)/(Gminimo) +Hg1

# plt.plot(Tl2-273,Hg2,'*')
# plt.plot(Tl1-273,Hg1,'*')
# plt.plot([Tl1-273,Tl2-273],[Hg1,Hg2],label='Retta di Lavoro')

Hl2=Hg2 #acqua in
Hl1=Hg1 #acqua out

R12= lambda H: (H-(-((Hg2-Hg1)/(Tl2-Tl1))*Tl2 + Hg2))/((Hg2-Hg1)/(Tl2-Tl1))

R= lambda T,T0,H0: (-hl(T)/(kgv(T)*ROg(T)))*T+((hl(T)/(kgv(T)*ROg(T))))*T0+H0

# R2= lambda T,T0,H0: (-hlt(T)/(kgv(T)*RO(T)*pma))*T+((hlt(T)/(kgv(T)*RO(T)*pma)))*T0+H0

R3= lambda T,T0,H0: (-hlt(T)/(kgv(T)*ROg(T)))*T+((hlt(T)/(kgv(T)*ROg(T))))*T0+H0 #T0 e H0 sono di interfaccia, questa funzione si usa per trovare H gas, T è del liquido: è la curva 'r'

hl= lambda T: (0.0169*a**0.83 *(L/(RO(T)*pma*Ac))**0.73*(MUl(T)/MUl(273.15+20))**0.25)/100 #hold up liquido

hlt= lambda T: (kgv(T)*KTER(T))/DIFFg(T)
K_aria=0.028

Uge= lambda T:(G/(ROg(T)*Ac*3600))/(eps*(1-hl(T))*np.sin(alphar))
Ule= lambda T:(L/(RO(T)*pma*Ac*3600))/((eps*hl(T))*np.sin(alphar))

kgv= lambda T: DIFFg(T)/Sr *0.054*((Uge(T)+Ule(T))*ROg(T)*Sr/MUg(T))**0.8 * (MUg(T)/(DIFFg(T)*ROg(T)))**0.33 #m/s

solveeq= lambda T,T0,H0: R3(T,T0,H0)- curvaTH(T)
[solT]=fsolve(solveeq,Tl1,args=(Tl1,Hl1))

# plt.plot([Tl1-273,solT-273],[Hl1,curvaTH(solT)])


# Tsp=np.linspace(10+273.15,55+273.15,70)
# Hsp=curvaTH(Tsp)
# Tsp=Tsp-273
# plt.plot(Tsp,Hsp,label='Curva equilibrio')

npunti=100
Hvet=np.linspace(Hl1,Hl2,npunti)
Tgvet=np.zeros(npunti)
Tgvet[0]=Tg1
Tlvet=np.zeros(npunti)
Tlvet[0]=Tl1
Tsatvet=np.zeros(npunti)
Hsatvet=np.zeros(npunti)
Hmed=np.zeros(npunti)
Hass=np.zeros(npunti)
dH=np.zeros(npunti)
dz=np.zeros(npunti)
Tgmed=np.zeros(npunti)
Tlmed=np.zeros(npunti)
Tsatgas=np.zeros(npunti)
Tsatmed=np.zeros(npunti)
Tsat2=np.zeros(npunti)
MH2O=np.zeros(npunti)
yAIR=np.zeros(npunti)
yH2O=np.zeros(npunti)
PparzH2O=np.zeros(npunti)
NTU=np.zeros(npunti)
# Tgvetprova=np.zeros(npunti)
# Tgvetprova[0]=Tg1
Tgmed[0]=Tg1
Tlmed[0]=Tl1
Itot=0

for i in range(0,npunti-1):
    Tsatvet[i]=fsolve(solveeq,Tlvet[i],args=(Tlvet[i],Hvet[i])) # uso retta 'r' e vado dalla retta di lavoro verso la curva eq. e trovo T
    Hsatvet[i]=curvaTH(Tsatvet[i])*0.95
    Tgvet[i+1]=(Hvet[i+1]-(-(Hsatvet[i]-Hvet[i])/(Tsatvet[i]-Tgvet[i])*Tgvet[i]+Hvet[i]))/((Hsatvet[i]-Hvet[i])/(Tsatvet[i]-Tgvet[i])) #originale
    # Tgvetprova[i+1]=Tgvet[i]+((Hvet[i+1]-Hvet[i])/((Hsatvet[i]-Hvet[i])/(Tsatvet[i]-Tgvet[i])))     #mi sembra piú chiaro
    Tlvet[i+1]=R12(Hvet[i+1])
    # plt.plot([Tlvet[i]-273,Tsatvet[i]-273],[Hvet[i],Hsatvet[i]])
    # plt.plot([Tsatvet[i]-273,Tgvet[i]-273],[Hsatvet[i],Hvet[i]])
    # plt.plot([Tgvet[i+1]-273,Tlvet[i+1]-273],[Hvet[i+1],Hvet[i+1]])
    dH[i]=1/(Hsatvet[i]-Hvet[i])
    funt= lambda T,H : curvaTH(T)- H
    Tsatgas[i]=fsolve(funt,Tgvet[i],args=(Hvet[i]))
    if i>0:
        I=((dH[i]+dH[i-1])/2)*(Hvet[i]-Hvet[i-1])
        Tgmed[i]=(Tgvet[i]+Tgvet[i-1])/2
        Tlmed[i]=(Tlvet[i]+Tlvet[i-1])/2
        Tsatmed[i]=(Tsatgas[i]+Tsatgas[i-1])/2
        NTU[i+1]=I+NTU[i]
        dz[i]=dz[i-1]+I*G/(Ac*3600*ROg(Tgmed[i])*a*kgv(Tgmed[i]))
        Tsat2[i]=(Tsatvet[i]+Tsatvet[i-1])/2
        Itot=Itot+I
        Hmed[i]=(Hvet[i]+Hvet[i-1])/2
        funHass= lambda Ha,T,H : H-(Ha*LAMBDA(T)/(1000*pma) + (CP(T,valcpgARIA)/(1000*pmar)+Ha*CP(T,valcpgACQUA)/(1000*pma))*(T-273.15))
        Hass[i]=fsolve(funHass,Hass[i-1],args=(Tgmed[i],Hmed[i]))
        
    else:
        Tsatmed[i]=Tsatgas[i]
        Tsat2[i]=Tsatvet[i]
        Hmed[i]=Hvet[i]
        Hass[i]=Hsin


AltezzaColonna=dz[-2]

efficienza = 1 - np.exp(-max(NTU))
print('Efficienza','%.2f' %efficienza)
print('Altezza colonna è:', '%.3f' %AltezzaColonna,'m')
Tsp=np.linspace(10+273.15,55+273.15,70)
Hsp=curvaTH(Tsp)
Tsp=Tsp-273

plt.figure(1)
plt.grid()
plt.plot(Tsp,Hsp,label='Curva equilibrio',color='b')    
plt.plot([Tl1-273,Tl2-273],[Hg1,Hg2],label='Retta di Lavoro',color='r')   
plt.plot([Tl1-273,Tl2-273],[Hg1,Hg2min],label='Retta di Lavoro con $G_{min}$',color='g')   
plt.legend()
plt.ylabel('H (kj/kg aria secca)')
plt.xlabel('Temperatura (°C)')

plt.figure(2)
plt.grid()
Tsp=np.linspace(10+273.15,55+273.15,70)
Hsp=curvaTH(Tsp)
Tsp=Tsp-273
plt.plot(Tsp,Hsp,label='Curva equilibrio',color='b')    
plt.plot([Tl1-273,Tl2-273],[Hg1,Hg2],label='Retta di Lavoro',color='r')   
Tgvet=Tgvet-273
# Tgvetprova=Tgvetprova-273
plt.plot(Tgvet,Hvet,label='Temperatura gas',color='g')
plt.ylabel('H (kj/kg aria secca)')
plt.xlabel('Temperatura (°C)')
plt.plot(Tg1-273.15,Hg1,'o',label='Ingresso aria',color='y')
plt.plot(Tgvet[-1],Hg2,'o',label='Uscita aria',color='c')
plt.plot(Tl2-273,Hg2,'o',label='Ingresso $H_2O$',color='b')
plt.plot(Tl1-273,Hg1,'o',label='Uscita $H_2O$',color='m')
plt.legend()
# plt.plot([Tl1-273,solT-273],[Hl1,curvaTH(solT)])


Tgmed=Tgmed-273
Tlmed=Tlmed-273
Tsatmed=Tsatmed-273

plt.figure(3)
plt.grid()
plt.plot(Tgmed[0:-1],dz[0:-1],color='r',label='Temperatura aria')
plt.plot(Tlmed[0:-1],dz[0:-1],color='b',label='Temperatura $H_2O$')
plt.plot(0.5*(Tgmed[0:-1]+Tlmed[0:-1]),dz[0:-1],linestyle='--',label='Temperatura interfaccia',color='k',linewidth='0.7')
plt.ylabel('altezza (m)')
plt.xlabel('Temperatura (°C)')
plt.legend()
#plt.plot(Tsat2[0:-1],dz[0:-1],'*')


plt.figure(4)
plt.grid()
# plt.title(label='Entalpia vs altezza')
plt.xlabel('H (kj/kg aria secca)')
plt.ylabel('altezza (m)')
plt.plot(Hmed[0:-1],dz[0:-1])

plt.figure(5)
plt.grid()
# plt.title(label='Umidita assoluta vs altezza')
plt.ylabel('altezza (m)')
plt.xlabel('Umidità assoluta (Kg H20/Kg aria secca)')
plt.plot(Hass[0:-1],dz[0:-1])

plt.figure(6)
plt.grid()
# plt.title(label='Umidita Relativa vs altezza')
plt.ylabel('altezza (m)')
plt.xlabel('Umidità relativa')
PsatAcqua=np.zeros(npunti)

for i in range (0, npunti-1):
    PsatAcqua[i]=PSAT(Tlmed[i]+273.15)
    MH2O[i]=Hass[i]/(Hass[i]+1)
    yAIR[i]=1/29/(1/29+MH2O[i]/18)  #29 peso molecolare aria
    yH2O[i]=1-yAIR[i]
    PparzH2O[i]=yH2O[i]*1 #1 atm
    
plt.plot(PparzH2O[0:-1]/PsatAcqua[0:-1],dz[0:-1])


from thermo.chemical import Chemical
from thermo import Mixture
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

H2O_humidity        = lambda T:(Chemical('H2O',T=T,P=101325))
y_water             = lambda T: H2O_humidity(T).Psat/101325
UA100               = lambda T: y_water(T)*H2O.MW/((1-y_water(T))*air.MW)                                       #kgH2O/kg dry air

air=Mixture(['air'],T=303.15,P=101325)
H2O=Chemical('H2O') #ho messo 30C perchè faccio stima conservativa, so che la T si alzerà mentre si assorbe HCl

# Psychrometric Diagram #########################################
abs_humidity_in=np.zeros(501)
R_humidity90=np.zeros(501)
R_humidity80=np.zeros(501)
R_humidity70=np.zeros(501)
R_humidity60=np.zeros(501)
R_humidity50=np.zeros(501)
R_humidity40=np.zeros(501)
R_humidity30=np.zeros(501)
R_humidity20=np.zeros(501)
R_humidity10=np.zeros(501)

Tempvec=np.linspace(273.15,273.15+50,501)
# Tempvec=np.linspace(0,50,501)
for i in range (len(Tempvec)):
    H2Ov=Chemical('H2O',T=Tempvec[i])
    y_H2O=H2Ov.Psat/101325
    mass_H2Oin=y_H2O*H2O.MW
    y_air=1-y_H2O
    mass_air=y_air*air.MW
    abs_humidity_in[i]=mass_H2Oin/mass_air*1000
    R_humidity90[i]=mass_H2Oin/mass_air*1000*0.9
    R_humidity80[i]=mass_H2Oin/mass_air*1000*0.8
    R_humidity70[i]=mass_H2Oin/mass_air*1000*0.7
    R_humidity60[i]=mass_H2Oin/mass_air*1000*0.6
    R_humidity50[i]=mass_H2Oin/mass_air*1000*0.5
    R_humidity40[i]=mass_H2Oin/mass_air*1000*0.4
    R_humidity30[i]=mass_H2Oin/mass_air*1000*0.3
    R_humidity20[i]=mass_H2Oin/mass_air*1000*0.2
    R_humidity10[i]=mass_H2Oin/mass_air*1000*0.1

# intersectionHumidity = lambda x: abs_humidity_in-x
# [provabla] = fsolve(intersectionHumidity(5),290)
Tempvec=Tempvec-273.15
plt.figure(7)
plt.plot(Tempvec,abs_humidity_in,color='k',linewidth='0.7')
plt.plot(Tempvec,R_humidity90,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity80,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity70,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity60,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity50,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity40,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity30,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity20,color='k',linewidth='0.5')
plt.plot(Tempvec,R_humidity10,color='k',linewidth='0.5')
plt.xlim(0,50)
plt.ylim(0,90)
plt.xlabel('Temperatura (C)')
plt.ylabel('Umidità assoluta (g H2O/kg aria secca)')

prova=np.array([10,20,30,40,50,60,70,80])
for i in range(len(prova)):
    intersection2= lambda T: UA100(T)*1000-prova[i]
    [Tprova]=fsolve(intersection2,310)
    plt.plot([Tprova-273.15,50],[prova[i],prova[i]],linewidth='0.4',color='k')
# prova=np.linspace(0,90,19)
# for i in range(19): 
    
plt.plot(Tgmed[0:-1],Hass[0:-1]*1000,color='b')
plt.title('Diagramma Psicrometrico')
