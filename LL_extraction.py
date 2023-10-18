# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 16:00:34 2021

@author: Andrea Milazzo
"""

import numpy as np
from scipy.optimize import fsolve
import scipy.integrate as integrate
import matplotlib.pyplot as plt
plt.close('all')
a=49.17
b=34.17
a1=37.88
b1=20.25
curva1=lambda x: -2.799*((x-a)/b)**5+1.292*((x-a)/b)**4+10.56*((x-a)/b)**3-28.27*((x-a)/b)**2-8.887*((x-a)/b)+54.56 #curva di equilibrio sono in %
curva2=lambda x: 0.7507*((x-a1)/b1)**5+2.717*((x-a1)/b1)**4+4.207*((x-a1)/b1)**3+6.619*((x-a1)/b1)**2+12.52*((x-a1)/b1)+11.38 #curva di coniugazione interna
x1=np.linspace(0,100,101)
x2=np.linspace(0,61.4657,61)
plt.figure(1)
plt.plot(x1,curva1(x1)*np.cos(np.pi/6),color='b') #il x cos è per fare il grafico in coordinate xy
plt.plot(x2,curva2(x2)*np.cos(np.pi/6),color='r')
plt.plot([0,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([100,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([0,100],[0,0],color='k')
int1=lambda x: curva1(x)-curva2(x)
sol=fsolve(int1,[60]) #Troviamo l'intersezione tra le due curve
# print(sol) 
npun=10
xvet=np.linspace(1,59,npun)
x1=np.zeros(npun)
x2=np.zeros(npun)
x1[-1]=95
x2[-1]=70
y1=np.zeros(npun)
y2=np.zeros(npun)
int1=lambda x,y: y-curva1(x)*np.cos(np.pi/6) #troviamo intersezione orizzontale ( retta che va a destra partendo dalla curva di coniugazione)
int2=lambda x,xp,yp: -np.sqrt(3)*(x-xp)+yp -curva1(x)*np.cos(np.pi/6) #troviamo intersezione obliqua ( retta che sale partendo da curva di coniugazione)
for i in range(0,npun):
    y1[i]=curva2(xvet[i])*np.cos(np.pi/6)
    x1[i]=fsolve(int1,x1[i-1],args=(y1[i]))
    x2[i]=fsolve(int2,x2[i-1],args=(xvet[i],y1[i]))
    y2[i]=curva1(x2[i])*np.cos(np.pi/6)
    plt.plot([x1[i],x2[i]],[y1[i],y2[i]],color='k',linewidth='0.5')

# plt.xlim([0,100])
# plt.ylim([0,87])



#singolo stadio
r1=lambda x: np.sqrt(3)*x #retta caprolattame acqua
r3=lambda x: -np.sqrt(3)*(x-100) #retta caprolattame benzene
r2=lambda x,x0,y0: -np.sqrt(3)*(x-x0)-y0 #retta per unire il punto dalla curva 2 a r1 x0 e y0 in % (non frazioni)
int3=lambda x,y0: y0-curva2(x)*np.cos(np.pi/6) #intersezione tra il punto del benzene e curva 2 andando orizzontali, valore della X
int4=lambda x,x0,y0: r1(x)-r2(x,x0,y0) #valore della X  sulla retta caprolattame acqua (da portare in valore Y)
htr=100*np.cos(np.pi/6)
fconv=lambda y:y/100 * htr #da tri a xy
finv=lambda y:y/htr * 100


F=3000 #portata da trattare
PMBenz=6*12+6
PMCap=113.16
F=F/(0.25*PMCap+0.75*PMBenz)
xf=0.25
#xe=0.14 #fatto graficamente
#xeH2O=0.86
xr=0.02
def trovaEQ(xr): #do la composizione di raffinato e mi rende la composizione di capro nell estratto in %
    xr=fconv(xr)
    P1x=fsolve(int3,10,args=(xr*100)) #trovo la x
    P2x=fsolve(int4,20,args=(P1x,xr*100)) #trovo la y che poi andrà convertita
    return finv(r1(P2x))

[xe]=trovaEQ(xr)/100

xeH2O=1-xe
#Singolo stadio
A=np.array([[xr,xe,0],
            [1,1,-1],
            [0,xeH2O,-1]])
B=np.array([[F*xf,F,0]]).T
x=np.linalg.solve(A,B)
print('SINGOLO STADIO DI ESTRAZIONE')
print('Singolo stadio di estrazione a partire da 3000 Kg/h di miscela da tratttare')
print('la portata di raffinato è di: {:.2f} [kg/h]'.format(x[0,0]*(xr*PMCap+(1-xr)*PMBenz)))
print('la portata di estratto è di: {:.2f} [kg/h]'.format(x[1,0]*(xeH2O*18+PMCap*(xe))))
print('la portata di solvente è di: {:.2f} [kg/h]'.format(x[2,0]*18))
print('Singolo stadio di estrazione a partire da {:.2f} kmoli/h di miscela da tratttare'.format(F))
print('la portata di raffinato è di: {:.2f} [kmol/h]'.format(x[0,0]))
print('la portata di estratto è di: {:.2f} [kmol/h]'.format(x[1,0]))
print('la portata di solvente è di: {:.2f} [kmol/h]'.format(x[2,0]))
print('---------------------------------')

# Parte 2 cross-flow

S=F #solvente (acqua)
def trovaCROSS(yr,yfeed):
    [ye]=trovaEQ(yr) #troviamo 
    yr=fconv(yr)
    ye=fconv(ye)
    yfeed=fconv(yfeed)
    xe=ye/np.sqrt(3)
    xr=(-yr/np.sqrt(3)) +100
    m4=(yr-ye)/(xr-xe)
    r4=lambda x,x0,y0,m: m*(x-x0)-y0
    xfeed=(-yfeed/np.sqrt(3)) +100
    int5=lambda x,x0,y0,m,xf,yf: (yf/xf) *x-r4(x,x0,y0,m)
    xincr=fsolve(int5,50,args=(xe,ye,m4,xfeed,yfeed))
    yincr=(yfeed/xfeed) *xincr
    l1=2*np.sqrt(xincr**2 + yincr**2)
    ltot=np.sqrt(xfeed**2 + yfeed**2)
    return l1-ltot

# zzzz=fsolve(trovaCROSS,10,args=(25))
# zzzz=finv(zzzz)


def ControCorrente(xfeed):
    freteq= lambda x,m,x0,y0: m*(x-x0)+y0
    ffeed= lambda x,m: m*x
    finter= lambda x,meq,mf,x0,y0: freteq(x,meq,x0,y0)- ffeed(x,mf)
    mfeed=(xfeed*np.cos(np.pi/6))/(100-xfeed*np.sin(np.pi/6))
    xvet=59
    xo1=95
    xo2=70
    xbuono=((xfeed*0.5)*np.cos(np.pi/6))/mfeed
    xp=100
    while xp-xbuono>0:
        y1=curva2(xvet)*np.cos(np.pi/6)
        x1=fsolve(int1,xo1,args=(y1)) #intersezione orizzontale
        x2=fsolve(int2,xo2,args=(xvet,y1)) #intersezione obliqua
        y2=curva1(x2)*np.cos(np.pi/6)
        xo1=x1
        xo2=x2
        meq=(y1-y2)/(x1-x2)
        xp=fsolve(finter,50,args=(meq,mfeed,x1,y1))
        xvet=xvet-0.1
    return x1,x2

xfeed=25
i=0
xp1=np.zeros(10)
xp2=np.zeros(10)
while xfeed/np.cos(np.pi/6)>2:
    [xp1[i],xp2[i]]=ControCorrente(xfeed)
    xfeed=curva1(xp1[i])
    i=i+1

ftox=lambda y:100-y*np.sin(np.pi/6)
    
plt.figure(2)
x1=np.linspace(0,100,101)
x2=np.linspace(0,61.4657,61)
plt.plot(x1,curva1(x1)*np.cos(np.pi/6),color='b')#il x cos è per fare il grafico in coordinate xy
plt.plot(x2,curva2(x2)*np.cos(np.pi/6),color='r')
plt.plot([0,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([100,50],[0,100*np.cos(np.pi/6)],color='k')

plt.plot([0,ftox(25)],[0,25*np.cos(np.pi/6)])
plt.plot([0,xp1[0]],[0,curva1(xp1[0])*np.cos(np.pi/6)])
plt.plot([xp1[0],xp2[0]],[curva1(xp1[0])*np.cos(np.pi/6),curva1(xp2[0])*np.cos(np.pi/6)])
plt.plot([xp1[1],xp2[1]],[curva1(xp1[1])*np.cos(np.pi/6),curva1(xp2[1])*np.cos(np.pi/6)])

plt.xlim([0,100])
plt.ylim([0,87])

#counterflow
npun=20
xvet=np.linspace(5,52.25,npun) #abbiamo messo 52.25 perchè abbiamo visto che ci dava il punto P più distante
x11=np.zeros(npun)
x22=np.zeros(npun)
x11[-1]=95
x22[-1]=70
y11=np.zeros(npun)
y22=np.zeros(npun)
m=np.zeros(npun)
xP=np.zeros(npun)
yP=np.zeros(npun)

ftoP=lambda x: x*(2*np.cos(np.pi/6)/ftox(2))
fintP= lambda x,m,x0,y0: ftoP(x)-rP(x,m,x0,y0)
rP= lambda x,m,x0,y0: m*(x-x0)+y0
x1=np.linspace(0,100,101)
x2=np.linspace(0,61.4657,61)
plt.figure(3)
plt.plot(x1,curva1(x1)*np.cos(np.pi/6),color='b',label='Curva di Equilibrio')#il x cos è per fare il grafico in coordinate xy
plt.plot(x2,curva2(x2)*np.cos(np.pi/6),color='r',label='Curva di coniungazione interna')
plt.plot([0,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([100,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([ftox(25)],[25*np.cos(np.pi/6)],'*',label='Feed')
plt.plot([0,200],[0,ftoP(200)],label='Retta $R_n$ - $E_{n+1}$')

for i in range(0,npun):
    y11[i]=curva2(xvet[i])*np.cos(np.pi/6)              #y della intersezione orizzontale
    x11[i]=fsolve(int1,x11[i-1],args=(y11[i]))          #x della intersezione orizzontale
    x22[i]=fsolve(int2,x22[i-1],args=(xvet[i],y11[i]))  #x intersezione obliqua
    y22[i]=curva1(x22[i])*np.cos(np.pi/6)               #y intersezione obliqua
    m[i]=(y11[i]-y22[i])/(x11[i]-x22[i])
    xP[i]=fsolve(fintP,130,args=(m[i],x11[i],y11[i]))
    yP[i]=ftoP(xP[i])
    plt.plot([x22[i],xP[i]],[y22[i],yP[i]],color='k',linewidth='0.5') #stiamo disegnando tutte le tie line
maxind=np.argmax(xP)
plt.plot(xP[maxind],yP[maxind],'*',label='$P_{min}$')
plt.xlim([0,170])
plt.ylim([0,87])
plt.legend()

mtest1=(y22[-1]-2*np.cos(np.pi/6))/(x22[-1]-ftox(2))    #il 2 sarebbe la % di capro nel raffinato, coefficiente angolare retta Estratto1-Raffinato_finale
mtest2=25*np.cos(np.pi/6)/ftox(25)                      #coefficente angolare retta operativa
ftest=lambda x,m1,x1,y1,m2,x2,y2: rP(x,m1,x1,y1)-rP(x,m2,x2,y2)

[Pm]=fsolve(ftest,50,args=(mtest1,x22[-1],y22[-1],mtest2,0,0))
Pym=Pm*mtest2


plt.figure(4)
plt.plot(x1,curva1(x1)*np.cos(np.pi/6),color='b',label='Curva di Equilibrio')
plt.plot(x2,curva2(x2)*np.cos(np.pi/6),color='r',label='Curva di coniungazione interna')
plt.plot([0,50],[0,100*np.cos(np.pi/6)],'k')
plt.plot([100,50],[0,100*np.cos(np.pi/6)],'k')
plt.plot([ftox(25)],[25*np.cos(np.pi/6)],'o',color='y',label='F') #composizione xF
plt.plot([x22[-1],ftox(2)],[y22[-1],2*np.cos(np.pi/6)],'g')
plt.plot([0,ftox(25)],[0,25*np.cos(np.pi/6)],'b')

plt.plot([Pm],[Pym],'o',color='r',label='$M_{min}$')
plt.xlim([0,140])
plt.ylim([0,87])
plt.plot(0.1,0.1,'o',color='m',label='S')


mMin=(np.sqrt((Pm-ftox(25))**2+(Pym-fconv(25))**2))/(np.sqrt(Pm**2+Pym**2)) #sarebbe il rapporto tra la distanza tra FM e  la distanza MS

SolventeMin=F*mMin

print('ESTRAZIONE IN CONTROCORRENTE')
print('la portata di solvente minima per controcorrente è di: {:.2f} [kg/h] - {:.2f} [kmol/h]'.format(SolventeMin*18,SolventeMin))

SolventeIn=1.33*SolventeMin
M=SolventeIn+F
yM=((F*0.25)/(M))*np.cos(np.pi/6)*100  #cartesiano caprolattame

xM=yM/mtest2
plt.plot([xM],[yM],'o b',label='$M_{lavoro}$') #composizione con solventemin*1.3

mtest3=(yM-2*np.cos(np.pi/6))/(xM-ftox(2)) #pendenza retta che passa per M e per composizione raffinato finale
int6 =lambda x,xM,yM,M: curva1(x)*np.cos(np.pi/6)-(M*(x-xM)+yM)
[Xest]=fsolve(int6,25,args=(xM,yM,mtest3))
Yest=curva1(Xest)*np.cos(np.pi/6)           #cartesiano
plt.plot([Xest],[Yest],'*r')                #estratto ennesimo

x_water_est=Yest/(2*np.sqrt(3))+Xest*0.5
y_water_est=x_water_est*np.sqrt(3)          #coordinate cartesiane
comp_ben_est = 100-Yest/np.cos(np.pi/6)-(100-y_water_est/np.cos(np.pi/6))

print('la composizione dell estratto è di: {:.2f} % di caprolattame, {:.2f} % di acqua, {:.2f} % di benzene'.format(y22[-1]/np.cos(np.pi/6),100-y_water_est/np.cos(np.pi/6),comp_ben_est))
# plt.plot([x_water_est],[y_water_est],'*')

Estratto=100*SolventeIn/(100-y_water_est/np.cos(np.pi/6))
RaffinatoOut=F+SolventeIn-Estratto

mtest4=(Yest-25*np.cos(np.pi/6))/(Xest-ftox(25))
[xinterR]=fsolve(fintP,150,args=(mtest4,Xest,Yest))
yinterR=ftoP(xinterR)
plt.plot([0,xinterR],[0,yinterR])
plt.plot([Xest,xinterR],[Yest,yinterR])
plt.plot([Xest],[Yest],'o',color='g',label='$E_1$')
plt.legend()


int7= lambda x,x0,y0: curva2(x)*np.cos(np.pi/6)+np.sqrt(3)*(x-x0)-y0
int8= lambda x,y0: curva1(x)*np.cos(np.pi/6)-y0
int9= lambda x,m : curva1(x)*np.cos(np.pi/6)-(yinterR+m*(x-xinterR))
xE1=np.zeros(20)
xR1=np.zeros(20)
i=0
yE1=np.zeros(20) #frazione molare caprolattame nell'estratto ----> ordinate in grafico 6
yR1=np.zeros(20) #frazione molare caprolattame nel raffinato ----> ascisse in grafico 6

xE1[0]=Xest
yE1[0]=Yest
yR1[-1]=100
xR1[-1]=90

plt.figure(5)
plt.plot(x1,curva1(x1)*np.cos(np.pi/6),color='b',label='Curva Equilibrio')#il x cos è per fare il grafico in coordinate xy
plt.plot(x2,curva2(x2)*np.cos(np.pi/6),color='r',label='Curva Coniugazione interna')
plt.plot([0,50],[0,100*np.cos(np.pi/6)],color='k')
plt.plot([100,50],[0,100*np.cos(np.pi/6)],color='k')
plt.xlim([0,100])
plt.ylim([0,87])
x7=30

while yR1[i-1]>2*np.cos(np.pi/6):
    [x7]=fsolve(int7,30,args=(xE1[i],yE1[i]))
    y7=curva2(x7)*np.cos(np.pi/6)
    
    xR1[i]=fsolve(int8,xR1[i-1],args=(y7))
    yR1[i]=curva1(xR1[i])*np.cos(np.pi/6)
    
    m=(yR1[i]-yinterR)/(xR1[i]-xinterR)
    
    xE1[i+1]=fsolve(int9,xE1[i],args=(m))
    yE1[i+1]=yinterR+m*(xE1[i+1]-xinterR)
    if yR1[i]>2*np.cos(np.pi/6):
        plt.plot([xE1[i],xR1[i]],[yE1[i],yR1[i]],color='m')           #tie line
        plt.plot([xE1[i+1],xR1[i]],[yE1[i+1],yR1[i]],color='g')       #costruzione
    
    i=i+1
plt.plot([ftox(25),xE1[0]],[25*np.cos(np.pi/6),yE1[0]],color='g')    
plt.legend()

plt.figure(6)
plt.plot([0,1],[0,1])


npun=100
xvet=np.linspace(1,sol,npun)
x1=np.zeros(npun)
x2=np.zeros(npun)
y1=np.zeros(npun)
y2=np.zeros(npun)
x1[-1]=95
x2[-1]=70

for i in range(0,npun):
    y1[i]=curva2(xvet[i])*np.cos(np.pi/6)               #raffinato
    x1[i]=fsolve(int1,x1[i-1],args=(y1[i]))             #raffinato
    x2[i]=fsolve(int2,x2[i-1],args=(xvet[i],y1[i]))     #estratto
    y2[i]=curva1(x2[i])*np.cos(np.pi/6)                 #estratto
    
plt.plot(np.append(0,y1[1:npun]/100)/np.cos(np.pi/6),np.append(0,y2[1:npun]/100)/np.cos(np.pi/6),color='b',label='Curva Equilibrio')
plt.plot(0.25,yE1[0]/100/np.cos(np.pi/6),'r o')

for i in range(0,4):
    plt.plot(yR1[i]/100/np.cos(np.pi/6),yE1[i+1]/100/np.cos(np.pi/6),'r o')
    plt.plot([yR1[i]/100/np.cos(np.pi/6),yR1[i+1]/100/np.cos(np.pi/6)],[yE1[i+1]/100/np.cos(np.pi/6),yE1[i+1]/100/np.cos(np.pi/6)],color='k')
    plt.plot([yR1[i]/100/np.cos(np.pi/6),yR1[i]/100/np.cos(np.pi/6)],[yE1[i+1]/100/np.cos(np.pi/6),yE1[i]/100/np.cos(np.pi/6)],color='k')
    

plt.plot(np.append(0.25,yR1[0:4]/100/np.cos(np.pi/6)),yE1[0:5]/100/np.cos(np.pi/6),color='r',label='Retta Operativa')

plt.plot([0.25,yR1[0]/100/np.cos(np.pi/6)],[yE1[0]/100/np.cos(np.pi/6),yE1[0]/100/np.cos(np.pi/6)],color='k')

plt.xlim([0,0.5])
plt.ylim([0,0.6])
plt.grid()
plt.legend()
plt.ylabel('Estratto')
plt.xlabel('Raffinato')

#valutare unita trasferimento
xE8=yE1[0]/np.cos(np.pi/6)
xE8=Estratto*xE8*PMCap/(Estratto*xE8*PMCap+Estratto*(1-xE8)*18)                     #### frazioni in peso in questa sezione
xR8=RaffinatoOut*0.02*PMCap/(RaffinatoOut*0.02*PMCap+RaffinatoOut*0.98*PMBenz)
xF8=0.25*PMCap/(0.25*PMCap+0.75*PMBenz)
K8=xE8/xR8
eps=SolventeIn*K8*18/(F*(0.25*PMCap+(1-0.25)*PMBenz))
pippo8= lambda n: ((eps**(n+1)-eps)/(eps**(n+1)-1))-(xF8-xR8)/(xF8)

# ciao=fsolve(pippo8,4)
# print(ciao)

# CURVE    #####################################################################################################################
Curva_Eq= lambda x: -642.6*x**6+1079*x**5-731.7*x**4+258.8*x**3-54.15*x**2+7.38*x-0.004038 #curva equilibrio ye=f(yr)

Curva_Op= lambda x: 39.95*x**4-5.913*x**3-7.6*x**2+3.967*x-0.07612 #curva operativa ye=f(yr)
integ= lambda x: 1/(Curva_Eq(x)-Curva_Op(x)) #stessa equazione di pippo100 ma integrata sulle y perche era piu comodo

#Siccome ci sono dei piccoli errori di precisione in come sono stati presi i punti dalla Curva_Op bisogna stare attenti che Curva_Op(x) non sia mai minore di 0
#altrimenti mi combina un casino quando fa integrale e logaritmo
#cerco quando becco lo zero in Curva_Op(x)

[x_0]=fsolve(Curva_Op,0.05) #zero della curva operativa

vector=np.array([yR1[0],yR1[1],yR1[2],yR1[3],x_0*100])/100
vector=vector/np.cos(np.pi/6)
NTU2=0

# for i in range (len(vector)-1):
#     [integrale,blabla]=integrate.quad(integ,Curva_Op(vector[i]),Curva_Op(vector[i+1]))
#     NTU2=integrale+NTU2+0.5*np.log((1-Curva_Op(vector[i]))/(1-Curva_Op(vector[i+1])))
#     print(NTU2)

# [integrale,blabla]=integrate.quad(integ,Curva_Op(vector[0]),Curva_Op(vector[-1]))
# NTU2=integrale+0.5*np.log((1-Curva_Op(vector[0]))/(1-Curva_Op(vector[-1])))


#Qui sotto è stato fatto insieme ########################################################################

pippo100= lambda x: 1/(Curva_Eq(x)-Curva_Op(x))
# [pippoint,zzzzzzz345a4edhse]=integrate.quad(pippo100,Curva_Op(0.25),0)

# NTU=pippoint+0.5*np.log((1-Curva_Op(0.25))/(1-0))

# print('Number of Theoretical Units  =',NTU)


####### TEST 3

vecxr=np.linspace(0.02,0.25,50)

vecyeq=np.zeros(50)
vecxeq=np.zeros(50)

vecyop=Curva_Op(vecxr)

for i in range (0,len(vecxr)):
    int14 = lambda x,x0,y0: Curva_Eq(x)-(-x+x0+y0)
    [a]=fsolve(int14,vecxr[i]*0.8,args=(vecxr[i],vecyop[i]))
    vecyeq[i]=Curva_Eq(a)
    vecxeq[i]=a

vecint=1/(vecyeq-vecyop)
vecint2=1/(vecxr-vecxeq)

trap_Integrale=np.trapz(vecint,vecxr)
trap_Integrale2=np.trapz(vecint2,vecxr)

NTU=trap_Integrale
print('NTU = %.2f' %NTU)

ComposizioniE=np.zeros([4,3]) # Caprolattame, acqua, benzene
ComposizioniR=np.zeros([4,3]) # Caprolattame, acqua, benzene
Portate=np.zeros([5,2]) #Estratto, raffinato
Portate[0,1]=F
Portate[4,0]=SolventeIn



#   estratto
for i in range(4):
    ComposizioniE[i,0]=finv(yE1[i]) #caprolattame
    interc = xE1[i]+yE1[i]/np.tan(np.pi/3)
    ComposizioniE[i,1]=100-interc #acqua
    ComposizioniE[i,2]=100-ComposizioniE[i,0]-ComposizioniE[i,1] #benz
    
    ComposizioniR[i,0]=finv(yR1[i]) #caprolattame
    interc = xR1[i]+yR1[i]/np.tan(np.pi/3)
    ComposizioniR[i,1]=100-interc #acqua
    ComposizioniR[i,2]=100-ComposizioniR[i,0]-ComposizioniR[i,1] #benz

ComposizioniR=np.vstack([[xf*100,0,100-xf*100],ComposizioniR]) 
ComposizioniE=np.vstack([ComposizioniE,[0,100,0]])   
#Bilancio Globale
A=np.array([[ComposizioniR[4,0]/100,ComposizioniE[0,0]/100],
            [1,1]])
B=np.array([[F*xf,F+SolventeIn]]).T
[Portate[4,1],Portate[0,0]]=np.linalg.solve(A,B)

    
for i in range(3): #numero 0 parte da ingresso Feed
        A=np.array([[-ComposizioniR[i+1,0]/100,+ComposizioniE[i+1,0]/100],
                    [-1,+1]])
        B=np.array([[Portate[i,0]*ComposizioniE[i,0]/100-Portate[i,1]*ComposizioniR[i,0]/100,Portate[i,0]-Portate[i,1]]]).T
        [Portate[i+1,1],Portate[i+1,0]]=np.linalg.solve(A,B)
    
print('PRIMO STADIO')
i=0
print('Estratto: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniE[i,0],ComposizioniE[i,1],ComposizioniE[i,2],Portate[i,0]))
print('Raffinato: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniR[i,0],ComposizioniR[i,1],ComposizioniR[i,2],Portate[i,1]))

print('SECONDO STADIO')
i=1
print('Estratto: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniE[i,0],ComposizioniE[i,1],ComposizioniE[i,2],Portate[i,0]))
print('Raffinato: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniR[i,0],ComposizioniR[i,1],ComposizioniR[i,2],Portate[i,1]))

print('TERZO STADIO')
i=2
print('Estratto: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniE[i,0],ComposizioniE[i,1],ComposizioniE[i,2],Portate[i,0]))
print('Raffinato: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniR[i,0],ComposizioniR[i,1],ComposizioniR[i,2],Portate[i,1]))

print('QUARTO STADIO')
i=3
print('Estratto: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniE[i,0],ComposizioniE[i,1],ComposizioniE[i,2],Portate[i,0]))
print('Raffinato: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniR[i,0],ComposizioniR[i,1],ComposizioniR[i,2],Portate[i,1]))

print('QUINTO STADIO')
i=4
print('Estratto: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniE[i,0],ComposizioniE[i,1],ComposizioniE[i,2],Portate[i,0]))
print('Raffinato: \tComp = {:.2f} %C {:.2f} %W {:.2f} %B \tPortata = {:.2f} [kmol/h]'.format(ComposizioniR[i,0],ComposizioniR[i,1],ComposizioniR[i,2],Portate[i,1]))










