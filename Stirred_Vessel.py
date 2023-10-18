# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:03:45 2022

@author: Andrea Milazzo
"""

from thermo.chemical import Chemical
from thermo import Mixture
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.close('all')

########################################## DATI INIZIALI ##############################################
yO2in=0.3                       #frazione molare di O2 in ingresso
P=1                             #atm
Tin=273.15+25                   #T iniziale dell'acqua nel reattore
Tfin=273.15+110                 #T finale dell'acqua nel reattore 

water=Chemical('Water')
air=Mixture(['N2','O2'],[1-yO2in,yO2in])


#######################################################################################
############################ REATTORE PILOTA ##########################################
#######################################################################################

############################# GEOMETRIA ###############################################
Vol_s=1.5                           #m3 Reattore piccolo

#Ipotizziamo inizialmente Z = T
T_s=((4*Vol_s)/np.pi)**(1/3)        #m diametro interno reattore
#Si parte con Z=0.75*T
#H è altezza tank
#Z è altezza liquido

Z_s=T_s                             #m altezza liquido
H_s=Z_s/0.75                        #m altezza reattore

#Turbina Rushton RT_6
D_s=1/3*T_s                         #m diametro turbina

#Distanza palette da fondo rettaore
c_s=1/3*T_s                         #m

#mettiamo 4 baffle
B_s=T_s/11                          #m larghezza baffle

Vr_s=0.25*np.pi*T_s**2*H_s          #m3 volume del reattore 

#Imposizione di vvm portata volumentrica di gas su volume reattore m3/min/m3 
vvm=0.5                             #m3/min/m3 
Q_g=(vvm/60)*Vol_s                  #m3/s

#bisogna calcolare Ncd e Nr, noi impondiamo diffusore ad anello quindi alpha = 3 
alpha=3
N_cd=(alpha*Q_g**0.5*T_s**0.25*D_s**(-2))*60            #rpm velocità minima
N_r=(1.5*Q_g**0.2*T_s*D_s**(-2))*60                     #rpm velocità massima

#SI SCEGLIE UN VALORE DI N_s compreso tra N_cd e N_r
N_s=200/60                          #rps
v_p=np.pi*N_s*D_s                   #m/s velocità periferica

#VERIFICA DELLE CONDIZIONI
g=9.81                              #m/s2
Fr_s=N_s**2*D_s/g                   #numero Froude
N_a=Q_g/(N_s*D_s**3)                #numero di aerazione                    

#Le coordinate dei punti sono (Fr,N_a)

########################## CHECK CONDIZIONI ###############################################

#CHECK NUMERO FROUDE
if Fr_s<0.04:
    raise ValueError('Fr pilot reactor < 0.04')

#CHECK FLOODING condizione N_a_fl>=N_a
N_a_fl=30*Fr_s*(D_s/T_s)**(3.5)

if N_a_fl<N_a:
    raise ValueError('N_a_fl<N_a')

# CHECK RICIRCULATION condizione N_a>=N_a_ric
N_a_ric=13*Fr_s**2*(D_s/T_s)**5

if N_a<N_a_ric:
    raise ValueError('N_a<N_a_ric')

#CHECK CAVITY LINE condizione Na>=N_a_cl
N_a_cl=0.025*(D_s/T_s)**(-0.5)

if N_a<N_a_cl:
    raise ValueError('N_a<N_a_cl')

#Se sono tutti soddisfatti sono nella zona 5 pag 64 slide (mixing slide)

############################ PORTATA MASSIMA (per non incorrere in flooding)  ###########################
Q_gmax=N_a_fl*N_s*3600*D_s**3                       #m3/h
Fl_percent=(Q_g*3600/(Q_gmax))*100                  #%

############################ POTENZA RICHIESTA ########################################
Re_s=water.rho*N_s*D_s**2/water.mu
Np_s=5.5                                            #NUMERO DI POTENZA si ottiene dal grafico

Pl_s=Np_s*water.rho*N_s**3*D_s**5                   # watt POTENZA UTILE LIQUIDO NECESSARIA

#vediamo se N_a è più grande di 0.037, nel nostro caso è più grande (pag. 58 slide)

Pgl_s=Pl_s*(0.62-1.85*N_a)                          # watt POTENZA UTILE GAS-LIQUIDO
vsg_s=4*Q_g/(np.pi*T_s**2)                          #m/s
kla_s=2.3*((Pgl_s/1000)/Vol_s)**0.7*(vsg_s)**0.6    #s^-1 coefficiente di trasporto liquido per area di interfaccia 

#ADESSO BISOGNA PORTARE A CONVERGENZA LE DUE EQUAZIONI DEL BILANCIO DI MATERIA E DEL FLUSSO MATERIALE

H=43400                                             #atm/(mol/l) costante di Henry
c_O2=2.6/1000                                       #g/L di O2 richiesta dal problema
conc_O2=c_O2/Chemical('O2').MW                      #mol/L di O2
O2_tot=conc_O2*Vol_s*1000                           #moli totali reattore piccolo

massaH2O=1500                                       #kg
moliH2O=1000*massaH2O/18                            #moli
moli_tot=moliH2O+O2_tot                             #reattore piccolo

x_O2s=O2_tot/moli_tot                               #reattore piccolo
J_O2=7e-4                                           #kg O2/(s m3) flusso di O2 richiesto
Jm_O2=J_O2/32                                       #kmol O2/(s m3)

C_TL=water.rholm/1000                               #kmol/m3
xO2eq_in    = yO2in*P/H 

yO2out      = lambda epsi: epsi*yO2in               #inizialmente epsilon (grado di avanzamento) è 80% e poi lo faremo cambiare
xO2eq_out   = lambda epsi: yO2out(epsi)*P/H
DeltaXml_s  = lambda epsi: ((xO2eq_in-x_O2s)-(xO2eq_out(epsi)-x_O2s))/(np.log((xO2eq_in-x_O2s)/(xO2eq_out(epsi)-x_O2s)))

#Bilancio di materia
JO2_bilmat  = lambda epsi: Q_g/32*(yO2in-yO2out(epsi)) 
#Flusso di materia
JO2_flusmat = lambda epsi: kla_s*C_TL*DeltaXml_s(epsi)*Vol_s #kmol/s

diff        = lambda epsi: JO2_bilmat(epsi)-JO2_flusmat(epsi)
[epsilon]=fsolve(diff,0.9)

if np.abs(JO2_bilmat(epsilon)-Jm_O2)>1e-3:
    print('Flussi non convergenti')


################################################### DISTRIBUTORE GAS ######################################
#diametro fori tra 3 e 6-7 mm
#Re in uscita dai fori >10000
df_s=0.003                                          #m diametro dei fori del gas sparger
A_s=0.25*np.pi*df_s**2                              #m area singolo foro
d_s=0.5*D_s                                         #m distanza distributore da impeller
L_s=np.pi*0.75*D_s                                  #m circonferenza dello sparger (diametro) così le bolle investono l'impeller

# si impone velocità gas dentro il distributore (circa 10/20 m/s)

vg_s=10                                             #m/s primo tentativo
Apassaggio_aria_s=Q_g/vg_s                          #m2 sezione di passaggio nel tubo dove passa il gas
dsparger_s=(Apassaggio_aria_s*4/np.pi)**0.5         #m diametro tubo
[dsparger_s]=np.around([dsparger_s],decimals=2)     #lo approssimiamo al diametro disponibile più vicino nei cataloghi dei tubi
                                                    #approssima di un cm 
Apassaggio_aria_s=0.25*np.pi*dsparger_s**2          #m2
vg_s=Q_g/Apassaggio_aria_s

Nfori_s=30
distfori_s=0.0297                                   # m tra un centro del foro e l'altro misurato su autocad
                                                   
vfori_s=Q_g/(Nfori_s*A_s)                           #m/s velocità uscita fori gas sparger
Re_dist_s=air.rho*vfori_s*df_s/air.mu

h2o_hot=Chemical('H2O',T=Tfin)

V_libero_s=Vr_s-Vol_s                               #m3
dP_s=h2o_hot.Psat-P*101325                          #Pascal differenza pressione liquido a 110 e P del gas
molih2o_s=dP_s*V_libero_s/(8.314*(Tfin))            #moli H2O evaporata
massah2o_s=molih2o_s*18                             #grammi H2O evaporata

Ql_s=molih2o_s*h2o_hot.Hvapm                        #Joule calore latente
Qs_s=moliH2O*h2o_hot.Cplm*(Tfin-Tin)                #Joule calore sensibile
Qr_s=Ql_s+Qs_s                                      #Joule calore totale

################################################### DIMENSIONAMENTO SEMITUBO ESTERNO ##############################
#DN 90 desterno 101.6 acciaio inox ---->tubi acciaio inox pdf
de_semitubo_s=0.1016                                #m da catalogo
spessore_s=0.00305                                  #m spessore semitubo esterno
di_semitubo_s=de_semitubo_s-2*spessore_s            #m diametro interno semitubo
semitubo_passo_s=1.13*de_semitubo_s                 #m
A_passaggio_semitubo_s=(np.pi*di_semitubo_s**2)/8   #m2

#DIMENSIONAMENTO CALOTTA
h_calotta_s=T_s/6                                                   #m
V_calotta_s=np.pi*(h_calotta_s/6)*(3*(T_s/2)**2+h_calotta_s**2)     #m3

Nspire_s=(Z_s-h_calotta_s)/semitubo_passo_s
[Nspire_s]=np.around([Nspire_s-1],decimals=0)
spessore_reattore_s=0.005                           #m spessore reattore 
Dext_reat_s=2*spessore_reattore_s+T_s       #m2
Lsemitubo_s=Nspire_s*np.pi*Dext_reat_s      #m2

#SCAMBIO TERMICO
h2o=Chemical('Water',T=((Tin+Tfin)/2))
Ascambio_s=di_semitubo_s*Lsemitubo_s                #m2
vliquidoservizio_s=1.5                              #m/s
Qservizio_s=vliquidoservizio_s*A_passaggio_semitubo_s*1000

delta_T_s=10
deq_s=(4*np.pi*di_semitubo_s**2/8)/(di_semitubo_s)
dmedio=(Dext_reat_s+Dext_reat_s+di_semitubo_s)/2

#SEZIONE TUBO
t1=140+273.15

#FISSIAMO LA T DEL FLUIDO DI SERVIZIO E CON THERMO CI TROVIAMO LA PSAT
Psat_s=Chemical('H2O',T=((t1+(Tfin+delta_T_s))/2)).Psat
ServiceFluid=Chemical('H2O',T=((t1+(Tfin+delta_T_s))/2),P=Psat_s*1.01)          #ci calcoliamo le proprietà del fluido di servizio ad una T media tra
                                                                                #quella di ingresso nella camicia e finale dell'acqua nel reattore+deltaTmin
Re_tubo_s=ServiceFluid.rho*vliquidoservizio_s*deq_s/(ServiceFluid.mu)
Pr_tubo_s=ServiceFluid.Cpl*ServiceFluid.mu/(ServiceFluid.k)
Nu_tubo_s=0.027*(Re_tubo_s**0.8)*(Pr_tubo_s**0.33)*(1+3.5*(deq_s/dmedio))
h_tubo_s=Nu_tubo_s*ServiceFluid.k/deq_s                                         #W/m2*K

Re_reattore_s=h2o.rho*N_s*D_s**2/(h2o.mu)
Pr_reattore_s=h2o.Cpl*h2o.mu/h2o.k
Nu_reattore_s=0.74*(Re_reattore_s**0.67)*(Pr_reattore_s**0.33)
h_reattore_s=Nu_reattore_s*ServiceFluid.k/T_s                                   ##W/m2*K

#COEFFICIENTE SCAMBIO GLOBALE
R=5e-4                                                                          #k*m2/W
thermalconductivity_inox=17                                                     #W/mK AISI 304 NOI ABBIAMO SCELTO 304L PERO CHE HA UNA K LEGGERMENTE PIU BASSA -->16.2

U_s=(1/h_reattore_s+1/h_tubo_s*(Dext_reat_s/T_s)+2*R+spessore_reattore_s/thermalconductivity_inox)**(-1) 
K=np.exp(U_s*Ascambio_s/(Qservizio_s*ServiceFluid.Cpl))

# print('k',h2o.k)
# print('Cp',h2o.Cpl)
# print('mu',h2o.mu)
# ServiceFluid=Chemical('H2O',T=(140+(110+delta_T_s))/2+273.15,P=280086.2260437338)
# print('k',ServiceFluid.k)
# print('Cp',ServiceFluid.Cpl)
# print('mu',ServiceFluid.mu)

def TempReactor(tau,Tr,mH2O,U,Ascambio,K):
    dTdt = [0]
    dTdt = (U*Ascambio)/(mH2O*h2o.Cpl)*((t1-Tr-(t1-Tr)/K))/(np.log((t1-Tr)/(Tr+((t1-Tr)/K)-Tr)))
    return dTdt

tf=86400                                                                       #s tempo finale vettore
npunti=tf+1
tau_span=np.linspace(0,tf,npunti)
sol=solve_ivp(TempReactor,t_span=[0.,tf],y0=[298.15],t_eval=tau_span,args=[massaH2O,U_s,Ascambio_s,K])

[Tr]=sol.y
# t2=Tr+(t1-Tr)/K
# tempofinale=np.where(np.isclose(Tr,Tfin))[0]
# t_star=tempofinale[0]

t2=np.zeros(npunti) #queste righe di codice andavano bene per il nostro caso, con le righe sopra pero
for i in range(npunti):
    Temp=sol.y[0,i]
    t2[i]=Temp+(t1-Temp)/K
    if Temp >= Tfin:
        break
    t_star=tau_span[i]


################################################## PERDITE DI CARICO ######################################
Deqfl_s=4*np.pi*(di_semitubo_s**2/8)/(np.pi*di_semitubo_s/2+di_semitubo_s)                  #m diametro equivalente fluidodinamico
Re_fl_s=ServiceFluid.rho*vliquidoservizio_s*Deqfl_s/ServiceFluid.mu
f_s=0.0791*Re_fl_s**(-0.25)                                                                 #Coeff. fanning
deltaPdistribuite_s=4*f_s*(Lsemitubo_s/Deqfl_s)*0.5*ServiceFluid.rho*vliquidoservizio_s**2  #Pa Perdite di carico

################################################### GRAFICI ###############################################
plt.figure(1)
plt.title('Reattore pilota (1.5 $m^{3}$)')

#FROUDE
curva_Froude= lambda x: 0*x+0.04
x=np.linspace(0.01,1,1001) #Na
plt.plot(x,curva_Froude(x),linestyle='--',label='Minimum dispersion speed')

#CAVITY LINE
plt.axvline(N_a_cl,linestyle='-.',label='Cavity Line')

#FLOODING LINE
curva_flooding= lambda x: (x/30)*(D_s/T_s)**(-3.5)
plt.plot(x,curva_flooding(x),color='k',linestyle=':',label='Flooding Line')

#RECIRCULATION LINE
curva_recirculation= lambda x: ((x/13)*((D_s/T_s)**(-5)))**0.5
plt.loglog(x,curva_recirculation(x),color='b',linestyle=':',label='Recirculation Line')
plt.xlabel('$N_{A}$')
plt.ylabel('Froude Number')

plt.plot(N_a,Fr_s,'o',color='r')
plt.xlim(0.01,1)
plt.legend()

plt.figure(2)
plt.grid()
plt.plot(tau_span/60,Tr-273.15,label='$T_{Reattore}$')
plt.plot(tau_span/60,t2-273.15,label='$T_{Servizio}$')
plt.ylabel('T $(C^{°})$')
plt.xlabel('Tempo (minuti)')
plt.axhline(110,linestyle='--',color='k')
plt.xticks([0,15,30,45,60,t_star/60])
plt.yticks([40,60,80,100,120,140])
plt.xlim(0,t_star/60)
plt.ylim(25,155)
plt.legend(loc='best')
plt.title('Andamento Temperature Reattore Pilota')

#####################################################################################################
############################# REATTORE INDUSTRIALE ##################################################
#####################################################################################################

############################# GEOMETRIA #############################################################
Vol_b=15                            #m3 Reattore grande

#Ipotizziamo inizialmente Z = T
T_b=((4*Vol_b)/np.pi)**(1/3)        #m diametro

#Si parte con Z=0.75*T
#H è altezza tank
#Z è altezza liquido

Z_b=T_b                                 #m altezza liquido
H_b=Z_b/0.75                            #m altezza reattore

#Turbina Rushton RT_6
D_b=1/3*T_b

W_b=D_b/5 #m larghezza pale
R_b=D_b/5 #m altezza pale

#Distanza palette da fondo rettaore
c_b=1/3*T_b

#mettiamo 4 baffle
B_b=T_b/11

Vr_b=0.25*np.pi*T_b**2*H_b #volume del reattore m3

#Imposizione di vvm portata volumentrica di gas su volume reattore m3/min/m3 
vvm = 0.4
Q_G=(vvm/60)*Vol_b          #m3/s

#bisogna calcolare Ncd e Nr, noi impondiamo diffusore a tubo forato quindi alpha = 4 
N_CD=alpha*Q_G**0.5*T_b**0.25*D_b**(-2)           #rps

N_CD=N_CD*60                        #rpm

N_R=1.5*Q_G**0.2*T_b*D_b**(-2)      #rps

N_R=N_R*60          #rpm

N_b=100             #rpm imposti reattore grande
N_b=N_b/60

V_P=np.pi*N_b*D_b   #m/s velocità periferica

########################## CHECK CONDIZIONI ###############################################
    
Fr_b=N_b**2*D_b/g

if Fr_b<0.04:
    raise ValueError('Fr industrial reactor < 0.04')

N_A=Q_G/(N_b*D_b**3)   

#CHECK FLOODING condizione N_a_fl>=N_a
N_A_fl=30*Fr_b*(D_b/T_b)**(3.5)

if N_A_fl<N_A:
    raise ValueError('N_A_fl<N_A')

# CHECK RICIRCULATION condizione N_a>=N_a_ric
N_A_ric=13*Fr_b**2*(D_b/T_b)**5

if N_A<N_A_ric:
    raise ValueError('N_A<N_A_ric')

#CHECK CAVITY LINE condizione Na>=N_a_cl
N_A_cl=0.025*(D_b/T_b)**(-0.5)

if N_A<N_A_cl:
    raise ValueError('N_A<N_A_cl')


#Se sono tutti soddisfatti sono nella zona 5 pag 64 slide (mixing slide)

############################ PORTATA MASSIMA (per non incorrere in flooding)  ###########################
Q_GMAX=N_A_fl*N_b*3600*D_b**3
Fl_PERCENT=(Q_G*3600/(Q_GMAX))*100

############################ POTENZA RICHIESTA ########################################
h2o=Chemical('Water',T=(298.15))
Re_b=(water.rho*N_b*D_b**2)/water.mu
Np_b=5.5 #vedi grafico

Pl_b=Np_b*water.rho*N_b**3*D_b**5
#vediamo se N_a è più grande di 0.037, nel nostro caso è più grande (pag. 58 slide)
Pgl_b=Pl_b*(0.62-1.85*N_A) #watt
vsg_b=4*Q_G/(np.pi*T_b**2)

kla_b=2.3*((Pgl_b/1000)/Vol_b)**0.7*(vsg_b)**0.6

#ADESSO BISOGNA PORTARE A CONVERGENZA LE DUE EQUAZIONI DEL BILANCIO DI MATERIA E DEL FLUSSO MATERIALE
O2_TOT=conc_O2*Vol_b*1000                #moli totali reattore grande

MASSAH2O=15000                                     #kg
MOLIH2O=1000*MASSAH2O/18                           #moli
MOLI_TOT=MOLIH2O+O2_TOT                            #reattore grande
x_O2b=O2_TOT/MOLI_TOT                              #reattore grande

DeltaXml_b      = lambda epsi: ((xO2eq_in-x_O2b)-(xO2eq_out(epsi)-x_O2b))/(np.log((xO2eq_in-x_O2b)/(xO2eq_out(epsi)-x_O2b)))

#Bilancio di materia
JO2_bilmat_b    = lambda epsi: Q_G/32*(yO2in-yO2out(epsi)) 
#Flusso di materia
JO2_flusmat_b   = lambda epsi: kla_b*C_TL*DeltaXml_b(epsi)*Vol_b  #kmol/s

diff_b          = lambda epsi: JO2_bilmat_b(epsi)-JO2_flusmat_b(epsi)
[epsilonb]      = fsolve(diff_b,0.9)

if np.abs(JO2_bilmat_b(epsilonb)-Jm_O2)>1e-3:
    print('Flussi non convergenti')


################################################### DISTRIBUTORE GAS ######################################
#diametro fori tra 3 e 6-7 mm
#Re in uscita dai fori >10000
df_b=0.006                  #m diametro dei fori del gas sparger
A_b=0.25*np.pi*df_b**2      #m area singolo foro
d_b=0.5*D_b                 #m distanza distributore da impeller
L_b=np.pi*0.75*D_b          #m lunghezza dello sparger così le bolle investono l'impeller

# si impone velocità gas dentro il distributore (circa 10/20 m/s)
vg_b=10                                          #m/s primo tentativo

Apassaggio_aria_b=Q_G/vg_b                       #m2 sezione di passaggio nel tubo dove passa il gas
dsparger_b=(Apassaggio_aria_b*4/np.pi)**0.5      #m 
[dsparger_b]=np.around([dsparger_b],decimals=2)  #lo approssimiamo al diametro disponibile più vicino nei cataloghi dei tubi DN100
                                                 # approssima di un cm 
Apassaggio_aria_b=0.25*np.pi*dsparger_b**2       #m2
vg_b=Q_G/Apassaggio_aria_b

Nfori_b=100
distfori_b=0.01739                               # m tra un centro del foro e l'altro misurato su autocad
                                                  
vfori_b=Q_G/(Nfori_b*A_b)
Re_dist_b=air.rho*vfori_b*df_b/air.mu

# costruttivamente poi abbiamo scelto il diametro dello sparger (vista in pianta) come i 3/4 del diametro della girante D_b

V_libero_b=Vr_b-Vol_b                               #m3
dP_b=h2o_hot.Psat-P*101325                          #Pascal differenza pressione liquido a 110 e P del gas
MOLIH2OEV=dP_b*V_libero_b/(8.314*(Tfin))            #moli H2O evaporata
MASSAH2OEV=MOLIH2OEV*18                             #grammi H2O evaporata

Ql_b=MOLIH2OEV*h2o_hot.Hvapm                        #Joule calore latente
Qs_b=MOLIH2O*h2o_hot.Cplm*(Tfin-Tin)                #Joule calore sensibile
Qr_b=Ql_b+Qs_b                                      #Joule calore totale

############################################## DIMENSIONAMENTO SEMITUBO ESTERNO ##############################
#DN 125 desterno 141.3 acciaio inox ---->tubi acciaio inox pdf
de_semitubo_b=0.1413                                 #m da catalogo
spessore_b=0.0034                                    #m spessore semitubo esterno
di_semitubo_b=de_semitubo_b-2*spessore_b             #m diametro interno semitubo
semitubo_passo_b=1.13*de_semitubo_b
A_passaggio_semitubo_b=(np.pi*di_semitubo_b**2)/8

#DIMENSIONAMENTO CALOTTA
h_calotta_b=T_b/6                                                   #m
V_calotta_b=np.pi*(h_calotta_b/6)*(3*(T_b/2)**2+h_calotta_b**2)     #m3

VolCilindro=np.pi*0.25*Z_b**2*Z_b-V_calotta_b #m3 Volume cilindro senza calotta
#assumiamo la densità dell'acqua =1000 kg/m3
hLiqCilindro= VolCilindro*1000/(np.pi*0.25*Z_b**2) #m altezza del liquido nel rettore considerando solo la parte cilindrica
hLiqReale=h_calotta_b+hLiqCilindro

Nspire_b=(Z_b-h_calotta_b)/semitubo_passo_b
[Nspire_b]=np.around([Nspire_b-1],decimals=0)
spessore_reattore_b=0.007                       #m spessore reattore 
Desterno_reattore_b=2*spessore_reattore_b+T_b   #m2
Lsemitubo_b=Nspire_b*np.pi*Desterno_reattore_b  #m2

hspira_da_fondo=0.1                             #m
h_max_spira=hspira_da_fondo+(Nspire_b+1)*(semitubo_passo_b) #altezza della spira più alta

# CHECK SE ALTEZZA LIQUIDO è MAGGIORE DELL'ALTEZZA DELLA SPIRA PIU' ALTA
if h_max_spira>hLiqCilindro:
    print('Liquido è più basso della spira più alta')

#SCAMBIO TERMICO
h2o=Chemical('Water',T=((Tfin+Tin)/2))
Ascambio_b=di_semitubo_b*Lsemitubo_b #m2
vliquidoservizio_b=1.5 #m/s
Qservizio_b=vliquidoservizio_b*A_passaggio_semitubo_b*1000 #kg/s

delta_T_b=10
deq_b=(4*np.pi*di_semitubo_b**2/8)/(di_semitubo_b)
dmedio_b=(Desterno_reattore_b+Desterno_reattore_b+di_semitubo_b)/2

#SEZIONE TUBO
Re_tubo_b=ServiceFluid.rho*vliquidoservizio_b*deq_b/(ServiceFluid.mu)
Pr_tubo_b=ServiceFluid.Cpl*ServiceFluid.mu/(ServiceFluid.k)
Nu_tubo_b=0.027*(Re_tubo_b**0.8)*(Pr_tubo_b**0.33)*(1+3.5*(deq_b/dmedio_b))
h_tubo_b=Nu_tubo_b*ServiceFluid.k/deq_b                                         #W/m2*K

Re_reattore_b=h2o.rho*N_b*D_b**2/(h2o.mu)
Pr_reattore_b=h2o.Cpl*h2o.mu/h2o.k
Nu_reattore_b=0.74*(Re_reattore_b**0.67)*(Pr_reattore_b**0.33)
h_reattore_b=Nu_reattore_b*ServiceFluid.k/T_b                                   #W/m2*K

#COEFFICIENTE SCAMBIO GLOBALE
R=5e-4                                                                          #k*m2/W
thermalconductivity_inox=17                                                     #W/mK AISI 304
U_b=(1/h_reattore_b+1/h_tubo_b*(Desterno_reattore_b/T_b)+2*R+spessore_reattore_b/thermalconductivity_inox)**(-1) 
# print('k',h2o.k)
# print('Cp',h2o.Cpl)
# print('mu',h2o.mu)
# ServiceFluid=Chemical('H2O',T=(140+(110+delta_T_s))/2+273.15,P=280086.2260437338)
# print('k',ServiceFluid.k)
# print('Cp',ServiceFluid.Cpl)
# print('mu',ServiceFluid.mu)

K_b=np.exp(U_b*Ascambio_b/(Qservizio_b*ServiceFluid.Cpl))

# def TempReactor_b(tau,Tr):
#     dTdt=[0]
#     dTdt = (U_b*Ascambio_b)/(15000*h2o.Cpl)*((t1-Tr-(t1-Tr)/K_b))/(np.log((t1-Tr)/(Tr+((t1-Tr)/K_b)-Tr)))
#     return dTdt

tau_span_b=np.linspace(0,tf,npunti)
sol_b=solve_ivp(TempReactor,t_span=[0.,tf],y0=[Tin],t_eval=tau_span_b,args=[MASSAH2O,U_b,Ascambio_b,K_b])

[Tr_b]=sol_b.y
t2_b=np.zeros(npunti)
for i in range(npunti):
    Temp_b=sol_b.y[0,i]
    t2_b[i]=Temp_b+(t1-Temp_b)/K_b
    if Temp_b >= 383.15:
        break
    t_star=tau_span_b[i]


######### PERDITE DI CARICO #############################################################################
Deqfl_b=4*np.pi*(di_semitubo_b**2/8)/(np.pi*di_semitubo_b/2+di_semitubo_b)                  #m diametro equivalente fluidodinamico
Re_fl_b=ServiceFluid.rho*vliquidoservizio_b*Deqfl_b/ServiceFluid.mu
f_b=0.0791*Re_fl_b**(-0.25)                                                                 #Coeff. Fanning
deltaPdistribuite_b=4*f_b*(Lsemitubo_b/Deqfl_b)*0.5*ServiceFluid.rho*vliquidoservizio_b**2  #Pa Perdite di carico

################################################### GRAFICI ###############################################
plt.figure(3)
plt.title('Reattore industriale (15 $m^{3}$)')
#FROUDE
curva_Froude= lambda x: 0*x+0.04
x=np.linspace(0.01,1,1001) #Na
plt.plot(x,curva_Froude(x),linestyle='--',label='Minimum dispersion speed')

#CAVITY LINE
plt.axvline(N_A_cl,linestyle='-.',label='Cavity Line')

#FLOODING LINE
curva_flooding= lambda x: (x/30)*(D_b/T_b)**(-3.5)
plt.plot(x,curva_flooding(x),color='k',linestyle=':',label='Flooding Line')

#RECIRCULATION LINE
curva_recirculation= lambda x: ((x/13)*((D_b/T_b)**(-5)))**0.5
plt.loglog(x,curva_recirculation(x),color='b',linestyle=':',label='Recirculation Line')
plt.xlabel('$N_{A}$')
plt.ylabel('Froude Number')

plt.plot(N_A,Fr_b,'o',color='r')
plt.xlim(0.01,1)
plt.legend()

# x=np.linspace(0,0.00001,1001)
# plt.figure(1)
# plt.title('Solubilità')
# plt.grid()
# plt.plot(x,curva_eq(x))
# plt.ylabel('$y_{O2}$')
# plt.xlabel('$x_{O2}$')

plt.figure(4)
plt.grid()
plt.plot(tau_span_b/60,Tr_b-273.15,label='$T_{Reattore}$')
plt.plot(tau_span_b/60,t2_b-273.15,label= '$T_{Servizio}$')
plt.ylabel('T $(C^{°})$')
plt.xlabel('Tempo (minuti)')
plt.axhline(110,linestyle='--',color='k')
plt.xticks([0,20,40,60,80,100,120,140,160,t_star/60])
plt.yticks([40,60,80,100,120,140])
plt.xlim(0,t_star/60)
plt.ylim(20,150)
plt.title('Andamento Temperature Reattore Industriale')
plt.legend()

#######################CONFRONTO REATTORI ##########################################
#CASO JO2 COSTANTE
plt.figure(5)
V1=Vr_s #piccolo
V2=Vr_b #grande

P1=Pgl_s #piccolo
P2=Pgl_b #grande
plt.loglog([V1/V1,V2/V1],[(P1/V1)/(P1/V1),(P2/V2)/(P1/V1)],label='$J_{O2}$ costante') #questo è il caso studiato
plt.grid(True, which='both', axis='x')
plt.grid(True, which='both', axis='y')
# V TIP costante
plt.plot([(V1/V1),(V2/V1)],[((P1/V1)/(P1/V1)),(V2/V1)**(-1/3)],label='$V_{TIP}$ costante') #questo è il caso V_tip costante
#tau mix costante
plt.plot([(V1/V1),(V2/V1)],[((P1/V1)/(P1/V1)),(V2/V1)**(2/3)],label='$tau_{MIX}$ costante') #questo è il caso tau_mix costante
plt.xlabel('$V_2/V_1$')
plt.ylabel('$(P_2/V_2)/(P_1/V_1)$')
plt.title('Confronto caso studio e criteri scale-up')
plt.xlim(1,10)
plt.legend()

# DIMENSIONAMENTO BOCCHELLI

#BOCCHELLO CARICA LIQUIDO
#Ipotizziamo 30 minuti per la carica
taucarica=15                                    #minuti per la carica
PortataLiqCarica_A=(15/(60*taucarica))
vtentativo=1.1                                  #m/s
Areabocchello_A=PortataLiqCarica_A/vtentativo   #m2
Dbocchello_A=(4*Areabocchello_A/(np.pi))**0.5   #m
Dbocchello_A=0.13576                            # DN125 spessore 2.77 mm Desterno=0.13576+0.00277
Abocchello_A=0.25*np.pi*Dbocchello_A**2
v_bocchello_A=PortataLiqCarica_A/Abocchello_A

#BOCCHELLO ENTRATA GAS
#Ipotizziamo che sia dello stesso diemtro del tubo dello sparger
Dbocchello_B=dsparger_b

#BOCCHELLO USCITA GAS
#Ipotizziamo che sia identico a quello di entrata, se 
Dbocchello_C=Dbocchello_B

#BOCCHELLO IN FLUIDO SERVIZIO
Dbocchello_D=di_semitubo_b
v_bocchello_D=(Qservizio_b/1000)/(np.pi*0.25*Dbocchello_D**2)
#BOCCHELLO OUT FLUIDO SERVIZIO
Dbocchello_E=di_semitubo_b

#BOCCHELLI SENSORI

#TERMOCOPPIA DN 10

#PASSO D'UOMO DN 650 660.4 con 9.52 mm di spessore

#ALBERO MOTORE
# TRIFASE 6 POLI POTENZA 18.5 KW SIGLA MOTORE:200LSA 
#DIAMETRO ALBERO 55 mm
#DIAMETRO SEZIONE ACCOPPIAMENTO 59 mm
#DIAMETRO ASTA IMPELLER DN 125 135.75 mm con spessore viene 3.4 mm
  




