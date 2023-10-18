# -*- coding: utf-8 -*-
"""
Created on Sat May 28 09:14:17 2022

@author: Andrea Milazzo
"""

import numpy as np
import  matplotlib.pyplot as plt
from scipy.optimize import fsolve
from thermo.chemical import Chemical,Mixture
plt.close('all')

#Sistema KClO3 e Na2SO4 le solubilità sono in g/100g H2O
# T_vec   = np.linspace(0,100,11)
# c_KClO3 = np.array([3.3,5,7.4,10.5,14,19.3,24.5,31.5,38.5,47.75,57])  #tabella 2.20 pag2-85 Perry 9
# T_vec1  = np.array([-1.2,10,20.5,23.5,28.5,32.4,35,40,45,50,54,59,64,69,79,90,100,102.2])
# c_Na2SO4= np.array([4,9.8,20.2,25.2,36.5,51,50.3,48.5,47.5,46.5,46.1,45.7,45.3,45,44.3,43.8,43,42.9]) #Nasini
# c_K2SO4 = np.array([7.35,9.22,11.11,12.97,14.76,16.5,18.17,19.75,21.4,22.8,24.1]) Perry 9

T_vec   = np.linspace(0,100,1001)       #°C
T_vec1  = np.linspace(-1.2,32.4,1036)   #°C
T_vec2  = np.linspace(32.4,102.2,699)   #°C

#CURVE OTTENUTE DA INTERPOLAZIONI LINEARI SU CFTOOL (MATLAB)
c_KClO3     = lambda x: (-3.531e-13*x**8 + 1.27e-10*x**7 - 1.819e-08*x**6 + 1.324e-06*x**5 - 5.174e-05*x**4 + 0.001068*x**3 - 0.007249*x**2 + 0.176*x + 3.3)/100
c_Na2SO4_1  = lambda x: (7.724e-05*x**4 - 0.003352*x**3 + 0.06742*x**2 + 0.154*x +4.083)/100 #fino a 32.4°C
c_Na2SO4_2  = lambda x: (7.031e-07*x**4 - 0.0002377*x**3 + 0.02962*x**2 - 1.672*x + 81.51)/100
c_K2SO4     = lambda x: -1.756e-06*x**3 - 4.347e-05*x**2 + 0.1897*x + 7.346

plt.figure(1)
plt.plot(T_vec,c_KClO3(T_vec)*100,label='Solubilità $KClO_3$')
plt.plot(T_vec,c_K2SO4(T_vec),color='g',label='Solubilità $K_2SO_4$')
plt.plot(T_vec1,c_Na2SO4_1(T_vec1)*100,label='Solubilità $Na_2SO_4$',color='orange')
plt.plot(T_vec2,c_Na2SO4_2(T_vec2)*100,color='orange')
plt.grid()
plt.xlim(0,100)
plt.ylim(0,60)
plt.xlabel('T (°C)')
plt.ylabel('g soluto/100 g $H_2O$')
plt.legend()

plt.figure(2)
plt.plot(T_vec,c_KClO3(T_vec)*100,label='Solubilità $KClO_3$')
plt.plot(T_vec1,c_Na2SO4_1(T_vec1)*100,label='Solubilità $Na_2SO_4$',color='orange')
plt.plot(T_vec2,c_Na2SO4_2(T_vec2)*100,color='orange')
plt.grid()
plt.xlim(0,100)
plt.ylim(0,60)
plt.xlabel('T (°C)')
plt.ylabel('g soluto/100 g $H_2O$')
plt.text(15,10,'$Na_2SO_4$')
plt.text(25,5.5,'$KClO_3$')

plt.text(71,29,'$Na_2SO_4$')
plt.text(71,22,'$KClO_3$')


plt.text(63.7,39,'$Na_2SO_4$')
plt.text(76,33,'$KClO_3$')

plt.text(88,45,'$Na_2SO_4$')
plt.text(87,55,'$KClO_3$')

#CONDIZIONI FEED
cfa=0.085   #Solfato
cfb=0.08    #Clorato
Tf=35
plt.plot([Tf],[cfa*100],'o',color='b',label='Feed')
plt.plot([Tf],[cfb*100],'o',color='b')

################################# CONDIZIONI CONCENTRATORE IN USCITA
T0=70
# c0a=0.35
c0b=0.90*c_KClO3(T0)

################################# CONDIZIONI CRISTALLIZZATORE 1
Tc1=75
c1a=0.87*c_Na2SO4_2(Tc1) #siamo sotto la curva di saturazione
c1b=c_KClO3(Tc1)

################################# CONDIZIONI CRISTALLIZZATORE 2
Tc2=96
c2a=c_Na2SO4_2(Tc2)
c2b=0.55

################################# CRISTALLIZZATORE 1
R1=0.625          #R1=c0b/c1b=Fc1/F0
c0a=R1*c1a
R1_1=1-R1       # R1_1 è rapporto evaporato e quello che entra
R1_2=c0b-R1*c1b #rapporto Pc1/F0
Pc1=400
F0=Pc1/R1_2
E1=F0*R1_1
Fc1=R1*F0

plt.plot([T0],[c0a*100],'o',color='g',label='Concentratore')
plt.plot([T0],[c0b*100],'o',color='g')
plt.plot([Tc1],[c1a*100],'o',color='r',label='Cristallizzatore 1')
plt.plot([Tc1],[c1b*100],'o',color='r')

################################# SEPARATORE 1
Ui1=0.07                        # kg H2O/kg sale umidità in uscita con KClO3
Fs1=Fc1-Pc1*Ui1*1/(1+c1a+c1b)   # kg/h corrente diretta alla sezione di cristallizzazione C2

################################# CRISTALLIZATORE 2
Fc2=Fs1*c1b/c2b
Pc2=Fs1*c1a-Fc2*c2a
E2=Fs1-Fc2

################################# SEPARATORE 2
Ui2=0.07                        # kg H2O/kg saleumidità in uscita con Na2SO4
Fs2=Fc2-Pc2*Ui2*1/(1+c2a+c2b)   # kg/h riciclo

def BilanciMat(x):
    sol=np.zeros(5)
    
    F = x[0]
    Fc = x[1]
    E0 = x[2]
    cca = x[3]
    ccb = x[4]
    
    sol[0] = F + Fs2 - Fc
    sol[1] = F*cfa + Fs2*c2a - Fc*cca
    sol[2] = F*cfb + Fs2*c2b - Fc*ccb
    sol[3] = Fc - E0 - F0
    sol[4] = Fc*cca - F0*c0a
    
    return sol

[F,Fc,E0,cca,ccb]=fsolve(BilanciMat,x0=[1000,1000,100,0.5,0.5])
F=(F0*c0a-Fs2*c2a)/(cfa)    # kg/h Feed
Fc=F+Fs2                    # kg/h portata in ingresso al concentratore

rho_feed=(Chemical('H2O',T=273.15+Tf).rhol+cfa*Chemical('Na2SO4',T=273.15+Tf).rhos+cfb*Chemical('KClO3',T=273.15+Tf).rhos)/((1+cfa+cfb))

plt.plot([Tc2],[c2a*100],'o',color='purple',label='Cristallizzatore 2')
plt.plot([Tc2],[c2b*100],'o',color='purple')
plt.legend()

################################# CALCOLO CP
Cp_f=(Chemical('H2O',T=273.15+Tf).Cpl+cfa*Chemical('Na2SO4',T=273.15+Tf).Cpl+cfb*Chemical('KClO3',T=273.15+Tf).Cpl)/(1000*(1+cfa+cfb)) #kJ/kg C feed
Cp_0=(Chemical('H2O',T=273.15+T0).Cpl+c0a*Chemical('Na2SO4',T=273.15+T0).Cpl+c0b*Chemical('KClO3',T=273.15+T0).Cpl)/(1000*(1+c0a+c0b)) #kJ/kg C uscita concentratore
Cp_1=(Chemical('H2O',T=273.15+Tc1).Cpl+c1a*Chemical('Na2SO4',T=273.15+Tc1).Cpl+c1b*Chemical('KClO3',T=273.15+Tc1).Cpl)/(1000*(1+c1a+c1b)) #kJ/kg C uscita cristallizzatore c1
Cp_2=(Chemical('H2O',T=273.15+Tc2).Cpl+c2a*Chemical('Na2SO4',T=273.15+Tc2).Cpl+c2b*Chemical('KClO3',T=273.15+Tc2).Cpl)/(1000*(1+c2a+c2b)) #kJ/kg C uscita crist. 2
Tc_prova=59 #ipotizziamo una T per calcolare Cp_c e poi calcolare la Tc reale, reiteriamo fino a convergenza
Cp_c=(Chemical('H2O',T=273.15+Tc_prova).Cpl+cca*Chemical('Na2SO4',T=273.15+Tc_prova).Cpl+ccb*Chemical('KClO3',T=273.15+Tc_prova).Cpl)/(1000*(1+cca+ccb)) #kJ/kg C

#BILANCIO ENTALPICO MIXING POINT
Tc=(F*(1+cfa+cfb)*Cp_f*Tf+Fs2*(1+c2a+c2b)*Cp_2*Tc2)/(Fc*(1+cca+ccb)*Cp_c) #Temperatura all'ingresso al concentratore che deve essere minore di T0=70°C

#BILANCIO ENTALPICO DEL CONCENTRATORE
Q_c=(E0*Chemical('H2O',T=273.15+T0).Hvap/1000+Fc*(1+cca+ccb)*Cp_c*(T0-Tc))/3600 #kW
qca = -0.28     #kcal/mol solfato sono calori di cristallizzazione (assorbono calore) da Tabella 2-72 dal “Perry's Chemical Engineers' Handbook”, edizione 9
qcb = -10.31    #kcal/mol clorato

Q_c1=(E1*Chemical('H2O',T=273.15+Tc1).Hvap/1000+F0*(1+c0a+c0b)*Cp_0*(Tc1-T0)+4186*qcb*Pc1/(Chemical('KClO3').MW))/3600      #kW
Q_c2=(E2*Chemical('H2O',T=273.15+Tc2).Hvap/1000+Fs1*(1+c1a+c1b)*Cp_1*(Tc2-Tc1)+4186*qca*Pc2/(Chemical('Na2SO4').MW))/3600   #kW

#Delta T Ebullioscopici
ioni    = 5
kb      = 0.512
b0      = 1000*(cca/Chemical('Na2SO4').MW+ccb/Chemical('KClO3').MW) #molalità
b1      = 1000*(c0a/Chemical('Na2SO4').MW+c0b/Chemical('KClO3').MW) #molalità
b2      = 1000*(c1a/Chemical('Na2SO4').MW+c1b/Chemical('KClO3').MW) #molalità

deltaT_eb0=ioni*kb*b0 #°C deltaT eb. nella sezione del concentratore
deltaT_eb1=ioni*kb*b1 #°C deltaT eb. nella sezione del Cristallizzatore C1
deltaT_eb2=ioni*kb*b2 #°C deltaT eb. nella sezione del Cristallizzatore C2

#RAPPORTI S/L
SL1=Pc1/(Fc1*(1+c1a+c1b)) # rapporto tra S/L uscenti dal cristallizzatore C1, se <0.03 allora abbiamo bisogno di un addensatore 
SL2=Pc2/(Fc2*(1+c2a+c2b)) # rapporto tra S/L uscenti dal cristallizzatore C2, se <0.03 allora abbiamo bisogno di un addensatore

ResaC1=Pc1/(F0*c0b)
ResaC2=Pc2/(Fs1*c1a)
ResaG=(Pc1+Pc2)/(F*(cfa+cfb))

#Purezze
P1=1-Ui1*c1a/(1+c1a+c1b)
P2=1-Ui2*c2b/(1+c2a+c2b)

#SPLIT1 il valore di 0.35 lo abbiamo imposto noi
split1=0.35
Fa1=Pc1/(split1*(1+c1a+c1b))    # va al separatore S/L e si miscela con Fad1 per dare origine a Fs1
Fad1=Fc1-Fa1                    # va al C2
Fsep1=Fs1-Fad1                  # corrente in ingresso al Cristallizzatore 2

#SPLIT2 il valore di 0.35 lo abbiamo imposto noi
split2=0.35
Fa2=Pc2/(split2*(1+c2a+c2b))    # va al separatore S/L
Fad2=Fc2-Fa2                    # forma corrente di riciclo insieme a Fsep2
Fsep2=Fs2-Fad2                  # corrente che esce dal separatore e insieme a Fad2 forma la corrente di riciclo

#Bilanci Essiccatori Si fissano le umidità in uscita dagli essiccatori
Uo1=0.01                        # umidità in uscita del KClO3 dall'essiccatore
Uo2=0.01                        # umidità in uscita del KClO3 dall'essiccatore
W1=Pc1*(Ui1-Uo1)                # portata di acqua rimossa dal sale 
W2=Pc2*(Ui2-Uo2)                # portata di acqua rimossa dal sale

#Analisi concentratore
# Il nostro delta T eb non è trascurabile ma invece di alzare la T scendo di P
# Valutiamo corrente di riciclo M 
# Prevedo riciclo molto grande e impongo delta T è molto piccolo, quindi posso stimare che le cP, conc e T
# siano approssimabili a quelle della F0

deltaT_r0   = 5                                             #°C
R0          = (Q_c*3600)/(Cp_0*deltaT_r0*(1+c0a+c0b))       #kg/h riciclo
FR0         = R0-Fc*(1+cca+ccb)/(1+c0a+c0b)                 # uscita dal concentratore che si mixa con Fc
Psat_R0     = Chemical('H2O',T=273.15+T0-deltaT_eb0).Psat   #spingiamo di piu il vuoto in modo che possa bollire sempre a T0
P_R0        = Psat_R0*1.4                                   #P in uscita dallo scambiatore

rho_R0=(Chemical('H2O',T=273.15+T0+deltaT_r0).rhol+c0a*Chemical('Na2SO4',T=273.15+T0+deltaT_r0).rhos+c0b*Chemical('KClO3',T=273.15+T0+deltaT_r0).rhos)/((1+c0a+c0b)) #kg/m3
hmax_R0=(P_R0-Psat_R0)/(rho_R0*9.81)                        # m stima dell'altezza dallo scambiatore all'imbocco del tank

#TERMOCOMPRESSIONE
#Siamo sottovuoto quindi utilizzeremo la termocompressione invece che la compressione meccanica(che si usa a P superiori a Patm)
Tv3_C0=T0+20                                                    # T in uscita dall'eiettore ed in entrata nello scambiatore
Tv1_C0=180                                                      # T vapore motore
Pv1_C0=Chemical('H2O',T=273.15+Tv1_C0).Psat                     # Pa Pressione del vapore motore
V3_C0=(Q_c*3600)/(Chemical('H2O',T=273.15+Tv3_C0).Hvap/1000)    # kg/h di vapore in uscita dall'eiettore
Pv3_C0=Chemical('H2O',T=273.15+Tv3_C0).Psat                     # Pa Pressione del vapore in uscita dall'eiettore ed in entrata nello scambiatore

H1_C0=2800 #kj/kg vapore motore
H2_C0=2630 #kj/kg vapore aspirato
H3_C0=2670 #kj/kg vapore in uscita dall'eiettore

deltaH3_C0=H3_C0-H2_C0
deltaH1_C0=H1_C0-H2_C0

#Fissiamo eta: deve rimanere compreso tra 0.63-0.78
eta_C0=0.7
Reff_C0=((deltaH3_C0/(eta_C0*deltaH1_C0))**(0.5))/(1-(deltaH3_C0/(eta_C0*deltaH1_C0))**(0.5)) # Rapporto tra i vapori in ingresso (può arrivare fino a 2)
if Reff_C0>2:
    raise ValueError('Il rapporto tra i vapori in ingresso al termocompressore nella sezione di concentrazione è maggiore di 2!!')

V2_C0=V3_C0/(1+Reff_C0)     #kg/h vapore trascinato 
V1_C0=Reff_C0*V2_C0         #kg/h vapore motore
Vsplit_C0=E0-V2_C0          #kg/h di vapore che va al condensatore e poi alla raccolta delle condense

######################### DIMENSIONAMENTO DELLO SCAMBIATORE DI CALORE (HE1) ############################
U_HE1=2000                                                  # W/m2 K 637 Coulson 6
TR_0=T0+deltaT_r0                                           # °C gli diamo 3°C in più e quando flasha si prende questo deltaT 
deltaT_ML_HE1=(T0-TR_0)/np.log((Tv3_C0-TR_0)/(Tv3_C0-T0))   
Ft_HE1=1
As_HE1=(Q_c*1000)/(U_HE1*deltaT_ML_HE1*Ft_HE1)              # m2 Area di scambio 

# TABELLA TECNO STEEL PER I TUBI
Lt_HE1=4                                            # m
de_tubi_HE1=48.3/1000                               # m     DN40
s_tubi_HE1=2.77/1000                                # m
di_tubi_HE1=de_tubi_HE1-2*s_tubi_HE1                # m
Nt_HE1=As_HE1/(np.pi*de_tubi_HE1*Lt_HE1)            # calcolo numero dei tubi da approssimare
Nt_HE1=np.around(Nt_HE1)                            # numero di tubi reale
Np_HE1=2                                            # numero di passaggi dei tubi
R0_vol=(R0*(1+c0a+c0b))/rho_R0                      # m3/h
vtubi_HE1=(Np_HE1*R0_vol/(0.25*np.pi*di_tubi_HE1**2*Nt_HE1))/3600 #m/s velocità all'interno dei tubi

Pt_HE1=1.25*de_tubi_HE1                             # m passo tra i tubi
K1_HE1=0.249                                        # parametri della formula pag.648 Coulson 6
n_HE1=2.207                                         # parametri della formula pag.648 Coulson 6
Dshell_HE1=(de_tubi_HE1)*(Nt_HE1/K1_HE1)**(1/n_HE1) # m diametro del fascio di tubi
forma_HE1=Lt_HE1/(Dshell_HE1*1.05)                  # fattore di forma dello scambiatore di calore 6-10
if forma_HE1>10:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di concentrazione è maggiore di 10!!')
elif forma_HE1<6:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di concentrazione è minore di 6!!')

#CALCOLO DELLE PERDITE DI CARICO
mu0a=1.38/1000                      # (approssimazioni molto forti) Pa*s dati utili corso e-learning questi dati sono a 25°C
mu0b=0.89/1000                      # (approssimazioni molto forti) Pa*s dati utili corso e-learning questi dati sono a 25°C
mu0=(mu0a*c0a+mu0b*c0b)/(c0a+c0b)   # viscosità della soluzione nella sezione di concentrazione

epsi=0.00005                                                                #rugosità 
Re_HE1=(rho_R0*vtubi_HE1*di_tubi_HE1)/mu0
f_HE1=0.0791*Re_HE1**(-0.25)                                                # coefficiente di Fanning
deltaP_dist_HE1=4*f_HE1*Np_HE1*(Lt_HE1/di_tubi_HE1)*0.5*rho_R0*vtubi_HE1**2 # Pa perdite di carico distribuite
deltaP_conc_HE1=4*Np_HE1*0.5*rho_R0*vtubi_HE1**2                            # Pa perdite di carico concentrate
vbocc_HE1=1                                                                 # m/s velocità bocchelli sia ingresso che uscita 
deltaP_bocc_HE1=rho_R0*vbocc_HE1**2                                         # Pa perdite di carico bocchelli
deltaP_tot_HE1=deltaP_dist_HE1+deltaP_conc_HE1+deltaP_bocc_HE1              # Pa perdite di carico totali

########################## DIMENSIONAMENTO DELLA CAMERA DI FLASH (serbatoio cilindrico a fondi bombati)
K_F0=0.02                                           # m/s noi abbiamo ipotizzato di non utilizzare il demister tra 0.01 e 0.05 m/s
rhov_F0=Chemical('H2O',T=T0+273.15,P=Psat_R0).rhog  # kg/m3 vapore 
uv_F0=K_F0*((rho_R0-rhov_F0)/(rhov_F0))**0.5        # m/s
A_F0=(E0/(rhov_F0*3600))/(uv_F0)                    # m2
D_F0=(4*A_F0/(np.pi))**0.5                          # m diametro tank

l_F0=0.75*D_F0                                      #compreso tra 0.5D e 1D battente liquido nel flash
h_F0=2.5*D_F0                                       #compreso tra 2D e 3D h altezza libera
H_F0=l_F0+h_F0                                      #m altezza del tank
V_F0=A_F0*H_F0                                      #m3 volume del flash tank

################################################# DIMENSIONAMENTO BOCCHELLI ###########################################################
################################################# bocc_0A USCITA E0
#ipotizziamo una velocità
u_0A=15 #m/s
d_0A=((4/np.pi)*(E0/(rhov_F0*3600))/(u_0A))**0.5 #m
#Stiamo scegliendo un DN450 con diam esterno = 457.2 mm e spessore 4.19mm
d_0A=0.44882 #m
u_0A=(E0/(rhov_F0*3600))/(np.pi*0.25*d_0A**2) #m/s

################################################# bocc_0B INGRESSO FC
#ipotizziamo una velocità
u_0B=1 #m/s
d_0B=((4/np.pi)*((Fc*(1+cca+ccb))/(rho_R0*3600))/(u_0B))**0.5 #m
# Stiamo scegliendo un DN50 con diam esterno = 60.3 mm e spessore 2.77mm
d_0B=0.05476 #m
u_0B=(Fc*(1+cca+ccb)/(rho_R0*3600))/(np.pi*0.25*d_0B**2) #m/s

################################################# bocc_0C USCITA F0
#ipotizziamo una velocità
u_0C=1 #m/s
d_0C=((4/np.pi)*((F0*(1+c0a+c0b))/(rho_R0*3600))/(u_0C))**0.5 #m
#Stiamo scegliendo un DN40 con diam esterno = 48.3 mm e spessore 1.65mm
d_0C=0.045 #m
u_0C=(F0*(1+c0a+c0b)/(rho_R0*3600))/(np.pi*0.25*d_0C**2) #m/s

################################################# TUBAZIONE TRATTO 1-2
#Stiamo analizzando il tratto 1-2 (quello in uscita dal flash)
u0_12=1.0 #m/s
d0_12=((4/np.pi)*((FR0*(1+c0a+c0b))/(rho_R0*3600))/(u0_12))**0.5 #m
#Stiamo scegliendo un DN300 con diam esterno = 323.8 mm e spessore 4.57mm
d0_12=0.31466 #m
u0_12=(FR0*(1+c0a+c0b)/(rho_R0*3600))/(np.pi*0.25*d0_12**2) #m/s

################################################# TUBAZIONE TRATTO 3-4
#Stiamo analizzando il tratto 3-4 (quello dopo il nodo prima della pompa)
u0_34=1.0 #m/s
d0_34=((4/np.pi)*((R0*(1+c0a+c0b))/(rho_R0*3600))/(u0_34))**0.5 #m
#Stiamo scegliendo un DN300 con diam esterno = 323.8 mm e spessore 4.57mm
d0_34=0.31466 #m
u0_34=(R0*(1+c0a+c0b)/(rho_R0*3600))/(np.pi*0.25*d0_34**2) #m/s

################################################# TUBAZIONE TRATTO 5-6
#Stiamo analizzando il tratto 5-6 (quello dopo il nodo dopo della pompa fino all'ingresso nella camera di flash)
u0_56=3.0 #m/s
d0_56=((4/np.pi)*((R0*(1+c0a+c0b))/(rho_R0*3600))/(u0_56))**0.5 #m
# Stiamo scegliendo un DN125 con diam esterno = 141.3 mm e spessore 2.77 mm
d0_56=0.13576 #m
u0_56=(R0*(1+c0a+c0b)/(rho_R0*3600))/(np.pi*0.25*d0_56**2) #m/s

################################################# POMPA DI RICIRCOLO G1
Qv_G1=R0*(1+c0a+c0b)/(rho_R0)                   # m3/h
zab_G1=2                                        # m distanza tra pelo libero e entrata nel flash tank
z1_G1=3.0                                       # m altezza dal fondo del flash fino imbocco Tee dove entra Fc
z2_G1=1.5                                       # m da terra fino a tee
h0_G1=np.round(hmax_R0-0.1,1)                   # m approssimazione in difetto di hmax_R0
l0_G1=np.round(l_F0,1)                          # m
z3_G1=zab_G1+l0_G1+z1_G1+z2_G1-h0_G1-Lt_HE1     # m
l2_G1=3.5                                       # m
l1_G1=2.5+0.2                                   # m

#Coefficienti perdite di carico concentrate (da slide Pumping of gases and liquids pag. 13)
K_tee=1.8
K_90=0.9
K_in=0.5
K_exit=1

################ TRATTO IN ASPIRAZIONE #############
#TRATTO 12 (tratto in verticale dal fondo del cristallizzatore fino al tee con Fc)
Re_12=rho_R0*u0_12*d0_12/(mu0)
f_12=0.0037                                                     # dal grafico pag.6 'pumping of gases and liquids' slide prof.
deltaPd_12=4*f_12*(z1_G1/d0_12)*0.5*rho_R0*u0_12**2             # Pa perdite di carico distribuite
deltaPc_12=0.5*K_in*rho_R0*u0_12**2                             # Pa perdite di carico concentrate

#TRATTO 34 (tratto dal tee fino alla pompa G1)
Re_34=rho_R0*u0_34*d0_34/(mu0)
f_34=0.0037                                                     # dal grafico pag.6 'pumping of gases and liquids' slide prof.
deltaPd_34=4*f_34*((z2_G1+l1_G1)/d0_34)*0.5*rho_R0*u0_34**2     # Pa perdite di carico distribuite
deltaPc_34=0.5*(K_90+K_tee)*rho_R0*u0_12**2                     # Pa perdite di carico concentrate
deltaPasp_G1=deltaPd_12+deltaPd_34+deltaPc_12+deltaPc_34        # Pa Perdite di carico in aspirazione

################ TRATTO IN MANDATA #############
#TRATTO 5-6 (tratto dalla pompa fino all'ingresso nel concentratore)
Re_56=rho_R0*u0_56*d0_56/(mu0)
f_56=0.0038                                                                 # dal grafico pag.6 'pumping of gases and liquids' slide prof.
deltaPd_56=4*f_56*((z3_G1+h0_G1+l2_G1+Lt_HE1)/d0_56)*0.5*rho_R0*u0_56**2    # Pa perdite di carico distribuite
deltaPc_56=0.5*(3*K_90+K_exit)*rho_R0*u0_56**2                              # Pa perdite di carico concentrate

########################################## POMPA G1
Q_G1=R0*(1+c0a+c0b)/rho_R0                                                  # m3/h portata di riciclo al concentratore
deltaPmand_G1=deltaPd_56+deltaPc_56+deltaP_tot_HE1                          # Pa perdite di carico nella sezione di mandata
deltaPf_G1=deltaPasp_G1+deltaPmand_G1                                       # Pa perdite di carico totali su tutta la linea
Hp_G1=(deltaPf_G1+0.5*rho_R0*u0_56**2+rho_R0*9.81*zab_G1)/(rho_R0*9.81)     # m Prevalenza richiesta dal sistema, in questo caso è anche la prevalenza della pompa, visto che è stima teorica

# POMPA SCELTA TCT 50Hz 4poli 100-200 1450rpm
NPSHr_G1=2                                                                                  # m da grafico pompa
NPSHd_G1=(0.5*rho_R0*u0_34**2+rho_R0*9.81*(l0_G1+z1_G1+z2_G1)-deltaPasp_G1)/(rho_R0*9.81)   # m di fluido che stiamo pompando NPSH disponibile

#Si legge 3.8 kW con un phi di 189 (teorico)
eta_G1=0.38                     # efficienza
C_G1=rho_R0/1000                # fattore correttivo 1000 kg/m3 è la densità dell'acqua
PotN_G1=7.5                       # kW Potenza nominale
PotR_G1=PotN_G1*C_G1            # kW Potenza con fattore di correzione
PotRic_G1=PotR_G1/eta_G1        # kW Potenza richiesta effettivamente dalla pompa

Pasp_G1=NPSHd_G1*rho_R0*9.81 - 0.5*rho_R0*u0_34**2 + Psat_R0  # Pressione assoluta in aspirazione

########################################## POMPA G2 (di vuoto)
air_G2=Mixture(['air'],T=273.15+T0,P=Psat_R0)
Vol_sis_G2=V_F0*1.3                                         # m3
Vol_sisfeet_G2=Vol_sis_G2/(0.02831685)                      # ft3 Vol_sis_G2 riportato in ft3
P_sis_G2=(Psat_R0/101325)*760                               # mmHg considerando Psat_R0 in mmHg
Qasp_inc_G2_pound=7.6                                       # lb/h di incondensabili che entrano nel nostro sistema che ci siamo trovati nel nostro grafico vacuum equipment slide pag 23
Qasp_inc_G2=0.45359243*Qasp_inc_G2_pound/air_G2.rho         # m3/h di incondensabili che entrano nel nostro sistema che ci siamo trovati nel nostro grafico

#SCELGO UNA POMPA SUL GRAFICO E VEDO LA PORTATA CHE NECESSITA LA POMPA Pompa a due stadi TRH 32-20 2900 rpm
phiT_G2=0.87                    # coefficiente correttivo che bisogna applicare dato che il nostro fluido di servizio non è a 15°C ma è a 30°C
Qasp_n_G2=21                    # m3/h nominale trovato da grafico campi funzionamento
Qeff_G2=Qasp_n_G2*phiT_G2       # m3/h effettivi
Liq_es_G2=5.7                   # l/min acqua di servizio per il funzionamento dell'anello liquido
                                # la pompa funziona a 2900 giri
Pot_G2=0.8                      # kW potenza nominale assorbita dalla pompa
c_G2=air_G2.rho/Mixture(['air'],T=293.15,P=Psat_R0).rho
PotR_G2=Pot_G2*c_G2             #kW potenza reale richiesta dalla pompa

# #calcoliamo il tempo per mettere sottovuoto l'apparecchiatura
P1_G2=Psat_R0/(101325)                          # atm Pressione che la pompa deve raggiungere
P2_G2=1                                         # atm
t_G2=Vol_sis_G2/Qeff_G2*60*np.log(P2_G2/P1_G2)  # min tempo necessario a raggiungere il vuoto

# Valutazione Torre Raffreddamento E9 TMA-EU 08-55 Potenza=4kW DECSA
Qcond_E9=Vsplit_C0*Chemical('H2O',T=273.15+T0).Hvap/(3600*1000) # kW
deltaT_E9=20                                                    # °C
W_E9=Qcond_E9*3600/(deltaT_E9*(Chemical('H2O',T=273.15+T0).Cpl/1000)) #kg/h di acqua che entra nella torre di raffreddamento e passa da 50 a 30°C

########################## SEZIONE DI CRISTALLIZZAZIONE C1 ########################
deltaT_r1   = 4                                                                 # °C
R1          = (Q_c1*3600)/(Cp_1*deltaT_r1*(1+c1a+c1b))                          # kg/h riciclo
FR1         = R1-F0*(1+c0a+c0b)/(1+c1a+c1b)                                     # kg/h di soluzione in uscita dal C1
Psat_R1     = Chemical('H2O',T=273.15+Tc1-deltaT_eb1).Psat                      # Pa spingiamo di piu il vuoto in modo che possa bollire sempre a T1
P_R1        = Psat_R1*1.4                                                       # Pa pressione in uscita dallo scambiatore

rho_R1=(Chemical('H2O',T=273.15+Tc1+deltaT_r1).rhol+c1a*Chemical('Na2SO4',T=273.15+Tc1+deltaT_r1).rhos+c1b*Chemical('KClO3',T=273.15+Tc1+deltaT_r1).rhos)/((1+c1a+c1b)) #kg/m3
hmax_R1=(P_R1-Psat_R1)/(rho_R1*9.81)                                            # m stima dell'altezza dallo scambiatore all'imbocco del tank

#TERMOCOMPRESSIONE               
#Siamo sottovuoto quindi utilizzeremo la termocompressione invece che la compressione meccanica(che si usa a P superiori a Patm)
Tv3_C1=Tc1+20                                                   # T in uscita dall'eiettore ed in entrata nello scambiatore
Tv1_C1=180                                                      # T vapore motore
Pv1_C1=Chemical('H2O',T=273.15+Tv1_C1).Psat                     # Pa Pressione del vapore motore
V3_C1=(Q_c1*3600)/(Chemical('H2O',T=273.15+Tv3_C1).Hvap/1000)   # kg/h di vapore in uscita dall'eiettore
Pv3_C1=Chemical('H2O',T=273.15+Tv3_C1).Psat                     # Pa Pressione del vapore in uscita dall'eiettore ed in entrata nello scambiatore

H1_C1=2800 # kj/kg vapore motore
H2_C1=2640 # kj/kg vapore aspirato
H3_C1=2680 # kj/kg vapore uscita eiettore

deltaH3_C1=H3_C1-H2_C1
deltaH1_C1=H1_C1-H2_C1

#Fissiamo eta deve rimanere compreso tra 0.63-0.78
eta_C1=0.7
Reff_C1=((deltaH3_C1/(eta_C1*deltaH1_C1))**(0.5))/(1-(deltaH3_C1/(eta_C1*deltaH1_C1))**(0.5)) # Rapporto tra i vapori in ingresso

V2_C1=V3_C1/(1+Reff_C1)     # kg/h vapore trascinato
V1_C1=Reff_C1*V2_C1         # kg/h vapore motore
Vsplit_C1=E1-V2_C1          # kg/h di vapore che va al condensatore e poi alla raccolta delle condense

# ######################### DIMENSIONAMENTO DELLO SCAMBIATORE DI CALORE (HE2) ############################
U_HE2=1800                                                      # W/m2 K 637 Coulson 6
TR_1=Tc1+deltaT_r1                                              # °C
deltaT_ML_HE2=(Tc1-TR_1)/np.log((Tv3_C1-TR_1)/(Tv3_C1-Tc1))     # °C
Ft_HE2=1
As_HE2=(Q_c1*1000)/(U_HE2*deltaT_ML_HE2*Ft_HE2)                 #m2 Area di scambio

# TABELLA TECNO STEEL PER I TUBI
Lt_HE2=5                                                            # m
de_tubi_HE2=60.3/1000                                               # m     DN50
s_tubi_HE2=1.65/1000                                                # m
di_tubi_HE2=de_tubi_HE2-2*s_tubi_HE2                                # m
Nt_HE2=As_HE2/(np.pi*de_tubi_HE2*Lt_HE2)
Nt_HE2=np.around(Nt_HE2)                                            # numero di tubi reale 
Np_HE2=2                                                            # numero di passaggi
R1_vol=(R1*(1+c1a+c1b))/rho_R1                                      # m3/h
vtubi_HE2=(Np_HE2*R1_vol/(0.25*np.pi*di_tubi_HE2**2*Nt_HE2))/3600   # m/s
Pt_HE2=1.25*de_tubi_HE2                                             # m passo tra i tubi
K1_HE2=0.249                                                        # parametri della formula pag. 649 Coulson 
n_HE2=2.207                                                         # parametri della formula pag. 649 Coulson
Dshell_HE2=(de_tubi_HE2)*(Nt_HE2/K1_HE2)**(1/n_HE2)                 # m diametro del fascio tubiero
forma_HE2=Lt_HE2/(Dshell_HE2*1.05)

if forma_HE2>10:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di C1 è maggiore di 10!!')
elif forma_HE2<6:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di C1 è minore di 6!!')

# CALCOLO DELLE PERDITE DI CARICO
mu1=(mu0a*c1a+mu0b*c1b)/(c1a+c1b)                                               # viscosità della soluzione nella sezione di cristallizzazione 1
Re_HE2=(rho_R1*vtubi_HE2*di_tubi_HE2)/mu1                                       # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
f_HE2=0.0052                                                                    # pag 6 of Pumping of Gases and Liquids
deltaP_dist_HE2=4*f_HE2*Np_HE2*(Lt_HE2/di_tubi_HE2)*0.5*rho_R1*vtubi_HE2**2     # Pa 
deltaP_conc_HE2=4*Np_HE2*0.5*rho_R1*vtubi_HE2**2                                # Pa
vbocc_HE2=1                                                                     # m/s v bocchelli in e out
deltaP_bocc_HE2=rho_R1*vbocc_HE2**2                                             # Pa perdite di carico indotte dai bocchelli
deltaP_tot_HE2=deltaP_dist_HE2+deltaP_conc_HE2+deltaP_bocc_HE2                  # Pa perdite di carico totali

# Valutazione Torre Raffreddamento E10 TMA-EU 08-109 Potenza = 11 kW
Qcond_E10=Vsplit_C1*Chemical('H2O',T=273.15+Tc1).Hvap/(3600*1000)                # kW
deltaT_E10=20                                                                    # °C
W_E10=Qcond_E10*3600/(deltaT_E10*(Chemical('H2O',T=273.15+Tc1).Cpl/1000))          # kg/h di acqua che entra nella torre di raffreddamento e passa da 50 a 30°C

############################# SEZIONE DI CRISTALLIZZAZIONE C2 ########################
deltaT_r2   = 3                                                                 # °C
R2          = (Q_c2*3600)/(Cp_2*deltaT_r2*(1+c2a+c2b))                          # kg/h riciclo
FR2         = R2-Fs1*(1+c1a+c1b)/(1+c2a+c2b)                                    # kg/h di soluzione in uscita dal C1
Psat_R2     = Chemical('H2O',T=273.15+Tc2-deltaT_eb2).Psat                      # Pa spingiamo di piu il vuoto in modo che possa bollire sempre a T2                    
P_R2        = Psat_R2*1.4                                                       # Pa Pressione uscente dallo scambiatore

rho_R2      = (Chemical('H2O',T=273.15+Tc2+deltaT_r2).rhol+c2a*Chemical('Na2SO4',T=273.15+Tc2+deltaT_r2).rhos+c2b*Chemical('KClO3',T=273.15+Tc2+deltaT_r2).rhos)/((1+c2a+c2b)) #kg/m3
hmax_R2     = (P_R2-Psat_R2)/(rho_R2*9.81)                                      # m stima dell'altezza dallo scambiatore all'imbocco del tank
Tvap_HE3    = 120                                                               # vapore bassa P per il cristallizzatore C2
                                                             
######################### DIMENSIONAMENTO DELLO SCAMBIATORE DI CALORE (HE3) ############################
U_HE3= 2000                                                                     # W/m2 K 637 Coulson 6
W_HE3=Q_c2*3600/(Chemical('H2O',T=180+273.15,P=1013250).Hvap/1000)              # kg/h di vapore (fluido di servizio) entrante nello scambiatore
W_HE3vol=W_HE3/(Chemical('H2O',T=182+273.15,P=1013250).rhog)                    # m3/h abbiamo fatto 180+2 per avere densità vapore
TR_2=Tc2+deltaT_r2                                                              # °C
deltaT_ML_HE3=(Tc2-TR_2)/np.log((Tvap_HE3-TR_2)/(Tvap_HE3-Tc2))                 # °C
Ft_HE3=1
As_HE3=(Q_c2*1000)/(U_HE3*deltaT_ML_HE3*Ft_HE3)                                 # m2

# # # TABELLA TECNO STEEL PER I TUBI
Lt_HE3=5                                                                        # m
de_tubi_HE3=88.9/1000                                                           # m DN80
s_tubi_HE3=3.05/1000                                                            # m
di_tubi_HE3=de_tubi_HE3-2*s_tubi_HE3                                            # m
Nt_HE3=As_HE3/(np.pi*de_tubi_HE3*Lt_HE3)                                        # numero di tubi da approssimare
Nt_HE3=np.around(Nt_HE3)                                                        # numero di reale
Np_HE3=1                                                                        # numero di passaggi
R2_vol=(R2*(1+c2a+c2b))/rho_R2                                                  # m3/h di riciclo 
vtubi_HE3=(Np_HE3*R1_vol/(0.25*np.pi*di_tubi_HE3**2*Nt_HE3))/3600               # m/s
Pt_HE3=1.25*de_tubi_HE3                                                         # m passo tra i tubi
K1_HE3=0.319                                                                    # parametri della formula pag. 649 Coulson 
n_HE3=2.142                                                                     # parametri della formula pag. 649 Coulson 
Dshell_HE3=(de_tubi_HE3)*(Nt_HE3/K1_HE3)**(1/n_HE3)                             # m diametro del fascio tubiero
forma_HE3=Lt_HE3/(Dshell_HE3*1.05) 

if forma_HE3>10:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di C2 è maggiore di 10!!')
elif forma_HE3<6:
    raise ValueError('Il fattore di forma dello scambiatore nella sezione di C2 è minore di 6!!')             

#CALCOLO DELLE PERDITE DI CARICO

mu2=(mu0a*c2a+mu0b*c2b)/(c2a+c2b)                                               # viscosità della soluzione nella sezione di cristallizzazione 2
Re0_HE3=(rho_R2*vtubi_HE3*di_tubi_HE3)/mu2                                      # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
f_HE3=0.0048                                                                    # pag 6 of Pumping of Gases and Liquids
deltaP_dist_HE3=4*f_HE3*Np_HE3*(Lt_HE3/di_tubi_HE3)*0.5*rho_R2*vtubi_HE3**2     # Pa perdite di carico distribuite
deltaP_conc_HE3=4*Np_HE3*0.5*rho_R2*vtubi_HE3**2                                # Pa perdite di carico concentrazione
vbocc_HE3=1                                                                     # m/s v bocchelli in e out
deltaP_bocc_HE3=rho_R2*vbocc_HE3**2                                             # Pa perdite di carico bocchelli
deltaP_tot_HE3=deltaP_dist_HE3+deltaP_conc_HE3+deltaP_bocc_HE3                  # Pa perdite di carico totali

# Valutazione Torre Raffreddamento E11 TMA EU 08-109 Potenza=11 kW
Qcond_E11=E2*Chemical('H2O',T=273.15+Tc2).Hvap/(3600*1000)                       # kW
deltaT_E11=20                                                                    # °C
W_E11=Qcond_E11*3600/(deltaT_E11*(Chemical('H2O',T=273.15+Tc2).Cpl/1000))        # kg/h di acqua che entra nella torre di raffreddamento e passa da 50 a 30°C

#SCARICATORE DI CONDENSE SC-HE3
PA_SC=Chemical('H2O',T=273.15+Tvap_HE3).Psat                                    # a monte dell'SC
PC_SC=101325                                                                    # dentro il serbatoio di raccolta
vb_SC=1                                                                         # m/s velocità all'interno della tubazione
# la velocità pelo libero acqua = 0
h_SC=3                                                                          # m
L_SC=5                                                                          # m
D_SC=((4*W_HE3/(3600*Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol))/(np.pi*vb_SC))**(1/2)                          # m diametro dello scaricatore di condense
Re_SC=Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol*vb_SC*D_SC/(Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).mul)   # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
f_SC=0.0065                                                                                         # pumping of gases and liquids pag 6
deltaP_dist_SC=4*f_SC*(L_SC/D_SC)*0.5*Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol*vb_SC**2    # Pa perdite di carico distribuite scaricatore di condense
deltaP_conc_SC=(3*K_90+K_exit)*0.5*Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol*vb_SC**2       # Pa Perdite di carico concentrate scaricatore di condense
deltaP_tot_SC=deltaP_dist_SC+deltaP_conc_SC                                                         # Pa perdite di carico totali scaricatore di condense
PB_SC=deltaP_tot_SC+Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol*9.81*h_SC-0.5*Chemical('H2O',T=273.15+Tvap_HE3,P=3*101325).rhol*vb_SC**2+PC_SC
Pdiff_SC=PA_SC-PB_SC
 
# VALVOROBICA FT14/43 diametro DN50 Pdiff=0.7 bar scarica 2150 kg/h di condensa

########################### DIMENSIONAMENTO VALVOLA (siamo prima del concentratore)
# Abbiamo DN50 della tubazione in entrata ed in uscita
rho_C_V = (Chemical('H2O',T=273.15+Tc).rhol+cca*Chemical('Na2SO4',T=273.15+Tc).rhos+ccb*Chemical('KClO3',T=273.15+Tc).rhos)/((1+cca+ccb)) #kg/m3
Q_Vmax  = 2*(Fc*(1+cca+ccb)/rho_C_V)                                # m3/h doppio della portata di Fc
deltaP_valv = 0.3*deltaPf_G1/(1000*0.7)                             # kPa 
Cv_V    = 1.17*Q_Vmax*((rho_C_V/1000)/(deltaP_valv/100))**(0.5)     # Calcolato secondo criterio Conflow, Pdf Control valve sizing per formula, Pdf 5000l per valvola
#Nella formula del Cv la rho va in kg/dm3 e il delta P in kg/cm2
Pv_V    = Chemical('H2O',T=Tc).Psat                                 # Pa pressione di saturazione 
Pc_V    = Chemical('H2O').Pc                                        # Pa Critical Pressure 
Ff_V    = 0.96-0.28*(Pv_V/Pc_V)**(0.5)                              # Fattore pressione critica del liquido
Fl_V    = 0.9                                                       # ((p1-p2)/(p1-pvc))**0.5 ovc pressure at the vena contracta
P1_V    = 101325                                                    # Pa P atm
deltaPmax_V= Fl_V**2*(P1_V-Ff_V*Pv_V)/1000                          # deltaPmax > di deltaPvalv

########################################### DIMENSIONAMENTO POMPA PAPPA SALINA G5 (C1)
#STIMA VOLUMI CRISTALLIZZATORI
# in C1 entra F0: stimiamo un tempo di residenza tra 1 e 6 hr
# Considerando la grandezza della particelle 0.2 micron, si ipotizza un tempo di residenza medio di 1hr 
#(vedi pag.25 Industrial Crystallization)

############################ CRISTALLIZZATORE KClO3
tau_C1=1                                # hr tempo permanenza all'interno del reattore
Vl_C1=(F0*(1+c0a+c0b)/rho_R0)*tau_C1    # m3 volume di soluzione all'interno del Cristallizzatore
                                        # Facciamo l'ipotesi che sia pieno al 40%
Vc_C1=Vl_C1/0.4                         # m3 Volume del Cristallizzatore

#ipotizziamo un diametro Dc1
D_C1=2.0                                # m diametro del cristallizzatore
h_C1=Vc_C1/(0.25*np.pi*D_C1**2)         # m altezza del cristallizzatore
hL_C1=0.4*h_C1                          # m altezza di liquido all'interno del cristallizzatore

#A è in fondo alla tubazione verticale che scende dal cristallizzatore
#B è all'ingresso dell'addensatore
z1_G5=2.5       # m
l1_G5=2.5       # m
z2_G5=4         # m
l2_G5=1.5       # m
PB_G5=101325    # Pa
PA_G5=Psat_R1+rho_R1*9.81*z1_G5+rho_R1*9.81*hL_C1 # Pa Pressione nel punto A, dove abbiamo considerato Psat + i battenti di liquido del C1 + battente z1

##################### PERDITE DI CARICO
#ipotizziamo velocità
u_asp_G5=0.6                                                                    # m/s tentativo
D_asp_G5=((4/np.pi)*((Fc1*(1+c1a+c1b)+Pc1)/(rho_R1))/(u_asp_G5*3600))**(0.5)
#Stiamo scegliendo un DN40 con diam esterno = 48.3 mm e spessore 1.65mm
D_asp_G5=0.045                                                                  # m diametro della tubazione in aspirazione
u_asp_G5=((Fc1*(1+c1a+c1b)+Pc1)/(rho_R1*3600))/(np.pi*0.25*D_asp_G5**2)         # m/s

#ipotizziamo velocità
u_mand_G5=1.5                                                                   # m/s tentativo
D_mand_G5=((4/np.pi)*((Fc1*(1+c1a+c1b)+Pc1)/(rho_R1))/(u_mand_G5*3600))**(0.5)
# #Stiamo scegliendo un DN25 con diam esterno = 33.4 mm e spessore 1.65 mm
D_mand_G5=0.0301                                                                # m diametro della tubazione in mandata
u_mand_G5=((Fc1*(1+c1a+c1b)+Pc1)/(rho_R1*3600))/(np.pi*0.25*D_mand_G5**2)       # m/s

Q_G5=(Fc1*(1+c1a+c1b)+Pc1)/(rho_R1)                                             # portata volumetrica m3/h della pompa 
Re_asp_G5=rho_R1*u_asp_G5*D_asp_G5/(mu0)                                        # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
f_asp_G5=0.0791*Re_asp_G5**(-0.25)
deltaPdist_asp_G5=4*f_asp_G5*((z1_G5+l1_G5)/D_asp_G5)*0.5*rho_R1*u_asp_G5**2    # Pa Perdite di carico distribuite
deltaPconc_asp_G5=0.5*(K_90+K_in)*rho_R1*u_asp_G5**2                            # Pa perdite di carico concentrate
deltaPtotasp_G5=deltaPdist_asp_G5+deltaPconc_asp_G5                             # Pa perdite di carico totali nel tratto in aspirazione

Re_mand_G5=rho_R1*u_mand_G5*D_mand_G5/(mu0)                                     # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
f_mand_G5=0.0791*Re_mand_G5**(-0.25)
deltaPdist_mand_G5=4*f_mand_G5*((z2_G5+l2_G5)/D_mand_G5)*0.5*rho_R1*u_mand_G5**2# Pa perdite di carico distribuite
deltaPconc_mand_G5=0.5*(2*K_90+K_exit)*rho_R1*u_mand_G5**2                      # Pa perdite di carico concentrate
deltaPtotmand_G5=deltaPdist_mand_G5+deltaPconc_mand_G5                          # Pa perdite di carico totali mandata
deltaPtot_G5=deltaPtotasp_G5+deltaPtotmand_G5                                   # Pa perdite di carico totali aspirazione più mandata

Hp_G5=(PB_G5+rho_R1*9.81*z2_G5+deltaPtot_G5-PA_G5-0.5*rho_R1*u_asp_G5**2)/(rho_R1*9.81) #m prevalenza della pompa

# TCA 32-125 DA GRAFICO 1450 rpm 50Hz 4 poli phi=135 (teorico)
NPSHd_G5=1/(rho_R1*9.81)*(PA_G5-deltaPtotasp_G5+0.5*rho_R1*u_asp_G5**2-Psat_R1) # m
NPSHr_G5=0.55                                                                    # m

PotN_G5=0.23                                                                    # kW potenza nominale
c_G5=rho_R1/(Chemical('H2O',T=273.15+Tc1)).rhol                                 # fattore correttivo per la densità della soluzione
PotR_G5=c_G5*PotN_G5                                                            # kW Potenza con fattore correzione
eta_G5=0.28                                                                     # efficienza
PotRic_G5=PotR_G5/eta_G5 #kW                                                    # kW potenza richiesta dalla pompa
Pasp_G5=NPSHd_G5*rho_R1*9.81 - 0.5*rho_R1*u_asp_G5**2 + Psat_R1                 # Pa Pressione assoluta in aspirazione

############################ CRISTALLIZZATORE Na2SO4
# in C2 entra Fc1
tau_C2=1                                    # hr tempo di residenza della soluzione dentro il cristallizzatore
Vl_C2=(Fc1*(1+c1a+c1b)/rho_R1)*tau_C2       # m3 volume di liquido all'interno del cristallizzatore

#Facciamo l'ipotesi che sia pieno al 40%
Vc_C2=Vl_C2/0.4                             # m3 

#ipotizziamo un diametro Dc2
D_C2=1.8                                    # m diametro del cristallizzatore
h_C2=Vc_C2/(0.25*np.pi*D_C2**2)             # m altezza del cristallizzatore
hL_C2=0.4*h_C2                              # m altezza liquido nel cristallizzatore

############################ Pompa pappa salina G8 (Cristallizzatore 2)
#A è in fondo alla tubazione verticale che scende dal cristallizzatore
#B è all'ingresso dell'addensatore

z1_G8=2.5                                   # m
l1_G8=2.5                                   # m
z2_G8=4                                     # m
l2_G8=1.5                                   # m
PB_G8=101325                                # Pa
PA_G8=Psat_R2+rho_R2*9.81*z1_G8+rho_R2*9.81*hL_C2 #Pa

#perdite di carico

#ipotizziamo velocità
u_asp_G8=0.6                                                                    # m/s tentativo
D_asp_G8=((4/np.pi)*((Fc2*(1+c2a+c2b)+Pc2)/(rho_R2))/(u_asp_G8*3600))**(0.5)
#Stiamo scegliendo un DN32 con diam esterno = 42.2 mm e spessore 1.65mm
D_asp_G8=0.0389                                                                 # m
u_asp_G8=((Fc2*(1+c2a+c2b)+Pc2)/(rho_R2*3600))/(np.pi*0.25*D_asp_G8**2)         # m/s

#ipotizziamo velocità
u_mand_G8=1.5                                                                   #m/s tentativo
D_mand_G8=((4/np.pi)*((Fc2*(1+c2a+c2b)+Pc2)/(rho_R2))/(u_mand_G8*3600))**(0.5)
#Stiamo scegliendo un DN20 con diam esterno = 26.7 mm e spessore 1.65 mm
D_mand_G8=0.0234                                                                 #m
u_mand_G8=((Fc2*(1+c2a+c2b)+Pc2)/(rho_R2*3600))/(np.pi*0.25*D_mand_G8**2)        #m/s

Q_G8=(Fc2*(1+c2a+c2b)+Pc2)/(rho_R2)                                             #m3/h portata volumetrica pappa salina G8
Re_asp_G8=rho_R2*u_asp_G8*D_asp_G8/(mu0)                                        
f_asp_G8=0.0791*Re_asp_G8**(-0.25)                                              # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
deltaPdist_asp_G8=4*f_asp_G8*((z1_G8+l1_G8)/D_asp_G8)*0.5*rho_R2*u_asp_G8**2    # Pa perdite di carico distribuite
deltaPconc_asp_G8=0.5*(K_90+K_in)*rho_R2*u_asp_G8**2                            # Pa perdite di carico concentrate
deltaPtotasp_G8=deltaPdist_asp_G8+deltaPconc_asp_G8                             # Pa perdite di carico totali nel tratto di aspirazione

Re_mand_G8=rho_R2*u_mand_G8*D_mand_G8/(mu0)
f_mand_G8=0.0791*Re_mand_G8**(-0.25)                                            # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
deltaPdist_mand_G8=4*f_mand_G8*((z2_G8+l2_G8)/D_mand_G8)*0.5*rho_R2*u_mand_G8**2# Pa perdite di carico distribuite
deltaPconc_mand_G8=0.5*(2*K_90+K_exit)*rho_R2*u_mand_G8**2                      # Pa perdite di carico concentrate
deltaPtotmand_G8=deltaPdist_mand_G8+deltaPconc_mand_G8                          # Pa perdite di carico totali nel tratto di mandata
deltaPtot_G8=deltaPtotasp_G8+deltaPtotmand_G8                                   # Pa perdite di carico totali aspirazione + mandata

Hp_G8=(PB_G8+rho_R2*9.81*z2_G8+deltaPtot_G8-PA_G8-0.5*rho_R2*u_asp_G8**2)/(rho_R2*9.81) #m prevalenza della pompa

# TCA 32-125 DA GRAFICO 1450 rpm 50Hz 4 poli phi (122) teorico
NPSHd_G8=1/(rho_R2*9.81)*(PA_G8-deltaPtotasp_G8+0.5*rho_R2*u_asp_G8**2-Psat_R2)
NPSHr_G8=0.5

PotN_G8=0.18                                                                    # kW potenza nominale
c_G8=rho_R2/(Chemical('H2O',T=273.15+Tc2)).rhol                                 # fattore correttivo per la densità della soluzione
PotR_G8=c_G8*PotN_G8                                                            # kW Potenza con efficienza
eta_G8=0.19                                                                     # efficienza
PotRic_G8=PotR_G8/eta_G8                                                        # kW potenza richiesta dalla pompa
Pasp_G8=NPSHd_G8*rho_R2*9.81 - 0.5*rho_R2*u_asp_G8**2 + Psat_R2                 # Pa Pressione assoluta in aspirazione

########## ESSICCATORE PNEUMATICO B1 ##########################
# Dati Input Pc1, densità del solido, Uin, Uout, dp, tetaGa
# Tempi si aggirano tra i 30 minuti e 1 hr

teta_bu=35                                          # °C T bulbo umido
teta_gi=89                                          # °C T gas ingresso
NTU=1.5 #NTU 1.5 e 3                                # unità teoriche
eff=1-np.exp(-NTU)                                  # efficienza
teta_go=teta_gi-eff*(teta_gi-teta_bu)               # °C T gas uscita
Hgi=0.015                                           # kg H2O/kg aria umidità assoluta ingresso 
Hgo=0.032                                           # kg H2O/kg aria umidità assoluta uscita 

Cp_a=Mixture(['air'],T=273.15+teta_gi).Cpg/1000     #kj/kg K
Cp_w=Chemical('water',T=273.15+Tc1).Cpl/1000        #kj/kg K
Hvap_bu=Chemical('water',T=273.15+teta_bu).Hvap/1000#kj/kg
Cp_v=Chemical('water',T=273.15+Tc1).Cpg/1000        #kj/kg K

G=(W1*(Cp_w*(teta_bu-Tc1)+Hvap_bu+Cp_v*(teta_go-teta_bu))+Pc1*Cp_w*Uo1*(teta_bu-Tc1))/((Cp_a+Cp_v*Hgi)*(teta_gi-teta_go)) #kg/h di aria entrante nell'essiccatore
# controllare che 1<G/Pc1<20

# #SEZIONE ORIZZONTALE ESSICCATORE B1
dp_medio=0.00025                                                                    #m diametro cristalli medio
vG_min=132.4*(Chemical('KClO3').rhos/(Chemical('KClO3').rhos+998))*dp_medio**(0.4)  # m/s tratto orizzontale (velocità di saltation)
teta_media=(teta_gi+teta_go)/2                                                      # °C T media gas
air=Mixture(['air'],T=273.15+teta_media)
vinf=1.5                                                                            # m/s velocità caduta libera di tentativo --> calcoliamo Re --> andiamo su grafico Perry 9 6-52 --> troviamo CD
Re_part=air.rhog*vinf*dp_medio/(air.mug)                                            
CD=4.6                                                                              # fattore di attrito Drag coefficent --> ricalcoliamo vinf--> confrontiamo con vinf iniziale e reiteriamo fino a convergenza
vinf=(4*9.81*dp_medio*(Chemical('KClO3').rhos-air.rhog))/(3*air.rhog*CD)            # m/s
vslip=vinf/(1.535+(vinf/(4.882+Chemical('KClO3').rhos))**(0.5))                     # m/s velocità di trascinamento
vG=10*vinf                                                                          # m/s velocità del gas
vS=vG-vslip                                                                         # m/s velocità del solido

A_B1=G/(air.rhog*3600*vG)                                                           # m2 sezione trasversale essiccatore
d_B1=(4*A_B1/(np.pi))**(0.5)                                                        # m diametro sezione trasversale essiccatore--> andiamo su catalogo

#CALCOLIAMO IL TEMPO
deltaT_ML_B1=(teta_gi-teta_go)/(np.log((teta_gi-teta_bu)/(teta_go-teta_bu)))                # °C
tau_B1=(1000*Hvap_bu*Chemical('KClO3').rhos*(W1)*dp_medio**2)/(12*air.k*deltaT_ML_B1*Pc1)   # s
L_B1=tau_B1*vS                                                                              # m lunghezza dell'essiccatore
Vol_B1=L_B1*A_B1                                                                            # m3 volume essiccatore
Hs_B1=(Pc1/3600*tau_B1)/(Chemical('KClO3').rhos*Vol_B1)*100                                 # % di Hold-up

#PERDITE DI CARICO
deltaPacc_B1=((Pc1/3600)/A_B1)*vS/(9.81)                                                    #Pa perdite di carico accelerazione del solido
Re_B1=d_B1*air.rhog*vG/air.mug
f_B1=0.331/(np.log(epsi/(3.7*d_B1)+7/Re_B1))**2                                             # coefficiente di attrito
l_bends_B1=40*d_B1                                                                          # fattore legata alle perdite di carico concentrate
l_tot_B1=l_bends_B1+L_B1                                                                    # poi bisognerà stimare orizzontale e verticale
deltaPg_B1=(4*f_B1*l_tot_B1*air.rhog*vG**2)/(2*9.81*d_B1)                                   # Pa perdite di carico dovute al gas 

#Facciamo 27m in orizzontale e circa 3m in verticale
lo_B1=27                                                                                    # m lunghezza orizzontale dell'essiccatore
lv_B1=L_B1-lo_B1                                                                            # m lunghezza verticale dell'essiccatore
Gc_B1=9.807
deltaHg_B1=(lv_B1*air.rhog*9.81)/(Gc_B1)                                                    # Pa perdite di carico dovute all'elevazione
deltaHs_B1=(lv_B1*((Pc1/3600)/A_B1)*9.81)/(vS*Gc_B1)                                        # Pa perdite di carico dovute all'elevazione

# # Si valutano le deltaP misc vedi Theory and Design of Dilute Phase Pneumatic Conveying Systems Agarwal
deltaP_HE7=15000                                                                            # Pa sono le perdite di carico massime lato mantello scambiatore alettato 
deltaP_cyc=1000                                                                             # Pa 100mm H2O
deltaP_filtro=2000                                                                          # Pa 200 mm H2O
l_onlygas=15                                                                                # m decisa arbitrariamente, lunghezza del tratto che l'aria attraversa prima di incontrare il solido
deltaP_distG=4*f_B1*(l_onlygas/d_B1)*0.5*air.rho*vG**2                                      # Pa perdite di carico distribuite
deltaP_misc=deltaP_HE7+deltaP_cyc+deltaP_filtro+deltaP_distG                                # Pa perdite di carico di ogni apparecchiatura lungo la linea 
deltaP_tot_B1=deltaP_misc+deltaHg_B1+deltaHs_B1+deltaPacc_B1+deltaPg_B1                     # Pa perdite di carico totali

# VALUTAZIONE VENTILATORE
H_vent=(deltaP_tot_B1+0.5*vG**2*(Mixture(['air'],T=273.15+teta_go).rhog-Mixture(['air'],T=298.15).rhog))/(air.rhog*9.81) # m prevalenza del ventilatore
G_vent=G/(air.rhog)                                                                                                      # m3/h portata del ventilatore 

#MZ-aspiratori VC 630-N 4550rpm
PotN_V=11                   # kW potenza nominale
eta_V=0.55                  # efficienza
PotR_V=PotN_V/eta_V         # kW potenza reale del ventilatore

#VALUTAZIONE ALTEZZA 1 (RELATIVA TRA CONCENTRATORE E C1)
vb_H1=u_0C                          # m/s
Re_H1=(rho_R0*vb_H1*d_0C)/(mu0)
f_H1=0.0791*Re_H1**(-0.25)          # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
ha_H1=4                             #m
z1_H1=1                             #m
l1_H1=4                             #m
deltaPd_H1=4*f_H1*((ha_H1-l_F0+l1_H1+2*z1_H1)/(d_0C))*0.5*rho_R0*vb_H1**2   # Pa perdite di carico distribuite
deltaPc_H1=0.5*(3*K_90+K_in+K_exit)*rho_R0*vb_H1**2                         # Pa perdite di carico concentrate
deltaPtot_H1=deltaPd_H1+deltaPc_H1                                          # Pa perdite di carico totali
zbatt_H1=1                                                                  # m altezza battente di liquido data dalla tubazione
batt_H1=hL_C1+zbatt_H1                                                      # m (distanza tra fondo C1 e B)
Pbmin_H1=Psat_R1+rho_R1*9.81*batt_H1                                        # Pa Pressione minima affinchè possa entrare nella sezione di cristallizzazione 1
ha_H1=(Pbmin_H1+0.5*rho_R0*vb_H1**2+deltaPtot_H1-Psat_R0)/(rho_R0*9.81)     # m altezza del punto a della sezione H1 (alimentazione Cristallizzatore 1)

# #VALUTAZIONE ALTEZZA 2 (RELATIVA TRA SERBATOIO CONDENSE E C2)

Pc_H2=101325                                                        # Pa
vd_H2=1                                                             # m/s ipotesi
d_H2=((4/np.pi)*((Fs1*(1+c1a+c1b))/(rho_R1*3600))/(vd_H2))**0.5     # m
# Stiamo scegliendo un DN32 con diam esterno = 42.2mm e spessore 1.65mm
d_H2=0.0389                                                          # m
vd_H2=(Fs1*(1+c1a+c1b)/(rho_R1*3600))/(np.pi*0.25*d_H2**2)          # m/s velocità effettiva 

Re_H2=rho_R1*vd_H2*d_H2/(mu0)
f_H2=0.0791*Re_H2**(-0.25)                                          # controlla: se Re>10^5 --> tabella, altrimenti f=0.0791*Re^(-0.25)
z1_H2=1                                                             # m
l1_H2=3                                                             # m
deltaPd_H2=4*f_H2*((2*z1_H2+l1_H2)/(d_H2))*0.5*rho_R1*vd_H2**2      # Pa perdite di carico distribuite
deltaPc_H2=0.5*(3*K_90+K_in+K_exit)*rho_R1*vd_H2**2                 # Pa perdite di carico concentrate
deltaPtot_H2=deltaPd_H2+deltaPc_H2                                  # Pa perdite di carico totali
zbatt_H2=1                                                          # m altezza battente di liquido data dalla tubazione
batt_H2=hL_C2+zbatt_H2                                              # m (distanza tra fondo C1 e B)
Pdmin_H2=Psat_R2+rho_R2*9.81*batt_H2                                # Pa Pressione minima affinchè possa entrare nella sezione di cristallizzazione 2
hc_H2=(Pdmin_H2+0.5*rho_R1*vd_H2**2+deltaPtot_H2-Pc_H2)/(rho_R1*9.81)#  m altezza del punto c della sezione H1 (alimentazione Cristallizzatore 1)

########### ESSICCATORE ROTATIVO B2 ##########################
#Dati Input Pc1, densità del solido, Uin, Uout, dp, tetaGa
#Tempi si aggirano tra i 30 minuti e 1 hr

teta_bu_B2=35                                                       # °C T bulbo umido
teta_gi_B2=89                                                       # °C T gas ingresso
NTU_B2=1.5                                                          # Unità teoriche NTU 1.5 e 3
eff_B2=1-np.exp(-NTU_B2)                                            # efficienza
teta_go_B2=teta_gi_B2-eff_B2*(teta_gi_B2-teta_bu_B2)                # °C T gas uscita
Hgi_B2=0.015                                                        # kg H2O/ kg aria secca umidità assolute
Hgo_B2=0.032                                                        # kg H2O/ kg aria secca umidità assolute

Cp_a_B2=Mixture(['air'],T=273.15+teta_gi_B2).Cpg/1000               #kj/kg K
Cp_w_B2=Chemical('water',T=273.15+Tc2).Cpl/1000                     #kj/kg K
Hvap_bu_B2=Chemical('water',T=273.15+Tc2).Hvap/1000                 #kj/kg
Cp_v_B2=Chemical('water',T=273.15+Tc2).Cpg/1000                     #kj/kg K
G_B2=(W2*(Cp_w_B2*(teta_bu_B2-Tc2)+Hvap_bu_B2+Cp_v_B2*(teta_go_B2-teta_bu_B2))+Pc2*Cp_w_B2*Uo2*(teta_bu_B2-Tc2))/((Cp_a_B2+Cp_v_B2*Hgi_B2)*(teta_gi_B2-teta_go_B2)) #kg/h gas nell'essiccatore

Gs_B2=1.5                       # kg/s m2 portata di gas massica per unità di area
A_B2=(G_B2/3600)/(Gs_B2)        # m2 area di passaggio del gas
d_B2=(4*A_B2/np.pi)**(0.5)      # m tra 0.3 e 3

#Correlazioni valide solo con un certo range di hsa (hold up di solido)
Ua_B2=(237*Gs_B2**0.67)/(d_B2)                                                                  # W/m3 K coefficiente di scambio 
Cs_B2=Cp_a_B2+Cp_v_B2*Hgi_B2                                                                    # kj/kg C calore specifico aria 
deltaTML_B2=(teta_gi_B2-teta_go_B2)/np.log((teta_gi_B2-teta_bu_B2)/(teta_go_B2-teta_bu_B2))     # °C
V_B2=(G_B2/3.600)*Cs_B2*(teta_gi_B2-teta_go_B2)/(Ua_B2*deltaTML_B2)                             # m3 volume di scambio
l_B2=V_B2/A_B2                                                                                  # m lunghezza dell'essiccatore
forma_B2=l_B2/d_B2

if forma_B2>10:
    raise ValueError('Il fattore di forma dell essiccatore è maggiore di 10!!')
elif forma_B2<4:
    raise ValueError('Il fattore di forma dell essiccatore è minore di 4!!') 

HTU_B2=1000*Gs_B2*Cs_B2/(Ua_B2)                                                                 # m altezza unità teoriche
NTU_B2=l_B2/HTU_B2                                                                              # unità teoriche
vp_B2=35                                                                                        # m/min velocità periferica del cilindro
N_B2=vp_B2/(60*np.pi*d_B2)                                                                      # rps del cilindro
s_B2=0.04                                                                                       # inclinazione del tamburo tan(alpha)
p_B2=d_B2*s_B2                                                                                  # m passo (tratto di avanzamento delle particelle)
beta_B2=0.35                                                                                    # coefficiente
peff_B2=p_B2*beta_B2                                                                            # m passo effettivo (tratto avanzamento delle particelle)

va_B2=N_B2*peff_B2                                                                              # m/s velocità di avanzamento 
tau_B2=(l_B2/va_B2)/60                                                                          # min tempo di essiccamento

Hs_B2=Pc2/(3600*Chemical('Na2SO4').rhos)*(tau_B2*60/V_B2)                                       # hold up di solido
epsilon_B2=0.35                                                                                 # grado di vuoto
rhosa_B2=Chemical('Na2SO4').rhos*(1-epsilon_B2)                                                 # kg/m3 densità di solido apparente
Hsa=Hs_B2*Chemical('Na2SO4').rhos/rhosa_B2                                                      # hold up apparente

print('F0 =',F0)
print('Fc1 =',Fc1)
print('Ec1 =',E1)
print('Fs1 =',Fs1)
print('Fc2 =',Fc2)
print('Ec2 =',E2)
print('Fs2 =',Fs2)
print('Pc2 =',Pc2)
print('F =',F)
print('Fc =',Fc)
print('E0 =',E0)
print('ResaC1 =', ResaC1)
print('ResaC2 =', ResaC2)
print('Resa G =', ResaG)
print('Qc =',Q_c)
print('Qc1 =',Q_c1)
print('Qc2 =',Q_c2)
print('deltaTeb0 =', deltaT_eb0 )
print('deltaTeb1 =', deltaT_eb1 )
print('deltaTeb2 =', deltaT_eb2 )
print('Fa1 =', Fa1)
print('Fad1 =', Fad1)
print('Fsep1 =', Fsep1)
print('Fa2 =', Fa2)
print('Fad2 =', Fad2)
print('Fsep2 =', Fsep2)
print('R0 =',R0)
print('Prevalenza G1 =',Hp_G1)
print('NPSH richiesto G1 =',NPSHr_G1)
print('NPSH disponibile G1 =',NPSHd_G1)
print('Potenza Richiesta G1 =', PotRic_G1)
print('Q aspirazione nominale G2 =',Qasp_n_G2)
print('Q aspirazione effettiva G2 =',Qeff_G2)

# print('G/Pc1 =',G/Pc1)
# print('vG =',vG)
# print('tau_B1 =', tau_B1)
# print('Hs_B1 =',Hs_B1)
# print('Hvent = ',H_vent)
# print('G_vent =', G_vent)
# print('d_B2 =', d_B2)
# print('rapporto forma essiccatore B2 =',l_B2/d_B2)
# print('tau B2 =', tau_B2)
# print('Hsa =', Hsa)

