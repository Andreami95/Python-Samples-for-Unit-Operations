import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from thermo.chemical import Chemical, Mixture, VolumeLiquidMixture
from scipy import integrate
from Pipes_Steel import *
from Pipes_Inox import * 

plt.close('all')

#SHELL SIDE
Ts_in = 113 + 273.15
Ts_out = 113 + 273.15
P_Shell = 1.6e5
MShellFluid = 27.6/3 #kg/s
ShellFluid = Mixture(['CH3OH','H2O'],zs=[0.01,0.99],T=Ts_in,P=P_Shell)

#Tube SIDE
Tt = Chemical('H2O').Tsat(101325*3)
Tt = Chemical('H2O').Tsat(101325*3)
P_Tube = 101325*3
TubeFluid = Chemical('H2O',T =Tt,P=P_Tube)
# deltat = Tt_out - Tt_in

def deltaT(Ts_in,Ts_out,Tt_in,Tt_out,flow): #flow=co-->co-current flow=counter-->counter-current
    
    if flow=='co':
        #co-current
        deltaTML    = ((Ts_in-Tt)-(Ts_out-Tt))/np.log((Ts_in-Tt)/(Ts_out-Tt))
        r           = (Ts_out-Tt)/(Ts_in-Tt) #parameter for Tcal
    
    elif flow=='counter': 
        #counter-current
        deltaTML    = ((Ts_in-Tt)-(Ts_out-Tt))/(np.log((Ts_in-Tt)/(Ts_out-Tt)))
        r           = (Ts_out-Tt)/(Ts_in-Tt) #parameter for Tcal
        
    # R = (Ts_in-Ts_out)/(Tt_out-Tt_in) #parameters for Correction Factor
    # S = (Tt_out-Tt_in)/(Ts_in-Tt_in)  #parameters for Correction Factor
    # print('R =',R)
    # print('S =',S)
    # print('Insert Correction Factor:')
    # link to Perry Correction Factors
    # decomment line to open browser
    # webbrowser.open('https://www.slideshare.net/claumontero93/perrys-chemical-engineers-handbook-7ma-ed-chap-11')
    # CorrectionFactor=float(input())
    CorrectionFactor = 1 #insert num grafico C 3 shell passes 6 or more tube passes
    deltaTeff = deltaTML*CorrectionFactor
    return deltaTML,deltaTeff,r


deltaT      = deltaT(Ts_in,Ts_out,Tt,Tt,'counter')
deltaT_ML   = deltaT[0]
print('\u0394T ML =','%.2f' %deltaT_ML)

deltaT_eff  = 128-113
print('\u0394T eff =','%.2f' %deltaT_eff)

plt.figure(1)
plt.title('Temperature Profile')
plt.grid('on')
plt.xlim(0,1)
plt.plot([0,1],[Ts_in-273.15,Ts_out-273.15],color='r',label='Shell Fluid') #Process Fluid
plt.plot([0,1],[Tt-273.15,Tt-273.15],color='b',label='Tube Fluid') #Service Fluid
plt.ylabel('°C')
plt.xlabel('Normalized Tube Length')
plt.legend()

# # Q (HEAT EXCHANGED)
MEvap       = 6.38 #kg/s
Qnom        = MEvap*(ShellFluid.Cpl*(Ts_in-Ts_out)+ShellFluid.Hvaps[1])/1000 #kW
MTubeFluid  = Qnom/(TubeFluid.Hvap)*1000 #kg/s

U       = np.zeros(1000)
h_i     = np.zeros(1000)
h_e     = np.zeros(1000)
h_ie    = np.zeros(1000)

# U[0]=input('U first try = ') #W/m2*K (iteration start)
U[0]=600
h_ie[0]=100         # first try
h_e[0]=1100         # first try
h_i[0]=10000        # first try

#Caloric Temperature
for i in range(1000):
    # print(i)
    Kc  = (h_e[i]-h_ie[i])/(h_ie[i])               # ci sarà un processo iterativo Uc=coeffieciente scambio lato freddo Uh coefficiente di scambio lato caldo
    r   = deltaT[2]                                 # deltaT lato freddo/ deltaT lato caldo deve essere compreso tra 0.01 e 10 (serve per calcolare T calorica)
    # Fc  = 1#(1/Kc+(r/(r-1)))/(1+(np.log(Kc+1))/np.log(r))-1/Kc
    
    #T caloriche
    Ts_cal = Ts_out+0.5*(Ts_in-Ts_out)
    Tt_cal = Tt
    
    # # CALCOLO FLUSSO DI MASSA E NUMERO DI REYNOLDS LATO TUBI E MANTELLO
    SafetyFactor    = 1.07
    Qscambio        = Qnom*SafetyFactor #kW
    Ascambio        = 1000*Qscambio/(U[i]*deltaT_eff) #m2
        
    #si sceglie la geometria dei tubi ed il materiale
    De          = DN20i.d*1e-3 #m diametro esterno
    spessore    = DN20i.t[0]*1e-3 
    Di          = De-2*spessore #m diametro interno 
    
    #vediamo a che velocità passa il fluido all'interno dei tubi
    Tube_Length = 6 #m
    Ntubi       = np.round(0.9*Ascambio/(np.pi*De*Tube_Length),0)
    Ascambio    = 1.2*Ntubi*np.pi*Tube_Length*De
    Tube_passes = 2
    Apassaggiotubi = (np.pi*Di**2)/4*(Ntubi/Tube_passes) #m2 sezione di passaggio fluido lato tubi
    Shell_passes = 1

    
    TubeFluid       = Chemical('H2O',T=Tt_cal+0.01,P=P_Tube) # kg/m3 per Serv Fluid
    Vol_TubeFluid   = MTubeFluid/TubeFluid.rhog #m3/s portata volumica fluido lato tubi
    v_TubeFluid     = Vol_TubeFluid/Apassaggiotubi #m/s
    Re_tube         = TubeFluid.rhog*v_TubeFluid*Di/TubeFluid.mug
    K_Pitch         = 1.25
    Pitch           = K_Pitch*De
    
    # #Calcoliamo l'area laterale di passaggio lato mantello
    TubeConfig='tri'
    def D_Bundle(Npasses,Ntubi,config):
        if config=='tri':
            if Npasses==1:
                K_1=0.319
                n_1=2.142
            elif Npasses==2:
                K_1=0.249
                n_1=2.207
            elif Npasses==4:
                K_1=0.175
                n_1=2.285
            elif Npasses==6:
                K_1=0.0743
                n_1=2.499
            elif Npasses>=8:
                K_1=0.0365
                n_1=2.675
            D_b=De*(Ntubi/K_1)**(1/n_1)
        
        if config=='quad':
            if Npasses==1:
                K_1=0.215
                n_1=2.207
            elif Npasses==2:
                K_1=0.156
                n_1=2.291
            elif Npasses==4:
                K_1=0.158
                n_1=2.263
            elif Npasses==6:
                K_1=0.0402
                n_1=2.617
            elif Npasses>=8:
                K_1=0.0331
                n_1=2.643
            D_b=De*(Ntubi/K_1)**(1/n_1)
        return D_b
    
    Bundle_diameter=np.around(D_Bundle(Tube_passes,Ntubi,TubeConfig),4) #m 
         
    Shell_diameter  = Bundle_diameter*1.8
    FormFactor      = Tube_Length/Bundle_diameter
    
    Ntubes_centralrow   = Bundle_diameter/Pitch
    
    def V_cell(config):
        if config=='tri':
            Vcell       = 1*Pitch*Pitch*np.sin(np.pi/3)/2               #m3 volume cella (per metro di tubo)
            Vcell_free  = Vcell-(0.5*1*np.pi*De**2)/4
        elif config=='quad':
            Vcell       = Pitch**2
            Vcell_free  = Vcell-(1*np.pi*De**2)/4
        return Vcell, Vcell_free
        
    #CALCOLO COEFFICIENTI DI SCAMBIO
    #LATO TUBI
    # bisogna considerare se passa acqua all'interno dei tubi (abbiamo una correlazione apposita)
    # bisogna considerare anche il Re del fluido se siamo laminari o turbolenti
    
    #check per vedere se è acqua
    h_i[i]  = 10000 # W/m2*K (diametro tubo interno in mm e temperatura in C)
    Nu      = h_i[i]*Di/TubeFluid.k
    h_ie[i] = h_i[i]*Di/De    
    
    #FOULING FACTOR
    # Internal Resistance dirt (m^2°C/W)
    R_i = 3e-04 #River water, Cooling water (towers),Towns water (soft),Refrigerated brine
    # R_i = 1e-03 #Sea water
    # R_i = 2e-04 #Air and industrial gases
    # R_i = 2e-04 #Organic liquids and organic vapours
        
    #LATO MANTELLO
    
    # The baffle cut is the height of the segment removed to form the baffle, expressed as a percentage of
    # the baffle disc diameter.
    
    # External Resistance dirt (m^2°C/W)
    # R_e = 3e-04 #River water, Cooling water (towers),Towns water (soft),Refrigerated brine
    # R_e = 1e-03 #Sea water
    R_e = 2e-04 #Air and industrial gases
    # R_e = 2e-04 #Organic liquids and organic vapours
    
    # j_hmantello15 = lambda x: 0.8584*x**(-0.5942)+0.00279 #15% baffle cut  # COULSON PAG.672 (verificato) NON PRECISISSIMO SOPRA Re = 1e5
    # j_hmantello25 = lambda x: 0.7028*x**(-0.5814)+0.002585 #25% baffle cut
    # j_hmantello35 = lambda x: 0.6139*x**(-0.5692)+0.001789 #35% baffle cut
    # j_hmantello45 = lambda x: 0.5885*x**(-0.5733)+0.001621 #45% baffle cut
    
    rhog            = Mixture(['CH3OH','H2O'],zs=[0.01,0.99],T=115+273.15,P=1.6e5).rhog
    Tsat            = Mixture(['CH3OH','H2O'],zs=[0.01,0.99],T=115+273.15,P=1.6e5).Tbubble
    T_wall          = (h_e[i]*Tsat+h_ie[i]*Tt_cal)/(h_e[i]+h_ie[i])
    num = (ShellFluid.kl**0.79)*(ShellFluid.Cpl**0.45)*(ShellFluid.rhol**0.49)
    den = (ShellFluid.sigmas[1]**0.5)*(ShellFluid.mul**0.29)*(ShellFluid.Hvaps[1]**0.24)*(rhog**0.24)
    P_wall = Chemical('H2O',T=115+273.15,P=1.6e5).Psat
    Psat = 1.6e5
    h_e[i] = 0.00122*(num/den)*(T_wall-Tsat)**0.24*(P_wall-Psat)**0.75
    
    # print(h_e[i])
#     #CONDUCIBILITA' TERMICA METALLI k (W/m°C)
    k = 22  #Acciaio inox (20% Cr)
#     # k = 26  #Bronzo (75% Cu, 25% Sn)
#     # k = 36  #Acciaio al carbonio (1,5% C)
#     # k = 54  #Acciaio al carbonio (0,5% C)
#     # k = 73  #Ferro
#     # k = 111  #Ottone (70% Cu, 30% Zn)
#     # k = 204 #Alluminio
#     # k = 386 #Rame
    
    U[i] = (1/h_e[i]+1/h_ie[i]+spessore/k+R_i + R_e)**-1
    U[i] = '%.1f' %U[i]
    Q_eff=U[i]*Ascambio*deltaT_eff/1000 #kW
        
    tol=0.1
    
    if np.abs(U[i]-U[i-1])<tol:
        
        Uo      = U[i]          # Overall heat transfer coefficient 
        H_i     = h_i[i]        # Internal heat transfer coefficient
        H_ie    = h_ie[i]       # int --> ext heat transfer coefficient
        H_e     = h_e[i]        # External heat transfer coefficient
        break
    else:
        U[i+1]      = U[i]
        h_i[i+1]    = h_i[i]
        h_e[i+1]    = h_e[i]
        h_ie[i+1]   = h_ie[i]


if Q_eff>1.2*Qnom:
    print('Q_eff =', Q_eff)
    raise ValueError('Q eff > 120%!')
elif Q_eff>Qnom:
    print('---------------------------------------------------------------------------')
    print('Thermodynamic Design is verified')
    print('Q =','%.2f' %Q_eff, 'kW','Q =','%.2f' %(100*Q_eff/Qnom),'% of Q nom')
    print('U =',Uo,'W/(m^2 K)')
    print('Total Area =','%.2f'%(Ascambio),'m2')
    print('Hi =', '%.2f' %H_i, 'W/(m^2 K)')
    print('Hie =', '%.2f' %H_ie, 'W/(m^2 K)')
    print('He =', '%.2f' %H_e, 'W/(m^2 K)')
    print('Total tubes =','%.0f' %Ntubi)
    print('V fluid tube side =','%.2f' %v_TubeFluid,'m/s')
    # print('V fluid shell side =','%.2f' %v_ShellFluid,'m/s')
    print('Form factor =','%.1f' %FormFactor)
    print('Tube Length =','%.1f' %Tube_Length, 'm')
    print('Shell Diameter =', '%.3f' %Shell_diameter, 'm')
else:
    raise ValueError('Thermodynamic Design is NOT verified')

Heat_Flux = Qnom/Ascambio
Heat_Flux_Critical = (0.131*ShellFluid.Hvaps[1]*((ShellFluid.sigmas[1]*9.81*(ShellFluid.rhol-rhog)*rhog**2))**0.25)/1000
# Heat_Flux_Critical = (3.67e4*(ShellFluid.Pc/101325)*(P_Shell/ShellFluid.Pc)**0.35*(1-(P_Shell/ShellFluid.Pc))**0.9)/1000
# Heat_Flux_Critical = 38 #kW da letteratura
print('Heat Flux = ','%.1f' %Heat_Flux,'kW/m2')
print('Critical Heat Flux = ','%.1f' %Heat_Flux_Critical,'kW/m2')

alpha = 0.2
h_liq = Bundle_diameter+alpha
teta=np.arccos(2*alpha/Shell_diameter)              #rad

corda=Shell_diameter*np.sin(teta)                   #m corda sul mantello al livello di liquido
A_liqs=corda*Tube_Length                            #m2 sezione di liquido
Vol_ShellFluid_vap = MEvap/rhog                     #m3/s
v_vaps=Vol_ShellFluid_vap/(A_liqs)                  #m/s velocità superficiale del vapore
v_vapmax=0.2*((ShellFluid.rhol-rhog)/rhog)**0.5     #m/s velocità massima del vapore
if v_vaps>v_vapmax:
    raise ValueError('Velocità del vapore elevata')

MLiq = 8.4/3
Vol_ShellFluid_liq = MLiq/ShellFluid.rhol #m3/s
beta=0.1                            #m in aggiunta all'altezza del weir
h_weir=h_liq                        #m altezza del weir
# h_sc=0.75*h_weir                 #m livello liquido di scarico
A_Shell = 0.25*np.pi*Shell_diameter**2

alpha2 = 0.5
h_sc = Shell_diameter/2 - alpha2
teta2=np.arccos(2*alpha2/Shell_diameter)              #rad

A_scarico   = (((Shell_diameter/2)**2)*(np.pi-teta2)+(Shell_diameter/2-(Shell_diameter-h_sc))*Shell_diameter/2)
A_tubi      = (((Shell_diameter/2)**2)*(np.pi-teta)+(Shell_diameter/2-(Shell_diameter-h_weir))*Shell_diameter/2)-(0.25*Ntubi*np.pi*De**2)*1.2
Vol_liq_tubi= A_tubi*Tube_Length*ShellFluid.rhol
L_sc=1                           #m lunghezza dello scarico
Vol_liq_scarico = A_scarico*0.2*ShellFluid.rhol
Vol_tot         = Vol_liq_tubi+ Vol_liq_scarico
Residence_Time  = (Vol_tot/Vol_ShellFluid_liq)/3600
Residence_Time2 = (Vol_liq_scarico/Vol_ShellFluid_liq)/3600

L_tot=Tube_Length*1.2+L_sc           #m lunghezza totale del reboiler
Am = (np.pi*Shell_diameter**2)/4             #m2 area del mantello

# Perdite di carico
# Lato tubi
f_dist          = 0.0791*Re_tube**-0.25
deltaPdist_tubi = 4*f_dist*(Tube_Length*Tube_passes/Di)*(0.5*TubeFluid.rhog*v_TubeFluid**2) #Pa

D_bocchellotubi = (DN750i.d-DN750.t[2]*2)/1000 #m
K_curva         = 1.5
v_TubeFluid_nozzle    = Vol_TubeFluid/(0.25*np.pi*D_bocchellotubi**2) #m/s
K_in            = 0.5*(1-(0.25*np.pi*Di**2)/(0.25*np.pi*D_bocchellotubi**2))
K_out           = 0.5*v_TubeFluid**2*(1-(0.25*np.pi*Di**2)/(0.25*np.pi*D_bocchellotubi**2))**2

deltaPconc_tubi = 0.5*K_curva*Tube_passes*(TubeFluid.rhog*v_TubeFluid**2)+0.5*(K_in+K_out)*(TubeFluid.rhog*v_TubeFluid_nozzle**2) #Pa
deltaPtot_tubi  = deltaPdist_tubi+deltaPconc_tubi

Vol_ShellFluid_vap = MEvap/rhog #m3/s
D_bocch_vap = (DN850.d-DN850.t[2]*2)/1000 #m
v_ShellFluid_vap_nozzle = Vol_ShellFluid_vap/(0.25*np.pi*D_bocch_vap**2)


MLiq_in = 8.4/3 + MEvap
Vol_ShellFluid_liq_in = MLiq_in/ShellFluid.rhol #m3/s
D_bocch_liq_in = (DN100i.d-DN100i.t[2]*2)/1000 #m
v_ShellFluid_liq_nozzle_in = Vol_ShellFluid_liq_in/(0.25*np.pi*D_bocch_liq_in**2)



MLiq_out = 8.4/3
Vol_ShellFluid_liq_out = MLiq_out/ShellFluid.rhol #m3/s
D_bocch_liq_out = (DN65i.d-DN65i.t[2]*2)/1000 #m
v_ShellFluid_liq_nozzle_out = Vol_ShellFluid_liq_out/(0.25*np.pi*D_bocch_liq_out**2)

print('---------------------------------------------------------------------------')
print('Pressure Drop')
print('Tube side Pressure Drop =', '%.2f' %(deltaPtot_tubi/1000),'kPa')
# print('Shell side Pressure Drop =', '%.2f' %(deltaPtot_Shell/1000),'kPa')