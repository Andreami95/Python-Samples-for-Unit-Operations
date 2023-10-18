import numpy as np
from thermo import Chemical
from thermo import Mixture
import matplotlib.pyplot as plt
# from pipes import Steel
# from Pipes_Steel import *
from Pipes_Inox import * 
import webbrowser

plt.close('all')

#HEAT EXCHANGER
# SHELL FLUID
Ts_in           = 130.3+273.15          #SHELL INLET
Ts_out          = 52+273.15             #SHELL OUTLET
P_Shell         = 4e6 
ShellFluid      = Mixture(['H2','CO2'],Vfgs=[0.75, 0.25],T=Ts_in,P=P_Shell)
MShellFluid     = 90532/3600 #kg/s

# TUBE FLUID
Tt_in           = 20+273.15             #TUBE INLET
Tt_out          = 110.3+273.15          #TUBE OUTLET
P_tube          = 1e6
TubeFluid       = Chemical('H2O',T=(Tt_out+Tt_in)/2,P=P_tube)

def deltaT(Ts_in,Ts_out,Tt_in,Tt_out,flow): #flow=co-->co-current flow=counter-->counter-current
    
    if flow=='co':
        #co-current
        deltaTML    = ((Ts_in-Tt_in)-(Ts_out-Tt_out))/np.log((Ts_in-Tt_in)/(Ts_out-Tt_out))
        r           = (Ts_out-Tt_out)/(Ts_in-Tt_in) #parameter for Tcal
    
    elif flow=='counter': 
        #counter-current
        deltaTML    = ((Ts_in-Tt_out)-(Ts_out-Tt_in))/(np.log((Ts_in-Tt_out)/(Ts_out-Tt_in)))
        r           = (Ts_out-Tt_in)/(Ts_in-Tt_out) #parameter for Tcal
        
    R = (Ts_in-Ts_out)/(Tt_out-Tt_in) #parameters for Correction Factor
    S = (Tt_out-Tt_in)/(Ts_in-Tt_in)  #parameters for Correction Factor
    print('R =',R)
    print('S =',S)
    # print('Insert Correction Factor:')
    # link to Perry Correction Factors
    # decomment line to open browser
    # webbrowser.open('https://www.slideshare.net/claumontero93/perrys-chemical-engineers-handbook-7ma-ed-chap-11')
    # CorrectionFactor=float(input())
    CorrectionFactor = 0.85 #insert num grafico C 3 shell passes 6 or more tube passes
    deltaTeff = deltaTML*CorrectionFactor
    return deltaTML,deltaTeff,r


deltaT      = deltaT(Ts_in,Ts_out,Tt_in,Tt_out,'counter')
deltaT_ML   = deltaT[0]
print('\u0394T ML =','%.2f' %deltaT_ML)

deltaT_eff  = deltaT[1]
print('\u0394T eff =','%.2f' %deltaT_eff)

plt.figure(1)
plt.title('Temperature Profile')
plt.grid('on')
plt.xlim(0,1)
plt.plot([0,1],[Ts_in-273.15,Ts_out-273.15],color='r',label='Shell Fluid') #Process Fluid
plt.plot([0,1],[Tt_out-273.15,Tt_in-273.15],color='b',label='Tube Fluid') #Service Fluid
plt.ylabel('°C')
plt.xlabel('Normalized Tube Length')
plt.legend()

# Q (HEAT EXCHANGED)
Qnom        = MShellFluid*ShellFluid.Cpg*(Ts_in-Ts_out)/1000 #kW
MTubeFluid  = Qnom/(TubeFluid.Cpl*(Tt_out-Tt_in))*1000 #kg/s

U       = np.zeros(1000)
h_i     = np.zeros(1000)
h_e     = np.zeros(1000)
h_ie    = np.zeros(1000)

# U[0]=input('U first try = ') #W/m2*K (iteration start)
U[0]=600
h_ie[0]=100         # first try
h_e[0]=110         # first try
h_i[0]=10000        # first try

#Caloric Temperature
for i in range(1000):
    Kc  = (h_e[i]-h_ie[i])/(h_ie[i])               # ci sarà un processo iterativo Uc=coeffieciente scambio lato freddo Uh coefficiente di scambio lato caldo
    r   = deltaT[2]                                 # deltaT lato freddo/ deltaT lato caldo deve essere compreso tra 0.01 e 10 (serve per calcolare T calorica)
    Fc  = (1/Kc+(r/(r-1)))/(1+(np.log(Kc+1))/np.log(r))-1/Kc
    
    #T caloriche
    Ts_cal = Ts_out+Fc*(Ts_in-Ts_out)
    Tt_cal = Tt_in+Fc*(Tt_out-Tt_in)
    
    # # CALCOLO FLUSSO DI MASSA E NUMERO DI REYNOLDS LATO TUBI E MANTELLO
    SafetyFactor    = 1.07
    Qscambio        = Qnom*SafetyFactor #kW
    Ascambio        = 1000*Qscambio/(U[i]*deltaT_eff) #m2
    
    #si sceglie la geometria dei tubi ed il materiale
    De          = DN15i.d*1e-3 #m diametro esterno
    spessore    = DN15i.t[1]*1e-3 
    Di          = De-2*spessore #m diametro interno 
    
    #vediamo a che velocità passa il fluido all'interno dei tubi
    Tube_Length = 6 #m
    Ntubi       = np.round(Ascambio/(np.pi*De*Tube_Length)+4,0)
    Ascambio    = Ntubi*np.pi*Tube_Length*De
    Tube_passes = 12
    Apassaggiotubi=(np.pi*Di**2)/4*(Ntubi/Tube_passes) #m2 sezione di passaggio fluido lato tubi
    Shell_passes = 3

    
    TubeFluid       = Chemical('H2O',T=Tt_cal,P=P_tube) # kg/m3 per Serv Fluid
    Vol_TubeFluid   = MTubeFluid/TubeFluid.rhol #m3/s portata volumica fluido lato tubi
    v_TubeFluid     = Vol_TubeFluid/Apassaggiotubi #m/s
    Re_tube         = TubeFluid.rhol*v_TubeFluid*Di/TubeFluid.mul
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
     
    def Bundle_clearance(method):
        # Coulson 6 pag. 636
        #Possible entries:
            # Fixed and U tube
            # Outside Packed Head
            # Split Ring Floating Head
            # Pull-Through Floating Head
        if method=='Fixed and U tube':
            clearance=10*Bundle_diameter+8.023
        elif method=='Outside Packed Head':
            clearance=0.376*100
        elif method=='Split Ring Floating Head':
            clearance=27.56*Bundle_diameter+44.66
        elif method=='Pull-Through Floating Head':
            clearance=9.534*Bundle_diameter+85.82
        return np.around(clearance,0)/1000 #m
         
    Shell_diameter  = Bundle_diameter+Bundle_clearance('Split Ring Floating Head')
    FormFactor      = Tube_Length/Bundle_diameter
    
    ShellFluid          = Mixture(['H2','CO2'],Vfgs=[0.75, 0.25],T=Ts_cal,P=P_Shell)  
    Vol_ShellFluid      = MShellFluid/ShellFluid.rhog                             #m3/s    
    N_Baffle            = 1
    Baffle_spacing       = Tube_Length/(1+N_Baffle)
    
    A_passaggio         = (1/Shell_passes)*Shell_diameter*(Pitch-De)*Baffle_spacing/Pitch  #area di passaggio nel punto piu largo (meta altezza)
    Ntubes_centralrow   = Bundle_diameter/Pitch
    v_ShellFluid        = Vol_ShellFluid/(A_passaggio)                         #m/s
       
    As                  = np.pi*De*1                                    #m2 area di scambio per un metro di tubo
    WettedPerimeter     = 2*(1+De)                                      #m perimetro bagnato per metro di tubo
    Deq                 = (1.1/De)*(Pitch**2-0.917*De**2)               #Coulson 6 pag 674
    Re_Shell            = ShellFluid.rhog*v_ShellFluid*Deq/ShellFluid.mug
    
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
    h_i[i]  = 4200*(1.35+0.02*(Tt_cal-273.15))*v_TubeFluid**0.8/((Di*1000)**0.2) # W/m2*K (diametro tubo interno in mm e temperatura in C)
    T_wall  = Tt_in+(h_e[i]/(h_e[i]+h_ie[i]))*(Tt_out-Tt_in)
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
    
    j_hmantello15 = lambda x: 0.8584*x**(-0.5942)+0.00279 #15% baffle cut  # COULSON PAG.672 (verificato) NON PRECISISSIMO SOPRA Re = 1e5
    j_hmantello25 = lambda x: 0.7028*x**(-0.5814)+0.002585 #25% baffle cut
    j_hmantello35 = lambda x: 0.6139*x**(-0.5692)+0.001789 #35% baffle cut
    j_hmantello45 = lambda x: 0.5885*x**(-0.5733)+0.001621 #45% baffle cut
    
    h_e[i]  = j_hmantello45(Re_Shell)*(ShellFluid.k/Deq)*Re_Shell*ShellFluid.Pr**(1/3)*(ShellFluid.mu/Mixture(['H2','CO2'],Vfgs=[0.75, 0.25],T=T_wall,P=P_Shell).mug)**(0.14)
    
    #CONDUCIBILITA' TERMICA METALLI k (W/m°C)
    k = 22  #Acciaio inox (20% Cr)
    # k = 26  #Bronzo (75% Cu, 25% Sn)
    # k = 36  #Acciaio al carbonio (1,5% C)
    # k = 54  #Acciaio al carbonio (0,5% C)
    # k = 73  #Ferro
    # k = 111  #Ottone (70% Cu, 30% Zn)
    # k = 204 #Alluminio
    # k = 386 #Rame
    
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
    print('V fluid shell side =','%.2f' %v_ShellFluid,'m/s')
    print('Form factor =','%.1f' %FormFactor)
    print('Tube Length =','%.1f' %Tube_Length, 'm')
    print('Shell Diameter =', '%.3f' %Shell_diameter, 'm')
else:
    raise ValueError('Thermodynamic Design is NOT verified')

# Perdite di carico
# Lato tubi
f_dist          = 0.0791*Re_tube**-0.25
deltaPdist_tubi = 4*f_dist*(Tube_Length*Tube_passes/Di)*(0.5*TubeFluid.rhol*v_TubeFluid**2) #Pa

D_bocchellotubi = DN100i.d*1e-3 -2*DN100i.t[1]*1e-3 #m
K_curva         = 1.5
v_TubeFluid_nozzle    = Vol_TubeFluid/(0.25*np.pi*D_bocchellotubi**2) #m/s
K_in            = 0.5*(1-(0.25*np.pi*Di**2)/(0.25*np.pi*D_bocchellotubi**2))
K_out           = 0.5*v_TubeFluid**2*(1-(0.25*np.pi*Di**2)/(0.25*np.pi*D_bocchellotubi**2))**2

deltaPconc_tubi = 0.5*K_curva*(Tube_passes-1)*(TubeFluid.rhol*v_TubeFluid**2)+0.5*(K_in+K_out)*(TubeFluid.rhol*v_TubeFluid_nozzle**2) #Pa
deltaPtot_tubi  = deltaPdist_tubi+deltaPconc_tubi

# if deltaPtot_tubi>100000:  # check perdite di carico massime consentite! Check che con T out fluido e Pin-perdite di carico totali non ci sia flash nel tubo!
    # raise ValueError('Pressure drop tube side is too high!')
    
# # #perdite di carico
# # #lato mantello

j_fmantello15 = lambda x: 37.95*x**(-1.021)+0.05488     #15% baffle cut COULSON PAG.673 (verificato)
j_fmantello25 = lambda x: 26.91*x**(-1.046)+0.04189     #25% baffle cut
j_fmantello35 = lambda x: 20.9*x**(-1.019)+0.03504      #35% baffle cut
j_fmantello45 = lambda x: 15.18*x**(-0.9707)+0.02376    #45% baffle cut

#la scelta del baffle cut deve essere coerente con quella fatta precedentemente per il jH mantello
deltaPdist_Shell = Shell_passes*8*j_fmantello45(Re_Shell)*(Shell_diameter/De)*(Tube_Length/Baffle_spacing)*(0.5*ShellFluid.rho*v_ShellFluid**2)

Kin = 1
Kout = 0.5
K90 = 1.5

Dbocch_Shell        = (DN450i.d-DN450i.t[2])/1000 #dovrei annotare le dimensioni dello scambiatore per capire quale bocchello effettivamente potrebbe starci
Abocch_Shell        = 0.25*np.pi*Dbocch_Shell**2
Velbocch_Shell      = (v_ShellFluid*A_passaggio)/Abocch_Shell 
deltaPconc_Shell    = (Kin+Kout)*0.5*ShellFluid.rho*Velbocch_Shell**2 + (2*Kout)*0.5*ShellFluid.rho*v_ShellFluid**2
deltaPtot_Shell     = deltaPconc_Shell+deltaPdist_Shell

if deltaPtot_Shell>80000:
    raise ValueError('Pressure Shell side is too high!')

print('---------------------------------------------------------------------------')
print('Pressure Drop')
print('Tube side Pressure Drop =', '%.2f' %(deltaPtot_tubi/1000),'kPa')
print('Shell side Pressure Drop =', '%.2f' %(deltaPtot_Shell/1000),'kPa')