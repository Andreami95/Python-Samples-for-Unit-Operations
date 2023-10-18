
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from thermo.chemical import Chemical, Mixture, VolumeLiquidMixture
from scipy import integrate
from Pipes_Steel import *
from Pipes_Inox import * 

plt.close('all')

#Shell
Ts_in = 85+273.15                                                   #K
Ts_out = 0 +273.15                                               #K
P_Shell = 1.9e5                                                           #Pa Pressure Shell Fluid
C1 = Chemical('CH3OH',T=Ts_in,P=P_Shell)                                     #Condensable compound
C2 = Chemical('H2O',T=Ts_in,P=P_Shell)                                     #Condensable compound
I = Chemical('CO2',T=Ts_in,P=P_Shell)                                       #Incondensable compound

#Tube
TubeFluid = Mixture(['Water','NaCl'],Vfls=[.85, .15], T=273.16, P=3E5)    #Tube Fluid
Tt_out = 273.15-10                                                          #K T output Tube Fluid
Tt_in = 273.15-20                                                    #K T input Tube Fluid
deltat = Tt_out-Tt_in

#TUBE GEOMETRY
De = DN25i.d/1000                                                   #m
s = DN25i.t[0]/1000                                                 #m
Di = De-2*s                                                         #m
Pitch = 1.25
D_pitch = Pitch*De                                                  #m

#FEED
kmol_ShellFluid = (250.6)                                                  #kmol/h
C_frac = 0.74                                                      #molar fraction methanol (condensable)
I_frac = np.around(1. -C_frac,1)                                   #molar fraction incondensable
Feed_C = (kmol_ShellFluid * C_frac)                                        #Condensable molar flow kmol/h
Feed_I = (kmol_ShellFluid * I_frac)                                        #Incondensable molar flow kmol/h

T = np.array([Ts_in,345,334.5,308,Ts_out])                              # Array of Temperatures to ideally divide Condenser
num = len(T)
T = np.append(T,Ts_out)

dL  = np.zeros(num)
Vc  = np.zeros(num)
Lc  = np.zeros(num)        
k   = np.zeros(num)
dQ  = np.zeros(num)
t   = np.zeros(num)
dQ  = np.zeros(num)
dT  = np.zeros(num)
Hie = np.zeros(num)
He  = np.zeros(num)
Hei = np.zeros(num)
Nu  = np.zeros(num)
yI  = np.zeros(num)
tw  = np.zeros(num)
v_ShellFluid = np.zeros(num)
A_s = np.zeros(num)
j_h = np.zeros(num)
j_H = np.zeros(num)
Hi  = np.zeros(num)
G   = np.zeros(num)

Heat_Liq_C  = np.zeros(num)
Heat_Vap    = np.zeros(num)
Heat_Vap_C  = np.zeros(num)
H_vap_C     = np.zeros(num)
T_wall      = np.zeros(num)
deltaT      = np.zeros(num)
deltaT_LN   = np.zeros(num)
Upulito     = np.zeros(num)
Qprova_     = np.zeros(num)
U_dT_LN     = np.zeros(num)
dQdT        = np.zeros(num)
T_film      = np.zeros(num)

Re_t = np.zeros(num)
Re_m = np.zeros(num)
Pr_m = np.zeros(num)

################## SECTION 0 ###################

C_i     = Chemical('CH3OH',T=T[0],P=P_Shell)
C_i1    = Chemical('CH3OH',T=T[1],P=P_Shell)

I_i     = Chemical('CO2',T=T[0],P=P_Shell)
I_i1    = Chemical('CO2',T=T[1],P=P_Shell)

k[0]    = C_i.Psat/P_Shell
Vc[0]   = Feed_C
#Condensable
Cpg_av_C        = (C_i.Cpgm + C_i1.Cpgm)/2                          #kJ/kmol K Gas Average Specific Heat
H_vap_C[0]      = C_i.Hvapm                                         #kJ/kmol Heat condensation
Heat_Vap_C[0]   = Vc[0]*(Cpg_av_C*(T[0]-T[1]) + H_vap_C[0])/3600    #kW (sens+lat) Q condensation 
Cpl_av_C        = (C_i.Cplm + C_i1.Cplm)/2                          #kJ/kmol K Liquid Average Specific Heat
Lc[0]           = 0                                                 #kmol/h Liquid Condensed Flow
Heat_Liq_C[0]   = 0                                                 #kW
 
# Incondensable
Cpg_av_I    = (I_i.Cpgm + I_i1.Cpgm)/2                              #kJ/kmol K Gas Average Specific Heat
Heat_Vap_I  = Feed_I*(Cpg_av_I*(T[0]-T[1]))/3600                    #kW Q 
Heat_Vap[0] = Heat_Vap_C[0] + Heat_Vap_I                            #kW Total Heat Vap
dQ[0]       = 0                                                     #kW Heat exchanged in SECTION 0
Q           = dQ[0]

#### ITERATION FROM SECTION 1 TO END SECTION (num)

for i in range(1,num):
    
    C_i     = Chemical('CH3OH',T=T[i],P=P_Shell)
    C_i1    = Chemical('CH3OH',T=T[i+1],P=P_Shell)
    I_i     = Chemical('CO2',T=T[i],P=P_Shell)
    I_i1    = Chemical('CO2',T=T[i+1],P=P_Shell)
    
    k[i]    = C_i.Psat/P_Shell
    x       = np.zeros(100)
    x[0]    = 1                                                     # guess Value of Vapor/Liquid in i-th section
    
    
    for j in range(100-1):                                          #Iteration on x to find the true Value
        
        dL[i] = (Vc[i-1]-Lc[i-1]*k[i]*x[j])/(1+x[j]*k[i])           #kmol/h condensed Liquid in i-th section
        Vc[i] = Vc[i-1]-dL[i]                                       #kmol/h of condensable Vapor in i-th section 
        
        #Condensable
        Cpg_av_C        = (C_i.Cpgm + C_i1.Cpgm)/2                          #kJ/kmol K Gas Average Specific Heat
        H_vap_C[i]      = C_i.Hvapm                                         #kJ/kmol Heat condensation
        Heat_Vap_C[i]   = Vc[i]*(Cpg_av_C*(T[i]-T[i+1]) + H_vap_C[i])/3600  #kW (sens+lat) Q condensation 
        Cpl_av_C        = (C_i.Cplm + C_i1.Cplm)/2                          #kJ/kmol K Liquid Average Specific Heat
        Lc[i]           = Lc[i-1]+dL[i]                                     #kmol/h Liquid Condensed Flow
        Heat_Liq_C[i]   = Lc[i]*Cpl_av_C*(T[i]-T[i+1])/3600                 #kW

        #Incondensable
        Cpg_av_I        = (I_i.Cpgm + I_i1.Cpgm)/2                         #kJ/kmol K Gas Average Specific Heat
        Heat_Vap_I      = Feed_I*(Cpg_av_I*(T[i]-T[i+1]))/3600             #kW
        yI[i]           = Feed_I/(Vc[i]+Feed_I)                            #molar fraction incondensables in Vapor

        Heat_Vap[i] = Heat_Vap_C[i] + Heat_Vap_I                            #kW Total Heat Vapor
        dQ[i] = Heat_Vap[i-1] + Heat_Liq_C[i-1] - Vc[i] * H_vap_C[i]/3600   #kW Heat exchanged in i-th section
                  
        Vc[i] = Vc[i-1] - dL[i]  
        x[j+1] = (Feed_I+Vc[i])/Lc[i]
        if np.abs(x[j]-x[j+1])<0.001:
            break
    dQdT[i]=dQ[i]/(T[i-1]-T[i])

Q = sum(dQ) #kW

t[0] = Tt_out
for i in range(1,num):
    t[i]    = t[i-1] - dQ[i]/Q*deltat                   #T cold utility per ogni sezione
    dT[i]   = 0.5*((T[i-1]- t[i-1])+(T[i]- t[i]))       #deltaT medio tra due sezioni consecutive

dQdT_ShellTube  = dQ[1:]/dT[1:]
dT_ML           = Q/sum(dQdT_ShellTube)
D_eq_s          = 4*((1.732/4)*D_pitch**2 - (np.pi/8)*De**2)/(0.5*np.pi*De) #diametro equivalente triangolare

MTubeFluid      = 1000*Q/(TubeFluid.Cpl*deltat) #kg/s
Vol_TubeFluid   = MTubeFluid/TubeFluid.rhol #m3/s

ShellFluid     = Mixture(['CH3OH','CO2'],Vfgs=[C_frac, I_frac], T=Ts_in, P=P_Shell)
G[0]    = kmol_ShellFluid*ShellFluid.MW/3600 #kg/s
Re_m[0] = G[0]*(D_eq_s)/ShellFluid.mug
Pr_m[0] = ShellFluid.Prg
v_ShellFluid[0]  = Re_m[0]*ShellFluid.mu/(ShellFluid.rhog*(D_eq_s))
Nu[0]   = 0.36*(Re_m[0])**0.55*(Pr_m[0])**0.33

U_tentativo = 350 #W/m2 K
Safety_Factor = 1.1
Ascambio = Safety_Factor*1000*Q/(U_tentativo*dT_ML) #m2
Tube_Length = 6 #m
Ntubi = np.around(Ascambio/(np.pi*(De)*Tube_Length)+4,0)


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
        elif Npasses==8:
            K_1=0.0365
            n_1=2.675
        D_b=(De)*(Ntubi/K_1)**(1/n_1)
    
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
        elif Npasses==8:
            K_1=0.0331
            n_1=2.643
        D_b=(De)*(Ntubi/K_1)**(1/n_1)
    

    return D_b

Tube_passes = 2
Bundle_diameter = np.around(D_Bundle(Tube_passes,Ntubi,TubeConfig),3) #m

Shell_diameter = Bundle_diameter+2*De
FormFactor = Tube_Length/Shell_diameter

l_spacing = 0.3*Shell_diameter #m 0.3 e 0.5 di Shell diameter best choice
A_pass_Shell = l_spacing*Bundle_diameter*(D_pitch-De)/D_pitch
v_ShellFluid[0] = (G[0]/ShellFluid.rhog)/A_pass_Shell
Re_m[0] = ShellFluid.rhog*v_ShellFluid[0]*(D_eq_s/1000)/ShellFluid.mug
Pr_m[0] = ShellFluid.Prg
Nu[0] = 0.36*(Re_m[0])**0.55*(Pr_m[0])**0.33


He[0] = Nu[0]*ShellFluid.k/(D_eq_s) #W/m2k esterno
Hei[0] = He[0]*De/Di

A_pass_tubi = 0.25*(Ntubi/Tube_passes)*np.pi*(Di)**2
vtubi = Vol_TubeFluid/A_pass_tubi
Re_t[0] = TubeFluid.rhol*vtubi*(Di)/TubeFluid.mul
T_wall[0] = 310 #K

jH_Kern = lambda Re: -2.566e-20*Re**4 + 2.655e-14*Re**3 - 9.973e-09*Re**2 + 0.003371*Re +2.667
j_H[0]  = jH_Kern(Re_t[0])
Hi[0]   = j_H[0]*TubeFluid.Pr**0.33*(TubeFluid.k/(Di))*(TubeFluid.mu/Mixture(['Water','NaCl'],Vfls=[.85, .15], T=T_wall[0], P=1E5).mul)**(0.14)
Hie[0]  = Hi[0]*Di/De

C_i     = Chemical('CH3OH',T=T[0],P=P_Shell)
C_i1    = Chemical('CH3OH',T=T[1],P=P_Shell)
I_i     = Chemical('CO2',T=T[0],P=P_Shell)
I_i1    = Chemical('CO2',T=T[1],P=P_Shell)

Qprova  = lambda T_wall_: Hie[0]*(T_wall_-t[0]) #W/m2
Pg_ML   = lambda T_wall_: (P_Shell-Chemical('CH3OH',T=0.5*(T[0]+T_wall_),P=P_Shell).Psat)/(np.log(P_Shell/Chemical('CH3OH',T=0.5*(T[0]+T_wall_),P=P_Shell).Psat)) #PA
Diffusivity = ((0.001*T[0]**1.75)*(1/(1/C_i.MW+1/I_i.MW))**0.5)/(P_Shell*(29.9**0.33 + 16.9**0.33)**2) #Fuller Schettler Giddings 5-48 Perry
kg      = He[j]/(ShellFluid.Cpg*ShellFluid.rhog)*ShellFluid.Prg**(2/3)*(ShellFluid.mug/(ShellFluid.rhog*Diffusivity))**(-2/3) #mass transfer coefficient 
Qprova2 = lambda T_wall_: He[0]*(T[0]-0.5*(T[0]+T_wall_))+(kg*ShellFluid.rhog/(Pg_ML(T_wall_)))*C_i.Hvap*(Chemical('CH3OH',T=T[0],P=P_Shell).Psat- Chemical('CH3OH',T=0.5*(T[0]+T_wall_),P=P_Shell).Psat)
diffQ   = lambda T_wall_: Qprova(T_wall_)-Qprova2(T_wall_)

T_wall[0]   = fsolve(diffQ,310)
tw[0]       = T_wall[0]

for i in range(1,100):  #CALCOLO T WALL
    j_H     = jH_Kern(Re_t[0])
    Hi[0]   = j_H*TubeFluid.Pr**0.33*(TubeFluid.k/(Di))*(TubeFluid.mu/Mixture(['Water','NaCl'],Vfls=[.85, .15], T=273.15, P=3E5).mul)**(0.14)
    Hie[0]  = Hi[0]*Di/De
    tw[i]   = fsolve(diffQ,330)
    if np.abs(tw[i-1]-tw[i])<1e-6:
        T_wall[0]=tw[i]
        break

Qprova_[0]  = Qprova(T_wall[0])
T_film[0]   = 0.5*(T[0]+T_wall[0])
deltaT[0]   = T[0]-t[0]
Upulito[0]  = Qprova2(T_wall[0])/deltaT[0]


# #-----------------------------------------------------------------------------------------------------------

for j in range(1,num):
    ShellFluid     = Mixture(['CH3OH','CO2'],zs=[1-yI[j], yI[j]], T=T[j], P=P_Shell)
    G[j]    = Vc[j]*ShellFluid.MW/3600 #kg/s
    v_ShellFluid[j]  = (G[j]/ShellFluid.rhog)/A_pass_Shell
    Re_m[j] = ShellFluid.rhog*v_ShellFluid[j]*(D_eq_s)/ShellFluid.mug
    Pr_m[j] = ShellFluid.Prg
    Nu[j]   = 0.36*(Re_m[j])**0.55*(Pr_m[j])**0.33

    He[j]   = Nu[j]*ShellFluid.kg/(D_eq_s) #W/m2k esterno
    Hei[j]  = He[j]*De/Di        

    T_wall[j] = 300 #K
    
    Re_t[j] = TubeFluid.rhol*vtubi*(Di)/TubeFluid.mul
    j_H     = jH_Kern(Re_t[j])
    Hi[j]   = j_H*TubeFluid.Pr**0.33*(TubeFluid.k/(Di))*(TubeFluid.mu/Mixture(['Water','NaCl'],Vfls=[.85, .15], T=T_wall[j], P=3E5).mul)**(0.14) #da kern
    Hie[j]  = Hi[j]*Di/De

    C_i     = Chemical('CH3OH',T=T[j],P=P_Shell)
    C_i1    = Chemical('CH3OH',T=T[j+1],P=P_Shell)
    I_i     = Chemical('CO2',T=T[j],P=P_Shell)
    I_i1    = Chemical('CO2',T=T[j+1],P=P_Shell)


    # #######RICONTROLLA UNITA DI MISURA!!!!
    Qprova      = lambda T_wall_: Hie[j]*(T_wall_-t[j]) #W/m2
    Pg_ML       = lambda T_wall_: (P_Shell-Chemical('CH3OH',T=0.5*(T[j]+T_wall_),P=P_Shell).Psat)/(np.log(P_Shell/Chemical('CH3OH',T=0.5*(T[j]+T_wall_),P=P_Shell).Psat)) #PA
    Diffusivity = ((0.001*T[j]**1.75)*(1/(1/C_i.MW+1/I_i.MW))**0.5)/(P_Shell*(29.9**0.33 + 16.9**0.33)**2) #Fuller Schettler Giddings 5-48 Perry
    kg          = He[j]/(ShellFluid.Cpg*ShellFluid.rhog)*ShellFluid.Prg**(2/3)*(ShellFluid.mug/(ShellFluid.rhog*Diffusivity))**(-2/3) #mass transfer coefficient 
    Qprova2     = lambda T_wall_: He[j]*(T[j]-0.5*(T[j]+T_wall_))+(kg*ShellFluid.rhog/(Pg_ML(T_wall_)))*C_i.Hvap*(Chemical('CH3OH',T=T[j],P=P_Shell).Psat- Chemical('CH3OH',T=0.5*(T[j]+T_wall_),P=P_Shell).Psat)
    diffQ       = lambda T_wall_: Qprova(T_wall_)-Qprova2(T_wall_)
    T_wall[j]   = fsolve(diffQ,300)
    tw[j]       = T_wall[j]
    
    for k in range(1,100):
        # print(i)
        jH_Kern = lambda Re: -2.566e-20*Re**4 + 2.655e-14*Re**3 - 9.973e-09*Re**2 + 0.003371*Re +2.667
        j_H     = jH_Kern(Re_t[j])
        Hi[k]   = j_H*TubeFluid.Pr**0.33*(TubeFluid.k/(Di))
        Hie[k]  = Hi[k]*Di/De
        tw[k]   = fsolve(diffQ,310)
        if np.abs(tw[k-1]-tw[k])<1e-6:
            T_wall[j]=tw[k]
            break
####---------------------------------------------------

    Qprova_[j] = Qprova(T_wall[j])
    T_film[j] = 0.5*(T[j]+T_wall[j])
    deltaT[j] = T[j]-t[j]
    deltaT_LN[j] = (deltaT[j-1]-deltaT[j])/np.log(deltaT[j-1]/deltaT[j])
    U_dT_LN[j] = (Qprova_[j]-Qprova_[j-1])/np.log(Qprova_[j]/Qprova_[j-1]) #W/m2
    Upulito[j] = (U_dT_LN[j]/deltaT_LN[j])
    A_s[j] = 1000*dQ[j]/U_dT_LN[j]

A_tot_scambio = sum(A_s) #m2 fare check con area effettiva (dovrebbe essere sempre minore)
R_sporco = 0.001
U_tot = ((sum(Upulito[1:]*A_s[1:])/A_tot_scambio)**-1 + R_sporco + (s/1000)/17)**-1

Verifica = dT_ML*Ascambio*U_tot/(1000*Q)        #Verifica termica
if Verifica>1:
    print('Thermodynamic Design is verified')
    print('Q eff/Q ric =','%.3f' %Verifica)
    print('Q eff = ''%.2f' %(dT_ML*Ascambio*U_tot/1e6),'MW')
    
# if FormFactor<6:
#     #CHECK FORM-FACTOR HEAT EXCHANGER
#     print ('Form Factor =','%.1f' %FormFactor)
#     raise ValueError('Form Factor is too low!')
# elif FormFactor>10:
#     print ('Form Factor =','%.1f' %FormFactor)
#     raise ValueError('Form Factor is too high!')


##### Perdite di carico LATO TUBI 
f_dist = 0.0791*Re_t[0]**-0.25
deltaPdist_tubi = 4*f_dist*(Tube_Length*(Tube_passes-1)/(Di))*(0.5*TubeFluid.rhol*vtubi**2) #Pa

#Concentrate
#Si devono aggiungere i fattori K per le perdite concentrate
D_bocchellotubi = (DN250i.d-2*DN250i.t[2])/1000 #m
K_curva = 1.5

#Dobbiamo calcolarci la velocità del fluido attraverso i bocchelli
velbocchTube = Vol_TubeFluid/(0.25*np.pi*D_bocchellotubi**2) #m/s
K_in = 0.5*(1-(0.25*np.pi*(Di)**2)/(0.25*np.pi*D_bocchellotubi**2))
K_out = 0.5*vtubi**2*(1-(0.25*np.pi*(Di)**2)/(0.25*np.pi*D_bocchellotubi**2))**2

deltaPconc_tubi = 0.5*K_curva*(Tube_passes-1)*(TubeFluid.rhol*vtubi**2)+0.5*(K_in+K_out)*(TubeFluid.rhol*velbocchTube**2) #Pa
deltaPtot_tubi = deltaPdist_tubi+deltaPconc_tubi

j_fmantello15 = lambda x: 37.95*x**(-1.021)+0.05488     #15% baffle cut PARTE LIBERA DI PASSAGGIO
j_fmantello25 = lambda x: 26.91*x**(-1.046)+0.04189     #25% baffle cut
j_fmantello35 = lambda x: 20.9*x**(-1.019)+0.03504      #35% baffle cut
j_fmantello45 = lambda x: 15.18*x**(-0.9707)+0.02376    #45% baffle cut

deltaPdist_Shell=(8*j_fmantello45(sum(Re_m)/num)*(Shell_diameter/(De))*(Tube_Length/l_spacing)*(0.5*ShellFluid.rhog*(sum(v_ShellFluid)/num)**2))

Kin = 1
Kout = 0.5
K90 = 1.5

#INLET
Gas_in         = Mixture(['CH3OH','CO2'],Vfgs=[C_frac, I_frac], T=Ts_in, P=P_Shell)
Dbocch_Shellin = (DN400i.d-DN400i.t[0]*2)/1000 #dovrei annotare le dimensioni dello scambiatore per capire quale bocchello effettivamente potrebbe starci
Abocch_Shellin = 0.25*np.pi*Dbocch_Shellin**2
Gas_inlet      = (1000/3600)*kmol_ShellFluid/Gas_in.rhogm
Velbocch_Shellin =(Gas_inlet)/Abocch_Shellin 
deltaPconc_Shellin = Kin*1.5*Gas_in.rhog*Velbocch_Shellin**2

#OUTLET
Gas_out         = Mixture(['CH3OH','CO2'],Vfgs=[1-yI[num-1], yI[num-1]], T=Ts_out, P=P_Shell)
Dbocch_Shellout = (DN400i.d-DN400i.t[0]*2)/1000 #dovrei annotare le dimensioni dello scambiatore per capire quale bocchello effettivamente potrebbe starci
Abocch_Shellout = 0.25*np.pi*Dbocch_Shellout**2
Gas_outlet      = (1000/3600)*kmol_ShellFluid/Gas_out.rhogm
Velbocch_Shellout = (Gas_outlet)/Abocch_Shellout 
deltaPconc_Shellout = Kout*1.5*Gas_out.rhog*Velbocch_Shellout**2

deltaPconc_Shell    = deltaPconc_Shellout + deltaPconc_Shellin
deltaPtot_Shell     = deltaPconc_Shell+deltaPdist_Shell

Di_Liq_Bocc = DN65.d-2*DN65.t[0]
vLiqOutlet = (1000*(Feed_C-Vc[num-1])/Chemical('methanol',T=273.15-10.82,P=1.2e5).rholm/3600)/(np.pi*0.25*(Di_Liq_Bocc/1000)**2) #m/s

# Tempo Residenza
Hmax            = 0.5*(Shell_diameter - Bundle_diameter) #livello max liquido
alpha           = np.arccos(1-(2*Hmax/Shell_diameter))
Volmax          = Tube_Length*(alpha*(Shell_diameter/2)**2-(Shell_diameter/2-Hmax)*(Hmax*Shell_diameter-Hmax**2)**0.5) #m3
ResidenceTime   = Volmax/(1000*(Feed_C-Vc[num-1])/Chemical('methanol',T=273.15-10.82,P=1.2e5).rholm/3600)

print( 'N tubi = ', '%.0f' %Ntubi)
print( 'Ascambio = ', '%.1f' %Ascambio,'m2')
print('Diameter Shell =','%.3f' %Shell_diameter, 'm')
print('V fluid Tube side =','%.2f' %vtubi,'m/s')
print('Average V fluid Shell side =','%.2f' %(np.average(v_ShellFluid)),'m/s')
print('U =','%.1f' %U_tot,'W/m2 K')
print('Qscambiato su Q ric =','%.3f' %Verifica)
print('Form Factor =','%.1f' %FormFactor)
print('---------------------------------------------------------------------------')
print('Pressure Drop')
print('Tube side Pressure Drop =','%.1f' %(deltaPtot_tubi/1000),'kPa')
print('Shell side Pressure Drop =' ,'%.1f'%(deltaPtot_Shell/1000),'kPa')

HenryWater = lambda T: 1/(np.exp(-159.854 + 8741.68/T + 21.6694*np.log(T) -1.10261e-3*T))

H = HenryWater(273)

xi = 0.2*(1.3)/H + 8*0.8*(1.3)/H