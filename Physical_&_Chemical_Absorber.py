# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 13:12:08 2022

@author: Andrea Milazzo
"""

from thermo.chemical import Chemical
from thermo import Mixture,VolumeLiquidMixture,VolumeGasMixture
import numpy as np
import sys
from scipy.optimize import fsolve,root
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import linregress
from fluids.packed_tower import separation_demister_ElDessouky,dP_demister_wet_ElDessouky,dP_demister_dry_Setekleiv_Svendsen
plt.close('all')

ratio = lambda x:x/(1-x)    #ottengo il rapporto molare dalla frazione molare
frac = lambda x:x/(1+x)     #ottengo la frazione molare dal rapporto molare

H2O_humidity        = lambda T:(Chemical('H2O',T=T,P=101325))
y_water             = lambda T: H2O_humidity(T).Psat/101325
UA100               = lambda T: y_water(T)*H2O.MW/((1-y_water(T))*air.MW)                                       #kgH2O/kg dry air
WaterCondensed      = lambda Tin,Tout,G: G*(UA100(Tin)-UA100(Tout))                                             #G=dry air in kg >>Gp_air
HeatCondensation    = lambda Tin,Tout,G: WaterCondensed(Tin,Tout,G)*H2O_humidity(0.5*(Tin+Tout)).Hvap/1000      #kj/h se Tout>Tin allora Heat<0
HeatAbsorption      = lambda Yin,Yout: Gs*(Yin-Yout)*HCl.MW*2100/1000                                           #kj/h
delta_T             = lambda Yin,Yout,Lp: (HeatAbsorption(Yin,Yout))/(Lp*H2O.Cpl/1000)     #calcola deltaT tra top e bottom di colonna che riferito al solvente liquido
                                                                                                                                                            #facciamo l'ipotesi che il gas scambi calore ed esca alla temperatura di entrata del liquido
X_1                 = lambda Y:(Y-q1)/m1 #--->trovo X da Y perchè uso la RETTA di lavoro che si ottiene dal grafico rapporti molari, uso X in ye
ye1                 = lambda Y: Solubility30(frac(X_1(Y)))
yeprova             = lambda Y: 0.5*(Solubility30(frac(X_1(Y)))+Solubility20(frac(X_1(Y)))) 
integrale           = lambda Y :1/((frac(Y))-ye1(Y)) #perchè serve la Y per X, poi passiamo a frazioni 
integraleprova      = lambda Y :1/((frac(Y))-yeprova(Y))

X_2                 = lambda Y:(Y-q_C2)/m_C2                                                 #--->trovo X da Y perchè uso la RETTA di lavoro che si ottiene dal grafico dei rapporti molari
ye2                 = lambda Y: Solubility30(frac(X_2(Y)))
integrale2          = lambda Y :1/((frac(Y))-ye2(Y))                                  #perchè serve la Y per X, poi passiamo a frazioni 

#FUNZIONI VALIDE PER GPDC RANDOM PACKING
GPDC_flooding       = lambda x: (1.896*x**2+23.35*x+5.929)/(x**5-4.259*x**4+15.8*x**3+23.05*x**2+15.81*x+0.7963)
GPDC_4mmH2O_r       = lambda x: (0.3624*x**2+0.06882*x+0.01766)/(x**3+2.002*x**2+0.6218*x+0.1137)
GPDC_8mmH2O_r       = lambda x: (15.65*x**2+1.054*x+0.03889)/(x**5-4.736*x**4+55.42*x**3+27.72*x**2+3.416*x+0.05792)
GPDC_21mmH2O_r      = lambda x: (0.3608*x**2-0.003868*x+0.0003302)/(x**3+0.2767*x**2+0.008821*x+7.255e-05)
GPDC_42mmH2O_r      = lambda x: (0.3976*x**2+0.1328*x-0.001566)/(x**3+0.3518*x**2+0.09089*x-0.001283)
GPDC_83mmH2O_r      = lambda x: (-1.296e+04*x**3+5.028e+04*x**2+8.31e+04*x+880.6)/(x**5-3.317e+04*x**4+1.516e+05*x**3+1.171e+05*x**2+3.413e+04*x+65.82)
#FUNZIONI VALIDE PER GPDC STRUCTURED PACKING FIG.14-56 PERRY 9 e PERRY 8
GPDC_8mmH2O_s       = lambda x: -0.09869*x**7+0.8262*x**6-2.794*x**5+4.88*x**4-4.717*x**3+2.64*x**2-1.134*x+0.7424
GPDC_21mmH2O_s      = lambda x: (-16*x**5+76.11*x**4+88.12*x**3-74.58*x**2+14.49*x+0.1311)/(x**5+299.9*x**4-72.09*x**3-38.93*x**2+12.98*x+0.1054)
GPDC_42mmH2O_s      = lambda x: (2.276*x**2+4.793*x+0.0717)/(x**5-3.447*x**4+4.99*x**3+8.059*x**2+3.674*x+0.03991)
GPDC_83mmH2O_s      = lambda x: (8122*x**3-(1.027e+04)*x**2+(1.713e+04)*x+472.8)/(x**5+9996*x**4-(1.034e+04)*x**3+(1.435e+04)*x**2+(1.279e+04)*x+200.2)
GPDC_125mmH2O_s     = lambda x: (3731*x**3+1177*x**2+6384*x+224.5)/(x**5+3346*x**4+2976*x**3+7787*x**2+4665*x+90.95)

# SOLUBILITA' FRAZIONI MOLARI ##############################################################################################################
Solubility30        = lambda x: 1.494e+04*x**6  - 6069*x**5 + 1018*x**4 - 87.16*x**3 + 3.957*x**2 - 0.0883*x + 0.0007383
Solubility20        = lambda x: -1.097e+04*x**7 + 6770*x**6 - 1372*x**5 + 135.4*x**4 - 7.055*x**3 + 0.1951*x**2 - 0.002578*x + 1.253e-05

# SOLUBILITA' RAPPORTI MOLARI
Solubility30r       = lambda x:ratio(Solubility30(x))                   #input frazioni molari >> output rapporti molari
Solubility20r       = lambda x:ratio(Solubility20(x))

op_line0            = lambda x: Yout + m_min*x                          #variabile x in RAPPORTI MOLARI
op_line_maxRecycle  = lambda x: Y0 + m_max*x                            #variabile x in RAPPORTI MOLARI

retta_op1           = lambda x: m1*x+q1                                 #la variabile è in rapporti molari returns in rapporti molari
intersection_in     = lambda x: Solubility30(x)-yin              #variabile in rapporti molari returns in rapporti molari
intersection_out    = lambda x: Solubility30(x)-yout                    #variabile in frazione molare returns in frazioni molari

#DEFINIZIONE CLASSE RIEMPIMENTO LISTE
class Packing:
    def __init__(self, size, area, void, Packing_Factor,alpha,S,B,h): 
#Area m^2/m^3, voids %, Packing Factor m^-1,
# alpha=corrugation angle °, S = Channel Side [mm], B=base channel [mm], h=crimp height [mm]
        self.size=size
        self.area=area
        self.void=void
        self.Packing_Factor=Packing_Factor
        self.alpha=alpha
        self.S=S
        self.B=B
        self.h=h

## MELLAPAK ---------------------------------------------------------
#Materials: Stainless Steel, Carbon Steel, Hastelloy, monel, aluminium, copper-bronze, brass, titanium and nickel    

Mellapak_list = []

Mellapak_list.append( Packing('125y',125,0.99,33,45,37,55,24.8))
Mellapak_list.append( Packing('170y',170,0.99,39,'none','none','none','none'))
Mellapak_list.append( Packing('2y',223,0.99,46,45,21.5,33,13.8))
Mellapak_list.append( Packing('250y',250,0.98,66,45,17,24.1,11.9))
Mellapak_list.append( Packing('350y',350,0.98,75,'none','none','none','none'))
Mellapak_list.append( Packing('500y',500,0.98,112,45,8.1,9.6,6.53))
Mellapak_list.append( Packing('750y',750,0.97,'none','none','none','none','none'))
Mellapak_list.append( Packing('125x',125,0.99,16,'none','none','none','none'))
Mellapak_list.append( Packing('170x',170,0.99,20,'none','none','none','none'))
Mellapak_list.append( Packing('2x',223,0.99,23,'none','none','none','none'))
Mellapak_list.append( Packing('250x',250,0.98,26,60,17,24.1,11.9))
Mellapak_list.append( Packing('350x',350,0.98,'none','none','none','none','none'))
Mellapak_list.append( Packing('500x',500,0.98,82,'none','none','none','none'))

Mellapak_125y=Mellapak_list[0]
Mellapak_170y=Mellapak_list[1]
Mellapak_2y=Mellapak_list[2]
Mellapak_250y=Mellapak_list[3]
Mellapak_350y=Mellapak_list[4]
Mellapak_500y=Mellapak_list[5]
Mellapak_750y=Mellapak_list[6]
Mellapak_125x=Mellapak_list[7]
Mellapak_170x=Mellapak_list[8]
Mellapak_2x=Mellapak_list[9]
Mellapak_250x=Mellapak_list[10]
Mellapak_350x=Mellapak_list[11]
Mellapak_500x=Mellapak_list[12]

# MELLAPAK PLASTIC------------------------------------------------------------------
# Materials: PP, PVC-C, PVDF, TEFLON PFA, PEEK
Mellapak_Plastic_list=[]

Mellapak_Plastic_list.append (Packing('125y',125,0.99,33,45,37,55,24.8))
Mellapak_Plastic_list.append (Packing('250y',250,0.96,72,45,17,24.1,11.9))
Mellapak_Plastic_list.append (Packing('125x',125,'none','none','none','none','none','none'))
Mellapak_Plastic_list.append (Packing('250x',250,'none','none','none','none','none','none'))

Mellapak_Plastic_125y=Mellapak_Plastic_list[0]
Mellapak_Plastic_250y=Mellapak_Plastic_list[1]
Mellapak_Plastic_125x=Mellapak_Plastic_list[2]
Mellapak_Plastic_250x=Mellapak_Plastic_list[3]

# MELLAPAK PLUS ----------------------------------------------------
#Materials: Stainless Steel, Carbon Steel, Hastelloy, monel, aluminium, copper-bronze, brass, titanium and nickel    
Mellapak_Plus_list =[]

Mellapak_Plus_list.append( Packing('202y','none',0.99,'none','none','none','none','none'))
Mellapak_Plus_list.append( Packing('252y',250,0.98,39,45,17,24.1,11.9))
Mellapak_Plus_list.append( Packing('352y','none',0.98,'none','none','none','none','none'))
Mellapak_Plus_list.append( Packing('452y',350,0.98,69,'none','none','none','none'))
Mellapak_Plus_list.append( Packing('752y',500,0.98,131,'none','none','none','none'))

Mellapak_Plus_202y=Mellapak_Plus_list[0]
Mellapak_Plus_252y=Mellapak_Plus_list[1]
Mellapak_Plus_352y=Mellapak_Plus_list[2]
Mellapak_Plus_452y=Mellapak_Plus_list[3]
Mellapak_Plus_752y=Mellapak_Plus_list[4]

class RandomPacking:
    def __init__(self, size,bulk_density, area, Packing_Factor): 
        self.size=size
        self.bulk_density=bulk_density
        self.area=area
        self.Packing_Factor=Packing_Factor

# RASCHIG RINGS ---------------------------------------------------------
Raschig_ceramic_list = []

Raschig_ceramic_list.append( RandomPacking(13,881,368,2100))
Raschig_ceramic_list.append( RandomPacking(25,673,190,525))
Raschig_ceramic_list.append( RandomPacking(38,689,128,310))
Raschig_ceramic_list.append( RandomPacking(51,651,95,210))
Raschig_ceramic_list.append( RandomPacking(76,561,69,120))

Raschig_13=Raschig_ceramic_list[0]
Raschig_25=Raschig_ceramic_list[1]
Raschig_38=Raschig_ceramic_list[2]
Raschig_51=Raschig_ceramic_list[3]
Raschig_76=Raschig_ceramic_list[4]

T_H2O=298.15            #T acqua ingresso al sistema
T_air=303.15            #T aria ingresso al sistema

air=Mixture(['air'],T=T_air,P=101325)
H2O=Chemical('H2O') #ho messo 30C perchè faccio stima conservativa, so che la T si alzerà mentre si assorbe HCl
HCl=Chemical('HCl')

# DATI INIZIALI #############################################################################################################################
V_in=10000              #Nm3/h
mHCl_in=5               #g/Nm3
ppmvout=150             #ppmv
yout=ppmvout/1e6        #frazione vapore di HCl in uscita dalla colonna (sezione di top)
Yout=yout/(1-yout)        #rapporto molare HCl in uscita dalla colonna

# UMIDITA' IN ARIA #########################################################################################################################
H2O_Psat=Chemical('H2O',T=T_air).Psat    #Pascal, frazione molare di acqua nell'aria

# DATI SOLUBILITA' PAG.2-132 PERRY 8 #########################################################################################################
P_parz20 = np.array([0.000044,0.00024,0.00178,0.0088,0.0428,0.205,1,4.9,23.5,105.5,399])/760    #mmHg
P_parz30 = np.array([0.000151,0.00077,0.00515,0.0234,0.106,0.48,2.17,9.9,44.5,188,627])/760     #mmHg
kgHClsuH2O = np.array([2.04,4.17,8.7,13.64,19.05,25,31.6,38.9,47,56.3,66.7])                    #g di HCl su 100g di soluzione

############## FRAZIONI MOLARI ############## Sono serviti per avere i punti in frazioni molari che poi abbiamo portato su MATLAB --->cftool
# xfrac=np.zeros(len(kgHClsuH2O))             #frazioni molari di HCl nella fase liquida
# y20=np.zeros(len(kgHClsuH2O))               #frazioni molari di H2O nella fase vapore a 20C
# y30=np.zeros(len(kgHClsuH2O))               #frazioni molari di H2O nella fase vapore a 30C
# for i in range (len(kgHClsuH2O)):
#     xfrac[i]=(kgHClsuH2O[i]/HCl.MW)/(kgHClsuH2O[i]/HCl.MW+100/H2O.MW)
#     y20[i]=P_parz20[i]
#     y30[i]=P_parz30[i]
    
# RAPPORTI MOLARI 
# X=np.zeros(len(kgHClsuH2O))                 #frazioni in peso di HCl nella fase liquida
# Y20=np.zeros(len(kgHClsuH2O))               #frazioni in peso di H2O nella fase vapore a 20C
# Y30=np.zeros(len(kgHClsuH2O))               #frazioni in peso di H2O nella fase vapore a 30C
# for i in range (len(kgHClsuH2O)):
#     X[i]=kgHClsuH2O[i]*(H2O.MW/HCl.MW)/100
#     Y20[i]=P_parz20[i]/(1-P_parz20[i])
#     Y30[i]=P_parz30[i]/(1-P_parz30[i])

# BILANCI DI MASSA ############################################################################################################################
R       = 8.314462618                                                       #costante gas
P       = 101325                                                            #P atmosferica in Pascal
nNm3    = P/(R*273.15)                                                      #moli/Nm3
yin     = (mHCl_in/HCl.MW)/(nNm3)                                           #frazione molare HCl all'ingresso il dato iniziale mHCl_in è in g/Nm3
Yin     = yin/(1-yin)                                                        #Rapporto molare HCl all'ingresso
Gv      = V_in*T_air/273.15                                                #m3/h aria totali della corrente di aria, acqua e HCl in ingresso l'acqua è a 30C
Gm      = P*Gv/(R*T_air)                                                    #moli/h totali della corrente di aria, acqua e HCl in ingresso l'acqua è a 30C
Gs      = Gm*(1-yin)                                        #moli/h aria secca 
Gp      = (Gm*(1-yin-y_water(303.15))*air.MW+Gm*yin*HCl.MW+Gm*y_water(303.15)*H2O.MW)/1000           #kg/h aria totali della corrente di aria, acqua e HCl in ingresso l'acqua è a 30C
Gp_air  = Gs*air.MW/1000                                                    #kg/h aria secca
nHCl_absorbed = Gs*(Yin-Yout)                                               #moli/h
mHCl_absorbed = nHCl_absorbed*HCl.MW/1000                                   #kg/h

#non sto considerando la condensazione dell'acqua 

[xout_max0]     = fsolve(intersection_in,0.122)                              
Xout_max       = xout_max0/(1-xout_max0)

#1) NO RICICLO (X2=0) : G*(yin-yout)=Gs*(Yin-Yout)=nHCl_absorbed=Ls*(Xout-X2)=L(xout-x2)
Xout0           = Xout_max
xout0           = frac(Xout0)
Ls_min          = nHCl_absorbed/(Xout0-0)                                   #moli/h di inerti
Lm_min_out      = nHCl_absorbed/(xout0-0)                                   #moli/h in uscita totali
Lm_min_in       = Lm_min_out-nHCl_absorbed                                  #moli/h in entrata totali
MW_out0         = (xout0)*HCl.MW+(1-xout0)*H2O.MW                           #Peso molecolare uscita
MW_in0          = H2O.MW                                                    #Peso molecolare entrata
Lp_min_out      = Lm_min_out*MW_out0/1000                                   #kg/h totali all'uscita
Lp_min_in       = Lm_min_in*MW_in0/1000                                     #kg/h totali entrata
Lp_min_H2O      = Ls_min*H2O.MW/1000                                        #kg/h di inerti (acqua) = kg/h totali all'entrata

m_min           = (Yin-Yout)/Xout0

############################################################## CALCOLO MOLARITA' #############################################################
massfrac_H2O    = Lp_min_in/(mHCl_absorbed+Lp_min_in)                       #frazione in peso di acqua all'uscita
massfrac_HCl    = mHCl_absorbed/(mHCl_absorbed+Lp_min_H2O)                  #frazione in peso di HCl all'uscita
Solution_Vol    = VolumeLiquidMixture(MW=[H2O.MW,HCl.MW],Tcs=[H2O.Tc,HCl.Tc],Pcs=[H2O.Pc,HCl.Pc],Vcs=[H2O.Vc,HCl.Vc],Zcs=[H2O.Zc,HCl.Zc],omegas=[H2O.omega,HCl.omega])
VL_Solution     = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-xout0,xout0],ws=[massfrac_H2O,massfrac_HCl],method='COSTALD_MIXTURE')
Sol_rho0        = MW_out0/(1000*VL_Solution)                                #kg/m3
M0              = nHCl_absorbed/(1000*Lp_min_out/Sol_rho0)                  #MOLARITA'

#VERIFICA FLG #i limiti del FLG si trovano a pag.494 del Kister, abbiamo fatto le medie tra le portate di ingresso e di uscita
flg_0         = ((Lp_min_out+Lp_min_in)/2)/((Gp+(Gp-mHCl_absorbed))/2)*((air.rho/H2O.rho)**0.5)
# if flg_0 <0.026:
#     print('FLG (without recycle) is too low')
# elif flg_0>2:
#     print('FLG (without recycle) is too high')

#VERIFICA CALORE ASSORBITO SENZA RICICLO
# https://www.ddpsinc.com/hcl-treatment si parla di circa 2100kj/kg assorbito
# https://www.dedietrich.com/en/solutions-and-products/halide-treatment/hcl-treament/absorption-hydrogen-chloride/adiabatic
Heat_absorption0= HeatAbsorption(Yin, Yout)                          #Potenza termica generata dall'assorbimento di HCl, in (kJ/h)
cp_mix0         = Mixture(['H2O', 'HCl'], zs=[1-xout0/2, xout0/2],T=303.15).Cpl/1000                                #kj/kg*K
delta_T0        = Heat_absorption0/((Lp_min_in+Lp_min_out)*cp_mix0)

#############################################################################################################################################################
#2) RICICLO MASSIMO (X2 = INTERSEZIONE SOLUBILITA' A 30°C CON Y2)
#Il massimo riciclo sarà individuato dalla retta che passa per l'intersezione tra la curva di equilibrio e l'ordinata yout
#Sapendo quindi due punti della retta è possibile calcolare l'equazione della retta
[xin_max]       = fsolve(intersection_out,0.088)
Xin_max         = ratio(xin_max)
m_max           = (Yin-Yout)/(Xout0-Xin_max)
Y0              = Yin-m_max*Xout0                                                   #termine noto retta operativa condizioni di massimo riciclo

# G*(yin-yout)=Gs*(Yin-Yout)=nHCl_absorbed=Ls*(Xout-X2)=L(xout-x2)
Ls_max          = nHCl_absorbed/(Xout0-Xin_max)                                     #mol/h di inerti (acqua)
Lm_max_in       = Ls_max*(1+Xin_max)                                                #mol/h in entrata in colonna totali
Lm_max_out      = Lm_max_in+nHCl_absorbed                                           #mol/h di in uscita  in colonna totali
MW_outREC       = (xout0)*HCl.MW+(1-xout0)*H2O.MW
MW_inREC        = (xin_max)*HCl.MW+(1-xin_max)*H2O.MW
Lp_max_in       = Lm_max_in*MW_inREC/1000                                           #kg/h totali
Lp_max_out      = Lm_max_out*MW_outREC/1000  
Lp_max_H2O      = Ls_max*H2O.MW/1000                                                #kg/h di inerti (acqua)

Ls_r            = Ls_max-Ls_min                                                     #mol/h inerti nel riciclo (che ha rapporto molare Xout0)
Lm_r            = Ls_r*(1+Xout0)
M0r=(Ls_max*Xout0)/(1000*Lp_max_out/Sol_rho0)                                       #ovviamente sarà usguale a M0 visto che escono entrambi allo stesso rapporto molare

#VERIFICA FLG
flg_0r  =((Lp_max_in+Lp_max_out)/2)/((Gp+(Gp-mHCl_absorbed))/2)*((air.rho/H2O.rho)**0.5)
# if flg_0r <0.026:
#     print('FLG (max recycle) is too low')
# elif flg_0r>2:
#     print('FLG (max recycle) is too high')

#VERIFICA CALORE ASSORBITO MAX RICICLO
cp_mix2         = Mixture(['H2O', 'HCl'], zs=[1-(xout0+xin_max)/2, (xout0+xin_max)/2],T=303.15).Cpl/1000              #kj/kg*K
delta_T0r       = Heat_absorption0/(((Lp_max_in+Lp_max_out)/2)*cp_mix2)
########################################################################################################################################

#################################### ACQUA CONDENSATA ########################################################################
#Supponiamo che la corrente gassosa condensi solo nella prima colonna e che il gas si raffreddi da 30C a 25C
# T_air_out=T_air-1.849
# T_airout1=T_air
# WaterCondensed1 = WaterCondensed(T_air,T_airout1,Gp_air)                                  #kg/h
# nWaterCondensed1= (WaterCondensed1/H2O.MW)*1000    
# T_airout2=T_airout1
# WaterCondensed2 = WaterCondensed(T_airout1,T_airout2,Gp_air)                                  #kg/h
# nWaterCondensed2= (WaterCondensed2/H2O.MW)*1000

####################################################################SISTEMA A DOPPIA COLONNA############################################################
#cominciamo dalla colonna dove esce il gas, dove supponiamo di lavorare a 25C con gas e solvente
#devo scegliere quanto assorbire nella colonna C2 le numeriamo in ordine seguendo il percorso del gas
###############################################COLONNA C2###############################################################################################
Yin_C1          = Yin
Yin_C2          = Yout*13.46
yin_c2          = Yin_C2/(1+Yin_C2)
nHCl_absorbed2  = Gs*(Yin_C2-Yout)
mHCl_absorbed2  = nHCl_absorbed2*HCl.MW/1000
Lpin_Global     = 158                                                           #kg/h
Lsin_Global     = Lpin_Global*1000/H2O.MW
Lsout_Global    = Lsin_Global                                                   #mol/h
Lpout_GlobalC2  = Lpin_Global+mHCl_absorbed2                                    #in uscita dalla colonna 2 e verso la colonna 1
m_Global        = Lsin_Global/Gs

rettaGlobal             = lambda x: +m_Global*(x)+Yout                                     #ridefinisco la retta operativa globale (sez.C1) con termine noto per far coincidere le due rette in Xout_2
intersezioneinC2        = lambda x: rettaGlobal(x)-Yin_C2
[Xout_C2]               = fsolve(intersezioneinC2,0.09)

Lph2o_C2        = 8000 #kg/h
Ls_C2           = Lph2o_C2*1000/H2O.MW
Ls_C2r          = Ls_C2-Lsin_Global #riciclo
m_C2            = Ls_C2/Gs
q_C2            = Yin_C2-m_C2*Xout_C2
retta_C2        = lambda x: m_C2*x+q_C2
intersezioneout_C2 = lambda x: retta_C2(x)-Yout
[Xin_C2]        = fsolve(intersezioneout_C2,0.05)


# COLONNA ASSORBIMENTO FISICO C2 (CONTROCORRENTE)
print('*************************** COLUMN 2*********************')
Yout_C2         = Yout                                         #Y out
yout_c2         = yout                                          #y out
xout_c2         = frac(Xout_C2)                                 #x out
xin_c2          = frac(Xin_C2)                                  #x in
Lm2             = Ls_C2+Ls_C2*Xin_C2                            #mol/h totali
Lp2in           = (Ls_C2*H2O.MW+Ls_C2*(Xin_C2)*HCl.MW)/1000
Lp2             = (Ls_C2*H2O.MW+0.5*Ls_C2*(Xin_C2+Xout_C2)*HCl.MW)/1000       #kg/h totali in
Lp2out          = Lp2in + mHCl_absorbed2
k2              = nHCl_absorbed2/nHCl_absorbed                  #frazione assorbita dalla seconda colonna

#medie
massfrac_H2O_2m    = (1-0.5*(xin_c2+xout_c2))*H2O.MW/((1-0.5*(xin_c2+xout_c2))*H2O.MW+0.5*(xin_c2+xout_c2)*HCl.MW)
massfrac_HCl_2m    = 0.5*(xin_c2+xout_c2)*HCl.MW/((1-0.5*(xin_c2+xout_c2))*H2O.MW+0.5*(xin_c2+xout_c2)*HCl.MW)
VL_Solution2m      = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-0.5*(xin_c2+xout_c2),0.5*(xin_c2+xout_c2)],ws=[massfrac_H2O_2m,massfrac_HCl_2m],method='COSTALD_MIXTURE')
Soluzione2m        = Mixture(['H2O','HCl'],zs=[1-0.5*(xout_c2+xin_c2),0.5*(xout_c2+xin_c2)],T=303.15)
rho_l2m            = (VL_Solution2m**-1)*Soluzione2m.MW/1000                            #densità liquida di miscela kg/m3

#in
massfrac_H2O_2in  = (1-(xin_c2))*H2O.MW/((1-(xin_c2))*H2O.MW+(xin_c2)*HCl.MW)
massfrac_HCl_2in  = (xin_c2)*HCl.MW/((1-xin_c2)*H2O.MW+(xin_c2)*HCl.MW)
VL_Solution2in    = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-(xin_c2),(xin_c2)],ws=[massfrac_H2O_2in,massfrac_HCl_2in],method='COSTALD_MIXTURE')
Soluzione2in      = Mixture(['H2O','HCl'],zs=[1-(xin_c2),(xin_c2)],T=303.15)
rho_l2in          = (VL_Solution2in**-1)*Soluzione2in.MW/1000                            #densità liquida di miscela kg/m3

#out
massfrac_H2O_2out  = (1-(xout_c2))*H2O.MW/((1-(xout_c2))*H2O.MW+(xout_c2)*HCl.MW)
massfrac_HCl_2out  = (xout_c2)*HCl.MW/((1-xout_c2)*H2O.MW+(xout_c2)*HCl.MW)
VL_Solution2out    = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-(xout_c2),(xout_c2)],ws=[massfrac_H2O_2out,massfrac_HCl_2out],method='COSTALD_MIXTURE')
Soluzione2out      = Mixture(['H2O','HCl'],zs=[1-(xout_c2),(xout_c2)],T=303.15)
rho_l2out          = (VL_Solution2out**-1)*Soluzione2out.MW/1000                            #densità liquida di miscela kg/m3

#CALCOLO MOLARITA #################################################################################################
Soluzione2out      = Mixture(['H2O','HCl'],zs=[1-(xout_c2),(xout_c2)],T=303.15)
MW_av2out          = Soluzione2out.MW
M2                 = nHCl_absorbed2/(1000*Lpout_GlobalC2/rho_l2out)                             #concentrazione M in uscita dalla corrente della colonna 2 e verso la colonna 1

##################################################################################################################
Soluzione2m      = Mixture(['H2O','HCl'],zs=[1-0.5*(xin_c2+xout_c2),0.5*(xin_c2+xout_c2)],T=303.15)
MW_av2          = Soluzione2m.MW

Lv2             = Lp2/rho_l2m
Lv2out          = Lp2out/rho_l2out                             #m3/h
Gp2             = Gp-(mHCl_absorbed-mHCl_absorbed2)
Gv2             = (Gs+Gs*(Yin_C2+Yout_C2)/2)*R*298.15/101325
#################################### ACQUA CONDENSATA ########################################################################
#ACQUA ENTRA SEMPRE A 25C, ALL'INIZIO IL GAS SI SCALDERA' ED EVAPORERA' ACQUA. L'ACQUA EVAPORANDO ASSORBE Q DAL LIQUIDO CHE SI RAFFREDDA.
# IL LIQUIDO A SUA VOLTA RAFFREDDA IL GAS, CHE USCIRA' ALLA STESSA T CON CUI E' ENTRATO, WaterCondensed2=0
T_air2          = T_H2O
# WaterCondensed2 = WaterCondensed(T_air2,T_air2,Gp_air)                               #kg/h

#################################### CHECK TEMPERATURA2 ########################################################################
# delta_T2        = np.zeros(1000)
delta_T2     = delta_T(Yin_C2,Yout_C2,Lp2) 


if delta_T2>5:
    print('Temperature rise is above 5°C !')
    print('Temperature rise is:','%.2f' %delta_T2,'°C')
else:
    print('Temperature rise is:','%.2f' %delta_T2,'°C')
    

#################################### CHECK FLG2 ########################################################################
GPDC2           = GPDC_21mmH2O_s
flg2            = Lp2/(Gp2)*(air.rho/rho_l2m)**0.5
if flg2<0.01:
    print('FLG is too low')
elif flg2>2:
    print('FLG is too high')
else:
    CP2=GPDC2(flg2)

############################ SCELGO RIEMPIMENTO2 #########################################################################
riemp2          = Mellapak_Plastic_250y
Fp2             = riemp2.Packing_Factor/3.28084                           # devo portarlo da 1/m a 1/feet

# pag 14-56 Perry deltaP flood, come scritto sopra c'è anche su Kister e se applichiamo intervallo di sicurezza del +15%
P_flood2        = (0.12*(Fp2)**0.7)*(25.4/0.3048)                                      # inch of H2O per foot of packing
CP_flood2       = GPDC_83mmH2O_s(flg2)                                                     #deltaP è circa 87mmH2O/m riempimento, noi abbiamo approssimato a 83 e abbiamo calcolato flood point velocity da li
usg2_flood      = CP_flood2/((Fp2)**(0.5)*(air.rho/(rho_l2m-air.rho))**0.5*(Soluzione2m.nul*1e6)**0.05)*0.3048        #m/s VELOCITA SUPERFICIALE GAS FLOOD

#Faccio Grafico GPDC con le curve per le perdite di carico per Packing Strutturati
########################################################################################################################
#   ATTENZIONE!!! lE UNITA DI MISURA DI CP(ORDINATE DEL GRAFICO) SONO IN UNITA' IMPERIALI!!!!!!
########################################################################################################################

# Ora vogliamo calcolare la velocità superficiale del gas, quindi usiamo la correlazione 14-140 per Perry 9
# La viscosita cinematica deve essere calcolata in centiStokes 1cS= 10^-6 m^2/s
# IL CP dal grafico di Kister&Gill per Structured Packing, formula è a pag 18 di Absorption3
usg2_           = CP2/((Fp2)**(0.5)*(air.rho/(rho_l2m-air.rho))**0.5*(Soluzione2m.nul*1e6)**0.05)*0.3048          #m/s VELOCITA SUPERFICIALE GAS
Ac2             = Gv2/(usg2_*3600)
Dc2             = (4*Ac2/np.pi)**0.5

[Dc2]           = np.around([Dc2],decimals=1) #metto -0.1 visto che ho già usato delle perdite di carico 15% inferiori, dovrei comunque essere bene al di sotto flooding
Ac2             = 0.25*np.pi*Dc2**2
usg2            = Gv2/(3600*Ac2)
F_Factor2       = usg2*air.rho**0.5
print('F-Factor =', '%.2f' %F_Factor2)

if usg2<0.801*usg2_flood and usg2>0.701*usg2_flood:                                         #questo margine ci consente di avere margine sufficiente per le incertezze associate al calcolo delle detaP
    print('usg2 is ', '%.2f' %(100*(usg2/usg2_flood)),' % of the flood point velocity')
else:
    print('ATTENTION!!!')
    print('usg2 is ', '%.2f' %(100*(usg2/usg2_flood)),' % of the flood point velocity')

usl2            = Lv2/(3600*Ac2)                                                                         #m/s VELOCITA' SUPERFICIALE LIQUIDO
# Ratio_Diam2     = Dc2/(riemp2.S/1000)
print('Column Area = ','%.3f' %Ac2,'m2', '& Diameter =','%.3f' %(Dc2*1000),'mm')
#CHECK RAPPORTI TRA DIAMETRI (prof consiglia di avere Dc>10*diametro riempimento)
# if Ratio_Diam2<10:
#     print('Diameter ratio is too low:','%.2f' %Ratio_Diam2)
# else:
#     print('Diameter ratio is:','%.2f' %Ratio_Diam2)
#CALCOLO WETTING RATE2
wr2             = Lp2/(Ac2*riemp2.area*rho_l2m) #volumetric liquid flowrate per unit wetted perimeter DALLE SLIDE ABSORPTION 3 PACKINGS PAG.9 di 22 si dice che 

#CALCOLO MIN WETTING RATE2
Lp_min_wr2      = 0.01*Gp2*((air.rho/rho_l2m)**-0.5)
mwr2            = Lp_min_wr2/(Ac2*riemp2.area*rho_l2m)

#CHECK MIN WETTING RATE2
if wr2<mwr2:
    print('wr is too low')

# Calcolo NTU int dy/(y-ye) tra y1(bottom colonna) e y2(top colonna) y1 e y2 sono frazioni molari

############################################ COLONNA 2 #############################################################

# ylm2=((yin_c2-Solubility30(xout_c2))-(yout_c2-Solubility30(xout_c2)))/np.log((yin_c2-Solubility30(xout_c2))/(yout_c2-Solubility30(xout_c2)))
# ntu2=(yin_c2-yout_c2)/(ylm2)
ntu2            = quad(integrale2,yout_c2,yin_c2)
NTU_oy2         = ntu2[0]

# Calcoliamo HTU2 con correlzione BRFT
L2              = Lp2/(rho_l2m*Ac2)
# si vede a pag.20 di Absorption 3 Packings
if L2<40:
    hL2=0.0169*riemp2.area**0.83*L2**0.37*(Soluzione2m.mul/Chemical('H2O',T=293.15).mu)**0.25
else:
    hL2=0.0075*riemp2.area**0.83*L2**0.59*(Soluzione2m.mul/Chemical('H2O',T=293.15).mu)**0.25

print('Liquid Hold up is :','%.2f'%hL2 ,'%') 

            
Uge2=usg2/(riemp2.void*(1-hL2/100)*np.sin(np.deg2rad(riemp2.alpha)))   #[m/s]
Ule2=usl2/(riemp2.void*(hL2/100)*np.sin(np.deg2rad(riemp2.alpha)))     #[m/s]

#Vedi pag.17 di Sizing Packed Columns

#GAS PHASE TRANSFER2 (Fuller)
Diff_2gas       = (1.013*(10**-7)*(T_air**1.75)*(1/air.MW + 1/HCl.MW )**0.5)/(P/1e5*((19.5+1.98)**(1/3)+(20.1)**(1/3)))**2         #m2/s
Sc_2gas         = air.mu/(air.rho*Diff_2gas)                                                            #Schmidt
Re_2gas         = (Uge2+Ule2)*air.rho*riemp2.S/(1000*air.mu)                                            #Reynolds
kg_brft2        = 0.054*(Diff_2gas/(riemp2.S/1000))*(Re_2gas**0.8)*(Sc_2gas**0.33)                      #m/s vedi pag.5-82 Perry8
Sh_2gas         = kg_brft2*(riemp2.S/1000)/Diff_2gas                                                    #Sherwood (mass transfer convettivo/diffusivo)

#LIQUID PHASE TRANSFER (Wilke-Chang)
Diff_2liq       = 1.173*(10**(-13))*(2.6*H2O.MW)**0.5*298/(1000*H2O.mu*(1000/H2O.rhom)**0.6)
Sc_2liq         = Soluzione2m.mul/(rho_l2m*Diff_2liq)
kL2             = 2*(Diff_2liq*0.9*Ule2/(np.pi*riemp2.S/1000))**0.5                                     #m/s

Q2              = Lv2/3600                                                                                  #m3/s volumetric flowrate 
Lp_2            = 4*(riemp2.S/1000)/((riemp2.B/1000)*(riemp2.h/1000))*Ac2                                   
sigma2          = Soluzione2m.sigma                                                                          
ae2             = riemp2.area*1.34*(((rho_l2m/sigma2)*9.81**(1/3)*(Q2/Lp_2)**(4/3)))**0.116

HTU_g2          = usg2/(kg_brft2*ae2)
HTU_l2          = usl2/(kL2*ae2)

slope_op2       = Ls_C2/Gs
vec2            = np.linspace(xout_c2,xin_c2,501)
slope_eq2       = linregress(vec2,Solubility30r(vec2))
slope_eq2       = slope_eq2[0]
LAMBDA2         = slope_eq2/slope_op2                                                                 #slope equilibrium/slope operative line

HTU_oy2         = HTU_g2+HTU_l2*LAMBDA2
                                                            #https://pubs.acs.org/doi/10.1021/ie940406i pag. 5-82 Perry8
print('NTU2 =','%.3f' %NTU_oy2)
print('HTU2 =','%.0f' %(HTU_oy2*1000),'mm')
Z2=HTU_oy2*NTU_oy2
print('Z2 =','%.0f' %(Z2*1000), 'mm')

################################## PERDITE DI CARICO DISTRIBUITE ############################################################################
if GPDC2==GPDC_8mmH2O_s:
    deltaP2=8
elif GPDC2==GPDC_21mmH2O_s:
    deltaP2=21
elif GPDC2==GPDC_42mmH2O_s:
    deltaP2=42
elif GPDC2==GPDC_21mmH2O_s:
    deltaP2=83
elif GPDC2==GPDC_21mmH2O_s:
    deltaP2=125

if F_Factor2>2.7:
    CorrectionP2=(usg2/usg2_)**2.2 #past loading point
else:
    CorrectionP2=(usg2/usg2_)**1.75

deltaP_colonna2=deltaP2*Z2*CorrectionP2 #mmH2O

print()
print('DeltaP Packing 2 =','%.2f' %(deltaP_colonna2/100) ,'kPa')
print('Delta P Flood2 =', '%.2f' %(Z2*P_flood2/100), 'kPa' )

# Dimensionamento bocchelli sia in ingresso che in uscita in quanto la portata cambia di poco
D_bocchelli2=0.5 #metri
A_bocchelli2=0.25*np.pi*D_bocchelli2**2
V_Bocchelli2=(Gv/3600)/(A_bocchelli2) #m/s
if V_Bocchelli2>20 or V_Bocchelli2<10:
    print('Modificare diametro bocchelli gas')

Nuscite2=100
Asingolauscita2=A_bocchelli2/Nuscite2
# Dprova2=0.08
# Aentrata2=0.25*np.pi*Dprova2**2
# Aentrata2tot=Aentrata2*Nuscite2

delta_P2_bocchelli=air.rho*V_Bocchelli2**2

# dimensionamento distributore liquido
l2=Lv2/3600
DripDensity2=160
Nfori2=DripDensity2*Ac2
diam_fori2=0.007 
A_fori2=Nfori2*0.25*np.pi*diam_fori2**2
Vfori2=l2/A_fori2
battente2=(l2/A_fori2)**2/(2*9.81) #metri battente nella canalina del distributore

#bocchelli entrata liquido
diam_bocchello_L2in=0.055 #2pollici
A_bocchello_L2in=0.25*np.pi*diam_bocchello_L2in**2
Vbocchello_L2in=(Lv2/3600)/A_bocchello_L2in

if Vbocchello_L2in<0.8 or Vbocchello_L2in>1.2:
    print('Modificare diam bocchelli liquido in')

#bocchelli uscita liquido
diam_bocchello_L2out=0.085
A_bocchello_L2out=0.25*np.pi*diam_bocchello_L2out**2
Vbocchello_L2out=(Lv2out/3600)/A_bocchello_L2out

if Vbocchello_L2out<0.3 or Vbocchello_L2out>0.5:
    print('Modificare diam bocchelli liquido out')

# Demister: tipo di demister T-01-P *2 https://vff.com/fileadmin/user_upload/Downloads/210323_RZ_Demister_E_WEB.pdf
#standard Sulzer KnitMesh demister polypropylene 9048 https://www.sulzer.com/-/media/files/products/separation-technology/feed-inlet-devices/gas_liquid_separation_technology.ashx
deltaP_demister2=dP_demister_dry_Setekleiv_Svendsen(S=360, voidage=.97, vs=usg2, rho=air.rho, mu=air.mu, L=0.2) #Pa

efficiency_demister2=separation_demister_ElDessouky(usg2, 0.97, 0.0002, 0.005)  

deltaP2tot=deltaP_colonna2*10+delta_P2_bocchelli+deltaP_demister2*2

print()
print('deltaP tot Column 2 =', '%.3f' %(deltaP2tot/1000), 'kPa')

#dimensionamento fondo colonna
hfondo2=0.6 #metri di battende di liquido fondo colonna
kgsolfondocolonna2= Ac2*rho_l2out*hfondo2
T_residenza2=(kgsolfondocolonna2)/(Lp2out/60)

if T_residenza2<1 or T_residenza2>8:
    print('Cambiare tempo residenza sul fondo colonna')

###############################################COLONNA 1###############################################################################################
print('*************************** COLUMN 1*********************')


Yout_C1                 = Yin_C2
yout_c1                 = Yout_C1/(1+Yout_C1)
nHCl_absorbed1          = Gs*(Yin_C1-Yout_C1)
mHCl_absorbed1          = nHCl_absorbed1*HCl.MW/1000
Lpout_Global            = Lpin_Global+mHCl_absorbed2 + mHCl_absorbed1
Lpin1_Global            = (Lsin_Global*H2O.MW+Lsin_Global*Xout_C2*HCl.MW)/1000

LpH2O_C1                = 8000 #kg/h
Ls_C1                   = LpH2O_C1*1000/H2O.MW
Ls_C1r                  = Ls_C1 - Lsin_Global
# Xin_C1calc              = (Ls_C1r*Xout_C1+Lsin_Global*Xout_C2)/(Ls_C1in)
Ls_C1medio              = Ls_C1
# Lp_C1medio              = Ls_C1*H2O.MW+Ls_C1*(Xin_C1+)


Yin_C1                  = Yin
intersezioneinC1        = lambda x: rettaGlobal(x)-Yin_C1
[Xout_C1]               = root(intersezioneinC1,0.15).x

m_C1                    = Ls_C1/Gs
q_C1                    = Yin_C1-m_C1*Xout_C1
retta_C1                = lambda x: m_C1*x+q_C1
intersezioneout_C1      = lambda x: retta_C1(x)-Yout_C1
[Xin_C1]                = fsolve(intersezioneout_C1,0.1)

Lp1                     = (Ls_C1*H2O.MW+0.5*(Xin_C1+Xout_C1)*Ls_C1*HCl.MW)/1000
Lp1in                   = (Ls_C1*H2O.MW+(Xin_C1)*Ls_C1*HCl.MW)/1000
Lp1out                  = Lp1in+mHCl_absorbed1
deltaT1                 = delta_T(Yin_C1,Yin_C2,Lp1)
if deltaT1>5:
    print('Temperature rise is above 5°C !')
    print('Temperature rise is:','%.2f' %deltaT1,'°C')
else:
    print('Temperature rise is:','%.2f' %deltaT1,'°C')

xin_c1                  = frac(Xin_C1)
xout_c1                 = frac(Xout_C1)

#medie
massfrac_H2O_1m    = (1-0.5*(xin_c1+xout_c1))*H2O.MW/((1-0.5*(xin_c1+xout_c1))*H2O.MW+0.5*(xin_c1+xout_c1)*HCl.MW)
massfrac_HCl_1m    = 0.5*(xin_c1+xout_c1)*HCl.MW/((1-0.5*(xin_c1+xout_c1))*H2O.MW+0.5*(xin_c1+xout_c1)*HCl.MW)
VL_Solution1m      = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-0.5*(xin_c1+xout_c1),0.5*(xin_c1+xout_c1)],ws=[massfrac_H2O_1m,massfrac_HCl_1m],method='COSTALD_MIXTURE')
Soluzione1m        = Mixture(['H2O','HCl'],zs=[1-0.5*(xout_c1+xin_c1),0.5*(xout_c1+xin_c1)],T=303.15)
rho_l1m            = (VL_Solution1m**-1)*Soluzione1m.MW/1000                            #densità liquida di miscela kg/m3

#in
massfrac_H2O_1in  = (1-(xin_c1))*H2O.MW/((1-(xin_c1))*H2O.MW+(xin_c1)*HCl.MW)
massfrac_HCl_1in  = (xin_c1)*HCl.MW/((1-xin_c1)*H2O.MW+(xin_c1)*HCl.MW)
VL_Solution1in    = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-(xin_c1),(xin_c1)],ws=[massfrac_H2O_1in,massfrac_HCl_1in],method='COSTALD_MIXTURE')
Soluzione1in      = Mixture(['H2O','HCl'],zs=[1-(xin_c1),(xin_c1)],T=303.15)
rho_l1in          = (VL_Solution1in**-1)*Soluzione1in.MW/1000                            #densità liquida di miscela kg/m3

#out
massfrac_H2O_1out  = (1-(xout_c1))*H2O.MW/((1-(xout_c1))*H2O.MW+(xout_c1)*HCl.MW)
massfrac_HCl_1out  = (xout_c1)*HCl.MW/((1-xout_c1)*H2O.MW+(xout_c1)*HCl.MW)
VL_Solution1out    = Solution_Vol.calculate(T=T_H2O,P=101325,zs=[1-(xout_c1),(xout_c1)],ws=[massfrac_H2O_1out,massfrac_HCl_1out],method='COSTALD_MIXTURE')
Soluzione1out      = Mixture(['H2O','HCl'],zs=[1-(xout_c1),(xout_c1)],T=303.15)
rho_l1out          = (VL_Solution1out**-1)*Soluzione1out.MW/1000                            #densità liquida di miscela kg/m3

#CALCOLO MOLARITA #################################################################################################
M1                 = nHCl_absorbed/(1000*Lpout_Global/rho_l1out) 

MW_av1          = Soluzione1m.MW
# rho_l1          = (VL_Solution1m**-1)*MW_av1/1000                    #densità liquida di miscela kg/m3

#################################### CHECK FLG1 ######################################################################## pag.493 Kister
flg1            = (Lp1)/(Gp)*(air.rho/rho_l1m)**0.5
GPDC1           = GPDC_21mmH2O_s
if flg1<0.001:                                                          #i limiti di flg sono stati presi a pag.648 Kister
    print('FLG is too low')
elif flg1>2:
    print('FLG is too high')
else:
    CP1=GPDC1(flg1)

############################ SCELGO RIEMPIMENTO1 #########################################################################
riemp           = Mellapak_Plastic_250y
Fp1             = riemp.Packing_Factor/3.28084                          #devo portarlo da 1/m a 1/feet

#Stando a cio che è scritto a pag 14-56 del Perry possiamo calcolare il deltaP flood
#prima pero dobbiamo scegliere il tipo di rimepimento
#KISTER PAG 482 P_flood +-15% quindi per essere conservativi possiamo prendere P_flood*0.85

P_flood1        = (0.12*(Fp1)**0.7)*(25.4/0.3048)                       #inch of H2O per foot of packing  mmH2O/m 1mmH2O=9.80665 Pa, 1m=3.28084feet
CP_flood        = GPDC_83mmH2O_s(flg1)                                  #deltaP è circa 87mmH2O/m riempimento, noi abbiamo approssimato a 83 e abbiamo calcolato flood point velocity da li
usg1_flood      = CP_flood/((Fp1)**(0.5)*(air.rho/(rho_l1m-air.rho))**0.5*(Soluzione1m.nul*1e6)**0.05)*0.3048        #m/s VELOCITA SUPERFICIALE GAS
########################################################################################################################
#   ATTENZIONE!!! lE UNITA DI MISURA DI CP(ORDINATE DEL GRAFICO) SONO IN UNITA' IMPERIALI!!!!!!
########################################################################################################################

# Ora vogliamo calcolare la velocità superficiale del gas, quindi usiamo la correlazione 14-140 per Perry 9
# La viscosita cinematica deve essere calcolata in centiStokes 1cS= 10^-6 m^2/s
# IL CP dal grafico di Kister&Gill per Structured Packing, formula è a pag 18 di Absorption3
usg1_           = CP1/((Fp1)**(0.5)*(air.rho/(rho_l1m-air.rho))**0.5*(Soluzione1m.nul*1e6)**0.05)*0.3048        #m/s VELOCITA SUPERFICIALE GAS
Ac1             = Gv/(usg1_*3600)
Dc1             = (4*Ac1/np.pi)**0.5
[Dc1]           = np.around([Dc1-0.1],decimals=1) #metto +0.3 visto per abbassare F-factor
Ac1             = 0.25*np.pi*Dc1**2
usg1            = Gv/(3600*Ac1)
F_Factor1       = usg1*air.rhog**0.5
print('F-Factor =', '%.2f' %F_Factor1)

if usg1<0.801*usg1_flood and usg1>0.701*usg1_flood:
    print('usg1 is ', '%.2f' %(100*(usg1/usg1_flood)),' % of the flood point velocity') #questo margine ci consente di avere margine sufficiente per le incertezze associate al calcolo delle detaP
else:
    print('ATTENTION!!!')
    print('usg1 is ', '%.2f' %(100*(usg1/usg1_flood)),' % of the flood point velocity')

Lv1             = Lp1/rho_l1m
Lv1out          = Lp1out/rho_l1out  
usl1            = Lv1/(3600*Ac1)                                                                      #m/s VELOCITA' SUPERFICIALE LIQUIDO
# Ratio_Diam1     = Dc1/(riemp.S/1000)
print('Column Area = ','%.3f' %Ac1,'m2', '& Diameter =','%.3f' %(Dc1*1000),'mm')

#CHECK RAPPORTI TRA DIAMETRI (prof consiglia di avere Dc>10*diametro riempimento)
# if Ratio_Diam1<10:
#     print('Diameter ratio is too low:','%.2f' %Ratio_Diam1)
# else:
#     print('Diameter ratio is:','%.2f' %Ratio_Diam1)
#CALCOLO WETTING RATE
wr1         = Lp1/(Ac1*riemp.area*rho_l1m) #volumetric liquid flowrate per unit wetted perimeter DALLE SLIDE ABSORPTION 3 PACKINGS PAG.9 di 22 si dice che 
# wr1top      = Lp1in/(Ac1*riemp.area*rho_l1)
# wr1bottom   = Lp1out/(Ac1*riemp.area*rho_l1)
#CALCOLO MIN WETTING RATE
flg_min_wr  = 0.01
Lp_min_wr   = flg_min_wr*Gp*((air.rho/rho_l1m)**-0.5)
mwr1        = Lp_min_wr/(Ac1*riemp.area*rho_l1m)

#CHECK MIN WETTING RATE
if wr1<mwr1:              #facciamo il check sul wr minomo, cioè in top alla colonna
    print('wr is too low')

# Calcolo NTU int dy/(y-ye) tra y1(bottom colonna) e y2(top colonna) y1 e y2 sono frazioni molari
# Per ye si può usare direttamente Solubility30
############################################ COLONNA 1 #############################################################
ylm=((yin-Solubility30(xout_c1))-(yin_c2-Solubility30(xin_c1)))/np.log((yin-Solubility30(xout_c1))/(yin_c2-Solubility30(xin_c1)))
ntu_approssimata=(yin-yin_c2)/(ylm)
# NTU_oy1=ntu_approssimata
m1=Ls_C1/Gs
q1=Yin_C1-m1*Xout_C1
ntu         = quad(integrale,yout_c1,yin)
NTU_oy1     = ntu[0]

# Calcoliamo HTU con correlzione BRFT
L1          = Lp1/(rho_l1m*Ac1)
# si vede a pag.20 di Absorption 3 Packings
if L1<40:
    hL1=0.0169*riemp.area**0.83*L1**0.37*(Soluzione1m.mul/Chemical('H2O',T=293.15).mu)**0.25
else:
    hL1=0.0075*riemp.area**0.83*L1**0.59*(Soluzione1m.mul/Chemical('H2O',T=293.15).mu)**0.25
            
print('Liquid Hold up is :','%.2f'%hL1 ,'%') 
    
Uge1        = usg1/(riemp.void*(1-hL1/100)*np.sin(np.deg2rad(riemp.alpha)))   #[m/s]
Ule1        = usl1/(riemp.void*(hL1/100)*np.sin(np.deg2rad(riemp.alpha)))     #[m/s]

#Vedi pag.17 di Sizing Packed Columns
# Equazione di Fuller Coulson 6 pag 331

############################### GAS PHASE TRANSFER ########################################################################################
Diff_1gas   = (1.013*(10**-7)*(T_air**1.75)*(1/air.MW + 1/HCl.MW )**0.5)/(P/1e5*((19.5+1.98)**(1/3)+(20.1)**(1/3)))**2         #m2/s
Sc_1gas     = air.mu/(air.rho*Diff_1gas)                                                      #Schmidt
Re_1gas     = (Uge1+Ule1)*air.rho*riemp.S/(1000*air.mu)                                       #Reynolds
kg_brft1    = 0.054*(Diff_1gas/(riemp.S/1000))*(Re_1gas**0.8)*(Sc_1gas**0.33)                 #m/s vedi pag.5-82 Perry8
Sh_1gas     = kg_brft1*(riemp.S/1000)/Diff_1gas                                               #Sherwood (mass transfer convettivo/diffusivo)

#################################### LIQUID PHASE TRANSFER equaione di wilke chang 333 ###################################################
Diff_1liq   = 1.173*(10**(-13))*(2.6*H2O.MW)**0.5*298/(1000*H2O.mu*(1000/H2O.rhom)**0.6)
Sc_1liq     = Soluzione1m.mul/(rho_l1m*Diff_1liq)
kL1         = 2*(Diff_1liq*0.9*Ule1/(np.pi*riemp.S/1000))**0.5

################################### EFFECTIVE AREA #######################################################################################
Q1          = Lv1/3600                                                                              #m3/s volumetric flowrate 
Lp_         = 4*(riemp.S/1000)/((riemp.B/1000)*(riemp.h/1000))*Ac1                                 #N/m Attenzione a riportare tutto in metri
sigma1      = Soluzione1m.sigma
ae1         = riemp.area*1.34*(((rho_l1m/sigma1)*9.81**(1/3)*(Q1/Lp_)**(4/3)))**0.116

#CALCOLO HTU GAS E LIQUIDA
HTU_g1      = usg1/(kg_brft1*ae1)
HTU_l1      = usl1/(kL1*ae1)

#CALCOLO ALTEZZA C2
slope_op1   = Ls_C1/Gs
slope_eq1   = linregress(np.linspace(xout_c1,xin_c1,501),Solubility30r(np.linspace(xout_c1,xin_c1,501)))
slope_eq1   = slope_eq1[0]
LAMBDA1     = slope_eq1/slope_op1                                                              #slope equilibrium/slope operative line
HTU_oy1     = HTU_g1+HTU_l1*LAMBDA1
Z1          = HTU_oy1*NTU_oy1                                                                       #https://pubs.acs.org/doi/10.1021/ie940406i pag. 5-82 Perry8
print('NTU1 =','%.3f' %NTU_oy1)
print('HTU1 =','%.0f' %(HTU_oy1*1000),'mm')
print('Z1 =','%.0f' %(Z1*1000), 'mm')

################################## PERDITE DI CARICO DISTRIBUITE ############################################################################
if GPDC1==GPDC_8mmH2O_s:
    deltaP1=8
elif GPDC1==GPDC_21mmH2O_s:
    deltaP1=21
elif GPDC1==GPDC_42mmH2O_s:
    deltaP1=42
elif GPDC1==GPDC_21mmH2O_s:
    deltaP1=83
elif GPDC1==GPDC_21mmH2O_s:
    deltaP1=125

if F_Factor1>3:
    CorrectionP1=(usg1/usg1_)**2.2 #past loading point
else:
    CorrectionP1=(usg1/usg1_)**1.75

deltaP_colonna1     = deltaP1*Z1*CorrectionP1
deltaP_flood1       = Z1*P_flood1

print()
print('Delta P column 1 =', '%.2f' %(deltaP_colonna1/100),'kPa')
print('Delta P Flood1 is =', '%.2f' %(deltaP_flood1/100), 'kPa' )


#Dimensionamento Bocchelli sia in ingresso che in uscita in quanto la portata cambia di poco

D_bocchelli1=0.5 #metri
V_Bocchelli1=(Gv/3600)/(0.25*np.pi*D_bocchelli1**2) #m/s
if V_Bocchelli1>20 or V_Bocchelli1<10:
    print('Modificare diametro bocchelli gas')

delta_P1_bocchelli=air.rho*V_Bocchelli2**2

# dimensionamento distributore liquido
l1=Lv1/3600
DripDensity1=160
Nfori1=DripDensity1*Ac1
diam_fori1=0.007 
A_fori1=Nfori1*0.25*np.pi*diam_fori1**2
Vfori1=l1/A_fori1
battente1=(l1/A_fori1)**2/(2*9.81) #metri battente nella canalina del distributore

#bocchelli entrata liquido
diam_bocchello_L1in=0.055 #2pollici
A_bocchello_L1in=0.25*np.pi*diam_bocchello_L1in**2
Vbocchello_L1in=(Lv1/3600)/A_bocchello_L1in

if Vbocchello_L1in<0.8 or Vbocchello_L1in>1.2:
    print('Modificare diam bocchelli liquido in')

#bocchelli uscita liquido
diam_bocchello_L1out=0.085
A_bocchello_L1out=0.25*np.pi*diam_bocchello_L1out**2
Vbocchello_L1out=(Lv1out/3600)/A_bocchello_L1out

if Vbocchello_L1out<0.3 or Vbocchello_L1out>0.5:
    print('Modificare diam bocchelli liquido out')
    
#Demister: tipo di demister T-01-P *2 https://vff.com/fileadmin/user_upload/Downloads/210323_RZ_Demister_E_WEB.pdf
#standard Sulzer KnitMesh demister polypropylene 9048 https://www.sulzer.com/-/media/files/products/separation-technology/feed-inlet-devices/gas_liquid_separation_technology.ashx
deltaP_demister1=dP_demister_dry_Setekleiv_Svendsen(S=360, voidage=.97, vs=usg1, rho=air.rho, mu=air.mu, L=0.2)

efficiency_demister1=separation_demister_ElDessouky(usg1, 0.97, 0.0002, 0.005)    

deltaP1tot=deltaP_colonna1*10+delta_P1_bocchelli+deltaP_demister1*2

#dimensionamento fondo colonna
hfondo1=0.6 #metri di battende di liquido fondo colonna
kgsolfondocolonna1= Ac1*rho_l1out*hfondo1
T_residenza1=(kgsolfondocolonna1)/(Lp1out/60)
print()
print('deltaP tot Column 1 =', '%.3f' %(deltaP1tot/1000), 'kPa')


if T_residenza1<1 or T_residenza1>8:
    print('Cambiare tempo residenza sul fondo colonna')

#ASSORBIMENTO CHIMICO
print('****************************CHEMICAL ABSORPTION COLUMN *********************')

kgNaOHsukgSol   = 0.00004                              #https://www.sigmaaldrich.com/IT/it/product/sigald/415413
NaOH            = Chemical('NaOH')
Yin_3           = Yout_C2                                       #Y in HCl
yin_3           = frac(Yin_3)                                   #y in HCl
intersezionechim2 = lambda x: Solubility30(x)-yin_3
[xout_3]        = fsolve(intersezionechim2,0.08)
xout_3          = 0                                    
Xout_3          = ratio(xout_3)
yout_3          = 1e-5                                    #y out HCl
Yout_3          = yout_3/(1-yout_3)                             #Y out HCl
intersezionechim1 = lambda x: Solubility30(x)-yout_3
[xin_3]         = fsolve(intersezionechim1,0.04)
Xin3_           = ratio(xin_3)
nHCl_reacted    = Gs*(Yin_3-Yout_3)                         #mol/h
mHCl_reacted    = nHCl_reacted*HCl.MW/1000                  #kg/h     
c1              = 1                                         #coefficiente stechiometrico (quante moli di NaOH reagiscono ogni mole di HCl consumata)
nNaOH_reacted   = c1*nHCl_reacted                   #mol/h
mNaOH_reacted   = nNaOH_reacted*NaOH.MW/1000        #kg/h
mSolution       = mNaOH_reacted/kgNaOHsukgSol           #kg/h
Solution        = Mixture(['H2O','NaCl'], ws=[(1-0.1),0.1],T=298.15)

NaCl            = Chemical('NaCl')
kgsale          = nNaOH_reacted*(NaCl.MW)/1000 #kg NaCl che si formano in un ora

M_NaOH          = (kgNaOHsukgSol*1000/NaOH.MW)/(1000/Solution.rhol)
print('NaOH is','%.3f' %M_NaOH,'M')
pOH             = -np.log10(M_NaOH)
print('pOH is','%.0f' %pOH)


Gp3             = Gp-mHCl_absorbed1-mHCl_absorbed2        #kg/h
Gm3             = Gm-nHCl_absorbed1-nHCl_absorbed2        #mol/h
Gv3             = Gp3/air.rho                             #m3/h
mwr3            = 0.1                                    #m3/h/m
riemp3          = Raschig_38

GPDC3=GPDC_42mmH2O_r

K4              = 1.2
intersezione_chemical= lambda x: GPDC3(x)-K4
[flg3]          = fsolve(intersezione_chemical,0.1)

K4_flooding     = GPDC_flooding(flg3)

Lp3             = flg3*Gp3*(air.rho/Solution.rhol)**(-0.5)              #kg/h totali
Lp3_H2O         = Lp3*(1-kgNaOHsukgSol)                                 #kg/h acqua
Ls3             = 1000*Lp3_H2O/(H2O.MW)                                 #mol/h
Lv3             = Lp3/Solution.rho                                      #m3/h

usg3_            = ((K4*(Solution.rhol-air.rho)/(13.1*riemp3.Packing_Factor*air.rho))*(Solution.mul/Solution.rhol)**(-0.1))**0.5               #m/s
usg3_flooding   = ((K4_flooding*(Solution.rhol-air.rho)/(13.1*riemp3.Packing_Factor*air.rho))*(Solution.mul/Solution.rhol)**(-0.1))**0.5  

Ac3             = (Gv3/3600)/usg3_                         #m2
Dc3             = (4*Ac3/np.pi)**0.5
[Dc3]           = np.around([Dc3-0.1],decimals=1)
Ac3             = 0.25*np.pi*Dc3**2
usg3            = Gv3/(Ac3*3600)
F_Factor3       = usg3*air.rho**0.5
print('F-Factor =', F_Factor3)

print('Ac3 =','%.2f' %(Ac3),'m', '& Dc3 =','%.2f' %(Dc3),'m')

DiamRatio3      = Dc3/(riemp3.size/1000)
#CHECK RAPPORTO DIAMETRI
if DiamRatio3<10:
    print('Diameter ratio is too small')
L3=Lv3/Ac3
usl3            = Lv3/(3600*Ac3)
epsilon         = 1-(riemp3.bulk_density/2500) #grado di vuoto
hL3             = 1.03/epsilon*(Solution.mul/Solution.rhol)**(1/6)*(usl3*riemp3.area)**0.5
wr3             = Lv3/(riemp3.area*Ac3) #m3/h/m

if wr3<mwr3:
    print("ATTENTION!!! wr3<mwr3")
########################## ONDA CORRELATIONS ################################################################################################
########################## LIQUID COEFFICIENTS ################################################################################################
# ATENZIONE UNITA DI MISURA!!
Lw3             = Lp3/(3600*Ac3)                    #kg/m2*s
sigmariemp      = 61e-3                             #N/m CERMICS DA PAG.528 KISTER

# Wilke-Chang
DiffL3          = 1.173*(10**(-13))*(1*Solution.MW)**0.5*298/(1000*Solution.mu*(1000/Solution.rhom))
Re_L3           = Solution.rhol*usl3/(riemp3.area*Solution.mul)
Fr_L3           = usl3**2*riemp3.area/9.81
We_L3           = Lw3**2/(Solution.rhol*Solution.sigma*riemp3.area)
aw3             = riemp3.area*(1-np.exp(-1.45*(sigmariemp/Solution.sigma)**(0.75)*(Re_L3)**0.1*(Fr_L3)**(-0.05)*(We_L3)**0.2))
Sc_L3           = Solution.nu/DiffL3
kL3             = 0.0051*(Lw3/(aw3*Solution.mu))**(2/3)*(Sc_L3)**(-0.5)*(riemp3.area*riemp3.size/1000)**0.4*(Solution.rhol/(Solution.mul*9.81))**(-1/3)
Hl3             = (Lw3)/(kL3*aw3*Solution.rhol)

############################### GAS COEFFICIENTS ###########################################################################################
if riemp3.size>15:
    K5=5.23
else:
    K5=2
Vw3             = Gp3/(3600*Ac3)                  #kg/m2*s
Vm3             = Vw3/air.MW                      #kmol/m2*s
DiffG3          = Diff_1gas
Sc_G3           = air.nu/DiffG3
kG3             = K5*(Vw3/(riemp3.area*air.mu))**0.7*(Sc_G3)*(1/3)*(riemp3.area*riemp3.size/1000)**(-2)*riemp3.area*DiffG3/(8.31*303)/1000 #kmol/[s*m2*(N/m2)]
Hg3             = Vm3/(kG3*aw3*1.01325e5)

#KINETICS
# https://kinetics.nist.gov/kinetics/ReactionSearch;jsessionid=425ED7347EA4C9C9F075073E46381923?r0=1310732&r1=7647010&r2=0&r3=0&r4=0&p0=7647145&p1=7732185&p2=0&p3=0&p4=0&expandResults=true&
Pa_in           = P*yin_3       #Pa
Pa_out          = P*yout_3      #Pa
cB              = M_NaOH        #mol/L

#Se usasassi una soluzione 1M per riportare la concentrazione a cB mi servirebbero NaOH_reint litri di soluzione e per mantenere la portata di liquido costante in colonna dovrei spurgarne altrettanti

#Spurgo NaCl e reintegro NaOH
Spurgo=kgsale/0.1                                       #kg/h soluzione spurgata supponiamo di operare al 10% peso di NaCl in soluzione
LSpurgo=1000*Spurgo/Solution.rhol                       #l/h soluzione spurgata
NaOH_reint=nNaOH_reacted+cB*LSpurgo                     #moli/h consumate
mNaOH_reint=NaOH.MW*NaOH_reint/1000                     #kg/h 
print('Purge is','%.2f' %LSpurgo,'l/h')                 #spurgo
print('NaOH Reint is', '%.2f' %mNaOH_reint,'kg/h')         
# ReintNaOH=mNaOH_reacted+Mspurgo*
if kG3*Pa_in<kL3*cB/c1:                                 #LEVENSPIEL 3RD ED. PAG. 532
    print('The reaction is instantaneous at the top')
else:
    print('The reaction is NOT instantaneous at the top')
    
#cB possiamo considerarlo costante visto che si consumano 2kg/h su 48mila kg/h di soluzione in entrata
if kG3*Pa_out<kL3*cB/c1:
    print('The reaction is instantaneous at the bottom')
else:
    print('The reaction is NOT instantaneous at the bottom')
  
slope_op3       = Ls3/Gs
x_3             = np.linspace(0.057,0.096,100)
slope_eq3       = linregress(x_3,Solubility30r(frac(x_3)))
slope_eq3       = slope_eq3[0]
lambda3         = slope_eq3/slope_op3

if kG3*Pa_in<kL3*cB/c1 and kG3*Pa_out<kL3*cB/c1:
    lambda3=0

HTU3=Hg3+lambda3*Hl3

#CHECK WETTING RATE3
if wr3<mwr3:
    print('Lp3 is too low')

#CHECK DIAMETER RATIO3
if Dc3/(riemp3.size/1000)<10:
    print('Packing Size is too small')

#CALCOLO NTU3
NTU3            = np.log(yin_3/yout_3) #supponiamo che reagisca tutto nel liquido
Z3              = NTU3*HTU3
print('NTU3 =','%.2f' %(NTU3))
print('HTU3 =','%.2f' %(HTU3*1000),'mm')
print('Z3 =','%.0f' %(Z3*1000),'mm')

if GPDC3==GPDC_4mmH2O_r:
    deltaP3=4
elif GPDC3==GPDC_8mmH2O_r:
    deltaP3=8
elif GPDC3==GPDC_21mmH2O_r:
    deltaP3=21
elif GPDC3==GPDC_42mmH2O_r:
    deltaP3=42
elif GPDC3==GPDC_83mmH2O_r:
    deltaP3=83

CorrectionP3=(usg3/usg3_)**1.75
deltaP_colonna3=deltaP3*Z3*CorrectionP3 #mmH2O
print()
print('DeltaP column 3 =','%.2f' %(deltaP_colonna3/100) ,'kPa')

#Dimensionamento Bocchelli sia in ingresso che in uscita in quanto la portata cambia di poco

D_bocchelli3=0.5 #metri
V_Bocchelli3=(Gv/3600)/(0.25*np.pi*D_bocchelli3**2) #m/s
if V_Bocchelli3>20 or V_Bocchelli3<10:
    print('Modificare diametro bocchelli gas')

delta_P3_bocchelli=air.rho*V_Bocchelli3**2

# dimensionamento distributore liquido
l3=Lv3/3600
DripDensity3=200
Nfori3=DripDensity3*Ac3
diam_fori3=0.006 
A_fori3=Nfori3*0.25*np.pi*diam_fori3**2
Vfori3=l3/A_fori3
battente3=(l3/A_fori3)**2/(2*9.81) #metri battente nella canalina del distributore

#bocchelli entrata liquido
diam_bocchello_L3in=0.14 #
A_bocchello_L3in=0.25*np.pi*diam_bocchello_L3in**2
Vbocchello_L3in=(Lv3/3600)/A_bocchello_L3in

if Vbocchello_L3in<0.8 or Vbocchello_L3in>1.2:
    print('Modificare diam bocchelli liquido in')

#bocchelli uscita liquido
diam_bocchello_L3out=0.2
A_bocchello_L3out=0.25*np.pi*diam_bocchello_L3out**2
Vbocchello_L3out=(Lv3/3600)/A_bocchello_L3out

if Vbocchello_L3out<0.3 or Vbocchello_L3out>0.5:
    print('Modificare diam bocchelli liquido out')


#Demister: tipo di demister T-01-P *2 https://vff.com/fileadmin/user_upload/Downloads/210323_RZ_Demister_E_WEB.pdf
#standard Sulzer KnitMesh demister polypropylene 9048 https://www.sulzer.com/-/media/files/products/separation-technology/feed-inlet-devices/gas_liquid_separation_technology.ashx
deltaP_demister3=dP_demister_dry_Setekleiv_Svendsen(S=360, voidage=.97, vs=usg3, rho=air.rho, mu=air.mu, L=0.2)

efficiency_demister3=separation_demister_ElDessouky(usg3, 0.97, 0.0002, 0.005)

deltaP3tot=deltaP_colonna3*10+delta_P3_bocchelli+deltaP_demister3*2
print()
print('deltaP tot Column 3 =', '%.3f' %(deltaP3tot/1000), 'kPa')


#dimensionamento fondo colonna
hfondo3=0.8 #metri di battende di liquido fondo colonna
kgsolfondocolonna3= Ac3*Solution.rhol*hfondo3
T_residenza3=(kgsolfondocolonna3)/(Lp3/60)

if T_residenza3<1 or T_residenza3>8:
    print('Cambiare tempo residenza sul fondo colonna')

############################################ VERIFICA TERMICA ###########################################################################
HeatReaction3   = nHCl_reacted*57 #57 sono kj/mol
delta_T3        = HeatReaction3/(Lp3*Solution.Cpl/1000)

q3=Yin_3-Ls3/Gs*Xout_3
retta_op3finale=lambda x: q3+(Ls3/Gs)*x

plt.figure(1)
# Grafici ottenuti dai dati del Perry8, i punti sono stati messi poi su MATLAB per ottenere delle curve di solubilità
# Sono riportate sia sia in x,y che in log10(x),log10(y)
# SONO TUTTI GRAFICI IN FRAZIONI MOLARI
plt.subplot(1,2,1)
plt.grid()
x=np.linspace(0.04,0.3,500)
plt.plot(x,Solubility30(x),color='r')
plt.plot(x,Solubility20(x),color='b')
plt.axhline(yout,linestyle='--',color='k',linewidth='0.8')
plt.axhline(yin,linestyle='--',color='k',linewidth='0.8')
# plt.axhline(yout/15,linestyle='--',color='k',linewidth='0.8')
plt.xlim(0,0.15)
plt.ylim(0,0.004)
plt.xlabel('$x_{HCl} $')
plt.ylabel('$y_{HCl} $')
plt.text(0.01,yout+0.00005,'150 ppmv',rotation=0)
plt.text(0.01,yin+0.00005,'5g HCl/$Nm^3$ - 3073.7 ppmv',rotation=0)
plt.legend(['Solubilità 30°C'],loc='upper left')

plt.subplot(1,2,2)
plt.grid(which='major')
plt.grid(which='minor',ls='-',linewidth=0.4)
plt.loglog(x,Solubility30(x),color='r')
plt.axhline(yout,linestyle='--',color='k',linewidth='0.8')
plt.axhline(yin,linestyle='--',color='k',linewidth='0.8')
plt.axhline(yout/15,linestyle='--',color='k',linewidth='0.8')
plt.xlim(0.04,0.2)
plt.ylim(1e-6,0.1)
plt.xlabel('$x_{HCl} $')
plt.ylabel('$y_{HCl} $')
plt.text(0.042,yout+0.00005,'150 ppmv',rotation=0)
plt.text(0.042,yin+0.00005,'5g HCl/$Nm^3$ - 3073.7 ppmv',rotation=0)
plt.text(0.1,yout/15+0.0000005,'10 ppmv',rotation=0)
plt.legend(['Solubilità 30°C'],loc='upper left')

plt.figure(2)
#Nei grafici di figura 2 vediamo sia in rapporti molari che in frazioni molari le curve di equilibrio a 20C e 30C
#e le rette operative nel caso 1 (SENZA RICICLO) che nel caso 2(MAX RICICLO)
# PRIMO GRAFICO ##################################################################################################################
plt.subplot(1,2,1)
plt.grid()
plt.plot(ratio(x),Solubility30r(x),color='r',label='Solubilità 30°C')

#SENZA RICICLO
x_=np.linspace(0,Xout0,500)
plt.plot(x_,op_line0(x_),'g',label='Curva Operativa (senza riciclo)')

#RICICLO
x_=np.linspace(Xin_max,Xout0,1000)
plt.plot(x_,op_line_maxRecycle(x_),'b',label='Curva Operativa (riciclo massimo)')

plt.plot([Xout_max],[Yin],'o',label='$X_{OUT,max}$')
plt.plot([Xin_max],[Yout],'o',label='$X_{IN,max}$')
plt.axhline(Yin, linestyle='-.',color='k',linewidth='0.7')
plt.axhline(Yout, linestyle='-.',color='k',linewidth='0.7')
plt.xlim(0,0.18)
plt.ylim(0,0.004)
plt.xlabel('$X_{HCl}$')
plt.ylabel('$Y_{HCl}$')
plt.title('Rapporti Molari')
plt.legend()

# SECONDO GRAFICO ##################################################################################################################
plt.subplot(1,2,2)
plt.grid()
plt.plot(x,Solubility30(x),color='r',label='Solubilità 30°C')

#SENZA RICICLO
x_=np.linspace(0,Xout0,1000)
plt.plot(frac(x_),op_line0(x_)/(1+op_line0(x_)),'g',label='Curva Operativa (senza riciclo)')

#RICICLO
x_=np.linspace(Xin_max,Xout0,500)
plt.plot(frac(x_),op_line_maxRecycle(x_)/(1+op_line_maxRecycle(x_)),'b',label='Curva Operativa (riciclo massimo)')

plt.plot([xout_max0],[yin],'o',label='$x_{OUT,max}$')
plt.plot([xin_max],[yout],'o',label='$x_{IN,max}$')
plt.axhline(yout, linestyle='-.',color='k',linewidth='0.7')
plt.axhline(yin, linestyle='-.',color='k',linewidth='0.7')
plt.xlim(0,0.18)
plt.ylim(0,0.004)
plt.xlabel('$x_{HCl}$')
plt.ylabel('$y_{HCl}$')
plt.text(0.135,yout+0.00005,'150 ppmv',rotation=0)
plt.text(0.02,yin+0.00005,'5g HCl/$Nm^3$ - 3073.7 ppmv',rotation=0)
plt.title('Frazioni Molari')
plt.legend(loc='upper left')

plt.figure(3)
vecprova=np.linspace(5e-3,2,1001)
plt.semilogx(vecprova,GPDC_125mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_83mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_42mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_21mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_8mmH2O_s(vecprova))
plt.axvline(flg1,color='r',ls='--',linewidth='0.6')
plt.axhline(CP1,color='r',ls='--',linewidth='0.6')
plt.legend(['$125mmH_{2}O/m$','$83mmH_{2}O/m$','$42mmH_{2}O/m$','$21mmH_{2}O/m$','$8mmH_{2}O/m$'])
plt.xlabel('$F_{LG}$')
plt.ylabel('CP')
plt.grid(which='major')
plt.grid(which='minor',ls='-',linewidth=0.4)
plt.title('Kister & Gill GPDC per riempimenti strutturati')
plt.suptitle('COLONNA 1')

plt.figure(4)
vecprova=np.linspace(5e-3,2,1001)
plt.semilogx(vecprova,GPDC_125mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_83mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_42mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_21mmH2O_s(vecprova))
plt.semilogx(vecprova,GPDC_8mmH2O_s(vecprova))
plt.axvline(flg2,color='r',ls='--',linewidth='0.6')
plt.axhline(CP2,color='r',ls='--',linewidth='0.6')
plt.legend(['$125mmH_{2}O/m$','$83mmH_{2}O/m$','$42mmH_{2}O/m$','$21mmH_{2}O/m$','$8mmH_{2}O/m$'])
plt.xlabel('$F_{LG}$')
plt.ylabel('CP')
plt.grid(which='major')
plt.grid(which='minor',ls='-',linewidth=0.4)
plt.title('Kister & Gill GPDC per riempimenti strutturati')
plt.suptitle('COLONNA 2')


plt.figure(5)
plt.grid()

xx=np.linspace(0.04,0.2,500)
plt.plot(ratio(xx),Solubility30r(xx),color='r',label='Curva di Solubilità 30°C')

x=np.linspace(0,Xout_C1,100)
plt.plot(x,rettaGlobal(x),color='b',label='Retta Operativa Globale')

x3=np.linspace(Xin_C2,Xout_C2,100)
plt.plot(x3,retta_C2(x3),color='g',label='Curva Operativa C2')

x4=np.linspace(Xin_C1,Xout_C1,100)
plt.plot(x4,retta_C1(x4),color='k',label='Curva Operativa C1')

plt.plot([Xout_C2],[Yin_C2],'o',color='r')
plt.plot([Xout_C1],[Yin_C1],'o',color='r')

plt.axhline(Yout,linestyle='--',color='k',linewidth='0.7')
plt.axhline(Yin_C1,linestyle='--',color='k',linewidth='0.7')
plt.axhline(Yin_C2,linestyle='--',color='k',linewidth='0.7')
plt.legend(loc='upper left')
plt.xlabel('$X_{HCl}$')
plt.ylabel('$Y_{HCl}$')
plt.xlim(0,0.16)
plt.ylim(0,0.004)

plt.figure(6)
plt.grid()

xx=np.linspace(0.04,0.2,500)
plt.plot(xx,Solubility30(xx),color='r')

xx=np.linspace(0.04,0.2,500)
plt.plot(xx,0.5*(Solubility20(xx)+Solubility30(xx)),color='r')


x=np.linspace(0,Xout_C2,100)
plt.plot(frac(x),(rettaGlobal((x)))/(1+rettaGlobal(x)),color='c')

x2=np.linspace(Xout_C2,Xout_C1,100)
plt.plot(frac(x2),(rettaGlobal(x2))/(1+rettaGlobal(x2)),color='b')

x3=np.linspace(Xin_C2,Xout_C2,100)
plt.plot(frac(x3),retta_C2(x3)/(1+retta_C2(x3)),color='g')

x4=np.linspace(Xin_C1,Xout_C1,100)
plt.plot(frac(x4),retta_C1(x4)/(1+retta_C1(x4)),color='k')

plt.plot([xout_c2],[yin_c2],'o',color='r')
plt.plot([xout_c1],[yin],'o',color='r')

plt.axhline(yout,linestyle='--',color='k',linewidth='0.7')
plt.axhline(yin,linestyle='--',color='k',linewidth='0.7')
plt.axhline(yin_c2,linestyle='--',color='k',linewidth='0.7')

plt.xlim(0,0.16)
plt.ylim(0,0.004)

plt.xlabel('$x_{HCl}$')
plt.ylabel('$y_{HCl}$')

plt.figure(7)
x=np.linspace(0.04,0.15,501)
plt.plot(ratio(x),Solubility30r(x))
# plt.plot([Xout_3],[Yin_3],'o')
# plt.plot(x,retta_op3finale(x))
plt.axhline(yin_3,linestyle='--',color='k')
plt.axhline(yout_3,linestyle='--',color='k')
plt.grid()
plt.ylim(0,0.0002)
plt.xlim(0,0.1)

plt.figure(8)
vecprova3=np.linspace(2e-2,3.5,1001)
plt.loglog(vecprova3,GPDC_flooding(vecprova3),'--')
plt.loglog(vecprova3,GPDC_83mmH2O_r(vecprova3))
plt.loglog(vecprova3,GPDC_42mmH2O_r(vecprova3))
plt.loglog(vecprova3,GPDC_21mmH2O_r(vecprova3))
plt.loglog(vecprova3,GPDC_8mmH2O_r(vecprova3))
plt.loglog(vecprova3,GPDC_4mmH2O_r(vecprova3))
plt.axvline(flg3,color='r',ls='--',linewidth='0.6')
plt.axhline(K4,color='r',ls='--',linewidth='0.6')
plt.legend(['Flooding','$83mmH_{2}O/m$','$42mmH_{2}O/m$','$21mmH_{2}O/m$','$8mmH_{2}O/m$','$4mmH_{2}O/m$'])
plt.xlabel('$F_{LG}$')
plt.ylabel('K4')
plt.grid(which='major')
plt.grid(which='minor',ls='-',linewidth=0.4)
plt.title('Norton Co. GPDC per riempimenti Random')
plt.suptitle('COLONNA DI ASSORBIMENTO CHIMICO')





