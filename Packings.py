# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:02:48 2022

@author: Andrea Milazzo
"""

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

# MELLAPAK ---------------------------------------------------------
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

Mellapak_Plastic_list.append (Packing('125y',125,'none','none','none','none','none','none'))
Mellapak_Plastic_list.append (Packing('250y',125,0.96,72,45,17,24.1,11.9))
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