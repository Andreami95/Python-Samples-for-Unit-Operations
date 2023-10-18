# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:25:29 2022

@author: Andrea
"""

#TUBI IN ACCIAIO INOX SENZA SALDATURA E SALDATE A NORME ANSI TECNOSTEEL
# le misure sono tutte in mm

#DEFINIZIONE CLASSE RIEMPIMENTO LISTE
class Size:
    def __init__(self, d, t, w): 

        self.d=d #mm
        self.t=t #mm
        self.w=w #kg/m

Steel_list = []


Steel_list.append( Size(10.3,[1.24,1.73,2.41],[0.290,0.360,0.460]))
Steel_list.append( Size(13.7,[1.65,2.24,3.02],[0.500,0.630,0.800]))
Steel_list.append( Size(17.2,[1.65,2.31,3.20],[0.640,0.850,1.100]))
Steel_list.append( Size(21.3,[1.65,2.11,2.77,3.73],[0.820,1.020,1.270,1.680]))
Steel_list.append( Size(26.7,[1.65,2.11,2.77,3.91],[1.050,1.310,1.680,2.190]))
Steel_list.append( Size(33.4,[1.65,2.77,3.38,4.55],[1.320,2.130,2.500,3.230]))
Steel_list.append( Size(42.2,[1.65,2.77,3.56,4.85],[1.690,2.750,3.380,4.470]))
Steel_list.append( Size(48.3,[1.65,2.77,3.68,5.08],[1.940,3.170,4.050,5.400]))
Steel_list.append( Size(60.3,[1.65,2.77,3.91,5.54],[2.450,3.990,5.440,7.480]))
Steel_list.append( Size(73.0,[2.11,3.05,5.16,7.01],[3.740,5.360,8.620,11.410]))
Steel_list.append( Size(88.9,[2.11,3.05,5.49,7.62],[4.510,6.580,11.280,15.250]))
Steel_list.append( Size(101.6,[2.11,3.05,5.74,8.08],[5.170,7.550,13.570,18.620]))
Steel_list.append( Size(114.3,[2.11,3.05,6.02,8.56],[5.840,8.520,16.070,22.310,28.270,33.540,41.200]))
Steel_list.append( Size(141.3,[2.77,3.40,6.55,9.52],[10.000,11.800,21.760,30.920]))
Steel_list.append( Size(168.3,[2,77,3.40,7.11,10.97],[11.310,14.110,28.260,42.560]))
Steel_list.append( Size(219.1,[2.77,3.76,8.18,12.70],[15.500,19.900,42.530,64.630]))
Steel_list.append( Size(273,[3.40,4.19,8.80,9.27,12.70],[25.900,30.850,58.200,60.290,81.460]))
Steel_list.append( Size(323.8,[3.96,4.57,9.52,12.70],[31.800,35.500,73.820,97.400]))
Steel_list.append( Size(355.6,[3.96,4.76,11.12,19.05],[34.990,42.100,94.300,157.940]))
Steel_list.append( Size(406.4,[4.19,4.76,9.52,12.70],[42.530,48.220,93.210,123.180]))
Steel_list.append( Size(457.2,[4.19,4.76,9.52,12.70],[47.700,54.300,105.140,139.070]))
Steel_list.append( Size(508,[4.76,5.54,9.52,12.70],[60.320,69.800,117.070,156.100]))
Steel_list.append( Size(558.8,[59.52,12.70],[129.010,171.010]))
Steel_list.append( Size(609.6,[5.54,6.35,9.52,12.70],[84.100,94.450,140.940,186.920]))
Steel_list.append( Size(660.4,[9.52,12.70],[152.87,202.83]))
Steel_list.append( Size(711.2,[9.52,12.70],[164.65,218.54]))
Steel_list.append( Size(762.0,[6.35,7.92,9.52,12.70],[111.200,136.370,164.340,218.730]))
Steel_list.append( Size(812.8,[9.52,12.70],[188.660,250.550]))
Steel_list.append( Size(863.6,[9.52,12.70],[200.590,266.460]))
Steel_list.append( Size(914.4,[9.52,12.70],[212.520,282.360]))

DN6i=Steel_list[0]
DN8i=Steel_list[1]
DN10i=Steel_list[2]
DN15i=Steel_list[3]
DN20i=Steel_list[4]
DN25i=Steel_list[5]
DN32i=Steel_list[6]
DN40i=Steel_list[7]
DN50i=Steel_list[8]
DN65i=Steel_list[9]
DN80i=Steel_list[10]
DN90i=Steel_list[11]
DN100i=Steel_list[12]
DN125i=Steel_list[13]
DN150i=Steel_list[14]
DN200i=Steel_list[15]
DN250i=Steel_list[16]
DN300i=Steel_list[17]
DN350i=Steel_list[18]
DN400i=Steel_list[19]
DN450i=Steel_list[20]
DN500i=Steel_list[21]
DN550i=Steel_list[22]
DN600i=Steel_list[23]
DN650i=Steel_list[24]
DN700i=Steel_list[25]
DN750i=Steel_list[26]
DN800i=Steel_list[27]
DN850i=Steel_list[28]
DN900i=Steel_list[29]
