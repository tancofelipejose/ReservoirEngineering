# -*- coding: utf-8 -*-
"""
    gasPVT.py - Reservoir Engineering Toolbox
        Copyright (C) 2015  Leytzher Muro

    This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
                (at your option) any later version.

    This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import math

def Ppc(spGr):
    return 756.8-131.07*spGr-3.6*spGr**2

def Tpc(spGr):
    return 169.2+349.5*spGr-74.0*spGr**2

def Pr(P,Ppc):
    return P/Ppc

def Tr(T,Tpc):
    return T/Tpc
    

def gasDensity(P,T,spGr):
    ''' pressure in psia, T in Rankine
    returns gas Density in g/cc'''
    Mg = 28.967*spGr 
    pseudoReducedPress = Pr(P,Ppc(spGr))
    pseudoReducedTemp = Tr(T,Tpc(spGr))
    z = zFactor(pseudoReducedPress,pseudoReducedTemp)
    return 0.00149406*P*Mg/(z*T)


def gasViscosity(P,T,spGr):
    ''' Lee et al. estimation of natural gas viscosity 
    P in psia, T in Rankine
    '''
    Mg = 28.967*spGr
    rhog = gasDensity(P,T,spGr)
    K1 = ((0.00094+(2e-6)*Mg)*T**1.5)/(209+19*Mg+T)
    X = 3.5+(986/T)+(0.01*Mg)
    Y = 2.4-0.2*X
    return K1*math.exp(X*rhog**Y)

def gasFVF(P,T,z,flag):
    '''
    
    flag:
    1 = Bg in rcf/scf
    2 = Bg in RB/scf
    3 = Bg in Rm3/Sm3

    Important:
    
    if flag=1 or flag=2 then P in psia and T in R
    if flag=3 then P is in Kpa and T in K
    '''
    if flag==1:
        return 0.0282793*z*T/P
    elif flag==2:
        return 0.00503676*z*T/P
    elif flag==3:
        return 0.350958*z*T/P


def zFactor(Pr,Tr):
    '''Dranchuk and Abou-Kassem fit of Standing and Katz gas compressibility factor '''
    A1 = 0.3265
    A2 = -1.0700
    A3 = -0.5339
    A4 = 0.01569
    A5 = -0.05165
    A6 = 0.5475
    A7 = -0.7361
    A8 = 0.1844
    A9 = 0.1056
    A10 = 0.6134
    A11 = 0.7210

    zguess=1.0
    delta=1.0
    while (delta>0.0001):
        rhor = 0.27*Pr/(zguess*Tr)
        z = 1+(A1+(A2/Tr)+(A3/Tr**3)+(A4/Tr**4)+(A5/Tr**5))*rhor+((A6+(A7/Tr)+(A8/Tr**2))*rhor**2)-(A9*((A7/Tr)+(A8/Tr**2))*rhor**5)+A10*(1+(A11*rhor**2))*((rhor**2)/(Tr**3))*(math.exp(-A11*(rhor*2)))
        delta=abs(z-zguess)
        zguess=z
    return z


def waterCompressibility(P,T,salinity):
    ''' Water compressibility in psi**-1 (Osif's laboratory measurements)
    P = pressure, psi
    T = temperature, F
    salinity in g/l of solution (good for 0 to 200 g/l NaCl)
    '''
    m1 = 7.033
    m2 = 541.5
    m3 = -537
    m4 = 404.4e3
    cwInv = m1*P+m2*salinity+m3*T+m4
    return 1/cwInv
    

def rockCompressibility(porosity):
    ''' Rock compressibility according to Hall's correlation (1/psi)
    '''
    return 1.87e-6*(porosity)**-0.415

