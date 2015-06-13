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

def tpss(phi,visc,ct,area,tDA,k):
    ''' Time required to reach  pseudosteady state flow in hours'''
    return (phi*visc*ct*area*tDA)/(0.0002637*k)

def tDA(phi,visc,ct,area,t,k):
    return (0.0002637*k*t/(visc*phi*ct*area))

def cfd(permMatrix,permFrac,widthFrac,halfLength):
    ''' Calcules fracture dimensionless conductivity (Cfd)
    permMatrix is the matrix permeability, md
    permFrac is the fracture permeability, md
    widthFrac is the fracture width, ft
    halfLength is the fracture half-length in ft
    '''
    return (widthFrac*permFrac)/(permMatrix*halfLength)

def F(cfd):
    ''' Calculates the Cinco Ley and Samaniego F
    cfd is teh dimensionless fracture conductivity
    '''
    u = math.log(cfd)
    return (1.65-0.328*u+0.116*u**2)/(1+0.18*u+0.064*u**2+0.005*u**3)

def sf(f,xf,rw):
    ''' Calculates pseudoskin
    f from Cinco-Ley & Samaniego plot or F function
    xf fracture half length, ft
    rw wellbore radius, ft
    '''
    return f-math.log(xf/rw)
    
def rwprime(rw,sf):
    ''' Calculates the effective wellbore radius
    rw is the wellbore radius
    sf is the fracture pseudoskin
    '''
    return rw*math.exp(-sf)


    


