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
def tpss(phi,visc,ct,area,tDA,k):
    ''' Time required to reach  pseudosteady state flow in hours'''
    return (phi*visc*ct*area*tDA)/(0.0002637*k)

def tDA(phi,visc,ct,area,t,k):
    return (0.0002637*k*t/(visc*phi*ct*area))




