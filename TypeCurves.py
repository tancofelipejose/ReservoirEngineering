# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 20:33:02 2015

@author: Leytzher Muro, DVM/11
"""

import math
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

      

