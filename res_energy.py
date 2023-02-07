"""
@author Donglai Ma
@email donglaima96@gmail.com
@create date 2023-01-25 21:16:29
@modify date 2023-01-25 21:16:29
@desc get resonant energy for particles
"""

from signal import raise_signal
import numpy as np
import constants as cst
import sympy as sym
from taichiphysics import *



def get_resonance_velocity(w, kpara, alpha, wce, n):
    v = sym.Symbol('v')

    lhs = w - kpara * v * np.cos(alpha)

    rhs = n * -1 * wce*(1 - (v**2)/cst.C**2)**0.5


    # print('lhs is ', lhs)
    # print('rhs is ', rhs)
    print(" resonance n is ", n )
    # print("the solution of v  are")
    #print(sym.solve(lhs - rhs, v))
    return sym.solve(lhs - rhs, v)


def get_resonance_p_whistler(w, wce, n, alpha, nres):

    # get k     
    Q = cst.Charge
    M = cst.Me
    
    wpe = np.sqrt(4 * np.pi * n*Q**2/M)
    RR = 1 - wpe**2 / ((w-wce) * w)
    k = w *np.sqrt(RR) / cst.C
    v = get_resonance_velocity(w, k, alpha, wce, nres)[0]
    if v< 0 :
        print('Check the pitch angle!')
        print('Changing the initial solution')
        v = get_resonance_velocity(w, k, alpha, wce, nres)[1]
    gamma = 1 / (1 - (v**2)/cst.C**2)**0.5
    p = gamma * v * cst.Me
    return p, k 

# if __name__ == "__main__" :
#     n = 10
#     B = 0.0014 # Gauss
#     wce = gyrofrequency(cst.Charge,cst.Me,B)
#     alpha = np.deg2rad(140)
#     w = 0.3 * wce
#     p = get_resonance_p_whistler(w,wce,n,alpha,nres = -1)
#     print(p)
    
#     #erg_particle  = p2e(p)
#     print('resonating frequency is ', w/(2 * np.pi))

#     print(p2e(p)) 
#     print(erg2ev(p2e(p))/1000, ' keV')
#     #print(np.sqrt(p**2 * C**2 + E0**2) - E0 )
#     #print(erg2ev(p2e(p)))