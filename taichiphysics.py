import numpy as np
import taichi as ti
import constants as cst
C = 3e10
M = 9.1094e-28 
Q = 4.8032e-10
def gyrofrequency(q, m ,B):
    return q * B / (m * cst.C)


def plasmafrequency(q,m,n):
    return np.sqrt(4 * np.pi *n / m) *q




def ev2erg(ev):
    return ev * 1.60218e-12
def erg2ev(erg):
    return erg / 1.60218e-12
# for electrons
def e2p(e, E0=M*C**2):
    return (e * (e + 2 * E0))**0.5 / C 

def p2e(p, E0=M*C**2):
    return (p**2 * C**2 + E0**2)**0.5 - E0 

def p2v(p,m= M):
    gamma_m = m * (1 + p**2 /(m * m*C**2)) **0.5
    return p/gamma_m

def v2p(v,m = M):
    gamma_m = m * 1 / (1 - v**2 / C**2)*0.5
    return gamma_m * v