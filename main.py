"""
@author Donglai Ma
@email donglaima96@gmail.com
@create date 2023-01-25 17:07:58
@modify date 2023-01-25 17:07:58
@desc taichi test particle of 
uniform background magnetic field
"""



import numpy as np
from sympy import Naturals0
import taichi as ti
import constants as cst
from res_energy import * 
from particle import Particle



###################################################
# parameters
id = 'test_xtao'
print('Welcom and start!')
print('id', id)
# background magnetic field strength (Gauss)
B0 = 0.0014

# density in cm^-3
n0 = 10

# number of charged particles
Np = 400

# initial pitch angle
pitch_angle_degree = 140 # to the background magnetic field!


# length of run (in gyroperiod)
t_total_num = 200

# record x (how many step make a record)
record_num = 10

# time step (in gyroperiod)
dt_num = 0.02

# resonant wave frequency (in gyrofrequency)
w_res = 0.3

# lowercutoff and uppercutoff
w_lc = 0.2
w_uc = 0.4

# number of wave frequency 
nw = 100

# wave amplitude in nT
Bw = 0.01

# initial energy in eV

print('Calculating resonance energy')
wce = gyrofrequency(cst.Charge,cst.Me,B0)
alpha = np.deg2rad(pitch_angle_degree)
w = w_res * wce
p0 = get_resonance_p_whistler(w,wce,n0,alpha,nres = -1)
print('momentum of resonating particles are:',p0)

#erg_particle  = p2e(p)
print('resonating frequency is (Hz) ', w/(2 * np.pi))

print('E0 is ', erg2ev(p2e(p0))/1000, ' keV')

gamma = (1 + p0**2 / (cst.Me**2*cst.C**2))**0.5; 
# end of parameter settings
###################################################
# functions that going to be used in taichi frame
wce = gyrofrequency(cst.Charge,cst.Me,B0)
#print('In xtao code')
#print(3.1415926 * 2 / (cst.Charge * B0/(cst.Me * cst.C)))

wce_rel = wce/gamma
T_gyro = 2 * np.pi/ wce_rel
dt = T_gyro * dt_num
print(' time step is:', dt)
print(' total time is', t_total_num * T_gyro)

# init the taichi 
# note: if init with gpu with float 64 might not be working 
ti.init(arch = ti.cpu,default_fp=ti.f64)

particles = Particle.field(shape = (Np,))
###################################################
# init the particle phase

###################################################
# Begin of the main loop


# End of the main loop
###################################################


###################################################

# plot

###################################################
