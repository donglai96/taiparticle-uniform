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
from wave_generate import Waves_generate
from wave import Wave
import matplotlib.pyplot as plt

###################################################
# parameters
id = 'test_xtao'
print('Welcom and start!')
print('id', id)
ti.init(arch = ti.cpu,default_fp=ti.f64)

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
record_num = 25

# time step (in gyroperiod)
dt_num = 0.02

# resonant wave frequency (in gyrofrequency)
w_res_num = 0.3

# lowercutoff and uppercutoff
w_lc_num = 0.2
w_uc_num = 0.4
# wave frequency width
w_width_num = 0.999

# number of wave frequency 
nw = 100

# z range in units of res wave
dz_num = 200000
# wave amplitude in Gauss
Bw = 1e-7

# init mass and charge
mass = cst.Me
charge = cst.Charge * -1 # -1 is electron

# wave direction
direction = 1


# initial energy in eV
# end of parameter settings
###################################################

print('Calculating resonance energy')
wce = gyrofrequency(cst.Charge,cst.Me,B0)
alpha = np.deg2rad(pitch_angle_degree) # pitch angle is here!
w = w_res_num * wce
w_lc = w_lc_num * wce
w_uc = w_uc_num * wce
w_width = w_width_num * wce


p0,k0 = get_resonance_p_whistler(w,wce,n0,alpha,nres = -1)
print('momentum of resonating particles are:',p0)
print('resonant wave number is :', k0)
#erg_particle  = p2e(p)
print('resonating frequency is (Hz) ', w/(2 * np.pi))

print('E0 is ', erg2ev(p2e(p0))/1000, ' keV')


# functions that going to be used in taichi frame
#print('In xtao code')
#print(3.1415926 * 2 / (cst.Charge * B0/(cst.Me * cst.C)))
gamma = (1 + p0**2 / (cst.Me**2*cst.C**2))**0.5; 

wce_rel = wce/gamma
T_gyro = 2 * np.pi/ wce_rel
dt = T_gyro * dt_num
Nt = int(t_total_num/dt_num)
print(' time step is:', dt)
print(' mass is ', mass)
print(' charge is ', charge)
print(' total time is', t_total_num * T_gyro)
print(' total time step is ', Nt)

###################################################

# init the taichi 
# note: if init with gpu with float 64 might not be working 

particles = Particle.field(shape = (Np,))
###################################################
# init the particle phase and location

# 1. get the wavelenght
# 2. z range (in units of the wave length of iw at iw_res) = 200000
# 
if nw > 1:

    dw = (w_uc - w_lc) / (nw - 1)
    ws = np.array([i * dw for i in range(nw)] ) + w_lc
    
else:
    ws = w
print('ws')
print(ws)
k_res = k0
wave_length = 2 * np.pi / k_res
if Np > 1:
    dz = dz_num * wave_length/(Np - 1)
    
else:
    dz = dz_num * wave_length

dphi = 2 * np.pi / Np

pperp = p0 * np.sin(alpha)

px_numpy = np.zeros(Np)
py_numpy = np.zeros(Np)
pz_numpy = np.zeros(Np)

for n in range(Np):
    phi = dphi * n
    px_numpy[n] = pperp * np.cos(phi)
    py_numpy[n] = pperp * np.sin(phi)
    pz_numpy[n] = p0 * np.cos(alpha)

px_init = ti.field(dtype = ti.f64,shape = (Np,))
py_init= ti.field(dtype = ti.f64,shape = (Np,))
pz_init= ti.field(dtype = ti.f64,shape = (Np,))

px_init.from_numpy(px_numpy)
py_init.from_numpy(py_numpy)
pz_init.from_numpy(pz_numpy)
###################################################
# init the wave info use numpy

waves_init = Waves_generate(direction = direction,
                            frequencies = ws,
                            B0 = B0,
                            ne = n0,
                            Bw = Bw,
                            w_width =w_width,
                            distribution ="Gaussian")
waves_init.generate_parallel_wave()
ws_taichi = ti.field(dtype = ti.f64,shape = (nw,))
phi0_taichi = ti.field(dtype = ti.f64,shape = (nw,))
Ewx_taichi = ti.field(dtype = ti.f64,shape = (nw,))
Bwy_taichi = ti.field(dtype = ti.f64,shape = (nw,))
k_taichi = ti.field(dtype = ti.f64,shape = (nw,))

ws_taichi.from_numpy(waves_init.ws)
phi0_taichi.from_numpy(waves_init.phi0)
Ewx_taichi.from_numpy(waves_init.Ewx)
Bwy_taichi.from_numpy(waves_init.Bwy)
k_taichi.from_numpy(waves_init.k)

###################################################
# init wave in taichi scope

waves = Wave.field(shape = (nw,))
# ts  =  ti.field(ti.f64, shape=())
dt_taichi = ti.field(ti.f64, shape=())
dt_taichi[None] = dt
###################################################

# init function
@ti.kernel
def init():
    for n in range(Np):
        particles[n].initParticles(mass, charge)
        particles[n].initPos(0.0,0.0,dz * n)
        particles[n].initMomentum(px_init[n], py_init[n],pz_init[n])
    for m in range(nw):
        waves[m].initialize(ws_taichi[m],phi0_taichi[m],Ewx_taichi[m],
                            Bwy_taichi[m],k_taichi[m])
        
@ti.kernel
def simulate():
    # one step

    for n in range(Np):
        # get field
        B =ti.Vector([0.0,0.0,B0])
        E = ti.Vector([0.0,0.0,0.0])
        for m in range(nw):
            
            waves[m].get_wavefield(particles[n].r, particles[n].t)
            B += waves[m].Bw
            E += waves[m].Ew
        
        particles[n].t += dt_taichi[None]
        particles[n].leap_frog(dt_taichi[None],E,B)
    #print('B')
        #particles[n].leap_frog(dt,E,B)
    # ts += dt
    # print('ts', ts)

###################################################
# Begin of the main loop
init()
p_record = np.zeros((Nt//record_num,Np,3))
r_record = np.zeros((Nt//record_num,Np,3))
for t_num in range(Nt):
    
    simulate()
    if t_num % record_num ==0:
        for n in range(Np):
            p_record[t_num //record_num,n,:] = particles[n].p.to_numpy()
                #print('this is p',particles[n].p.to_numpy())
            r_record[t_num //record_num,n,:] = particles[n].r.to_numpy()
print('finished')
# End of the main loop
###################################################


###################################################

# plot
# time , particles of index
px_result = p_record[:,:,0]  
py_result = p_record[:,:,1] 
pz_result = p_record[:,:,2] 
p_result = np.sqrt(px_result **2 + py_result **2 + pz_result **2)
energy_result = erg2ev(p2e(p_result))/1000
energy_0 = energy_result[0,:]
delta_energy = np.average((energy_result - energy_result[0,:])**2,axis = 1)
print(px_result.shape)
plt.plot(delta_energy)
plt.show()
#plt.plot()
###################################################
