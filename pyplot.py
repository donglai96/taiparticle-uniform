import numpy as np
import matplotlib.pyplot as plt
import os
from taichiphysics import *

id = input("Please enter the id:\n")
 
print(f'You entered {id}')
p_r_name = id + '/p_r_phi.npy'
with open(p_r_name,'rb') as f:
    p_results = np.load(f)
    r_results = np.load(f)
    phi_results = np.load(f) # phi is the phi_particle - phi_wave


p_total_result = np.sqrt(p_results[:,:,0] **2 + p_results[:,:,1] **2 +p_results[:,:,2] **2)
energy_result = erg2ev(p2e(p_total_result))/1000 #keV

px = p_results[:,:,0] 
py = p_results[:,:,1]
pz = p_results[:,:,2]
rx = r_results[:,:,0]
ry = r_results[:,:,1]
rz = r_results[:,:,2]


pper_sq = px**2+py**2
psq = pper_sq+pz**2
p = np.sqrt(psq)
pitch_angle = 180 - np.rad2deg(np.arccos((pz/p)))

# delta energy
energy_0 = energy_result[0,:]
print(energy_0.shape)
delta_energy = np.average((energy_result - energy_result[0,:])**2,axis = 1)
# delta pitch angle
#pitch_angle_0 = pitch_angle[0,:]
delta_angle = np.average((pitch_angle - pitch_angle[0,:])**2,axis = 1)


record_gyro = 0.5 # i set all the record number is 0.5 T_gyro per record
ips = np.random.choice(120, 24)


time_total_step = p_total_result.shape[0]
fig,axs = plt.subplots(2,2,sharex = 'col')
tplot = record_gyro * np.arange(time_total_step )
for ip in ips:
    #axs[0,1].plot(timev[:], energy[ip,:])
    axs[0,0].plot(tplot, pitch_angle[:,ip])
    axs[0,1].plot(tplot, energy_result[:,ip])
axs[1,1].plot(tplot,delta_energy)
axs[1,0].plot(tplot,delta_angle)

axs[1,1].set_xlabel(r'$t/T_{gyro}$')
axs[1,0].set_xlabel(r'$t/T_{gyro}$')
axs[1,0].set_ylabel(r'$<\Delta \alpha^2> (degree^2)$')
axs[1,1].set_ylabel(r'$<\Delta E^2> (keV^2)$')
axs[0,1].set_ylabel(r'$keV$')

axs[0,0].set_ylabel(r'$degrees$')


plt.show()
