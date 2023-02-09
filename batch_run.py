import numpy as np
from generate_input import save_input
pitch_angles = np.array([0.1,1,5,10,20,30,40,50,60,70,80,85,89,89.9])
pitch_angle_inputs = 180 - np.array([0.1,1,5,10,20,30,40,50,60,70,80,85,89,89.9])
B0 = 0.0014

Bw_B0 = np.array([10**-5,10**-4.5,10**-4,10**-3.5,10**-3,10**-2.5,10**-2,10**-1.5,10**-1,10**-0.5])
Bws = B0 * Bw_B0
for i,Bw in enumerate(Bws):
    print(i)
    Bw_round = format(Bw,'.4e')
    for j, pitch_angle in enumerate(pitch_angle_inputs):
        folder_name = 'oliver_results'
        id = 'Bw_' + str(i) + 'th_' + 'alpha_' + str(j) + 'th' 
        save_input(folder_name, id, 
                    B0 = 0.0014,
                    n0 = 10,
                    Np = 120,
                    pitch_angle_degree=pitch_angle,
                    t_total_num=1000,
                    record_num=250,
                    dt_num = 0.002,
                    w_res_num=0.3,
                    w_lc_num = 0.2,
                    w_uc_num=0.4,
                    w_width_num = 0.999,
                    nw = 100,
                    dz_num = 200000,
                    Bw = Bw_round,
                    wave_distribution='Constant')


import subprocess

for i,Bw in enumerate(Bws):
    for j, pitch_angle in enumerate(pitch_angle_inputs):
        folder_name = 'oliver_results'
        id = 'Bw_' + str(i) + 'th_' + 'alpha_' + str(j) + 'th' 
        run_arg = folder_name + '/' + id
        subprocess.run(["python", "main.py", run_arg])