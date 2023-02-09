import numpy as np
import constants as cst
import os

def save_input( folder,
               id,
               B0, 
               n0, 
               Np,
               pitch_angle_degree,
               t_total_num,
               record_num,
               dt_num,
               w_res_num,
               w_lc_num,
               w_uc_num,
               w_width_num,
               nw,
               dz_num,
               Bw,
               wave_distribution,
               ):
    print('generating input at folder: ', folder+ '/'  + id)
    isExist = os.path.exists(folder + '/' + id)
    if not isExist:

   # Create a new directory because it does not exist
        os.makedirs(folder + '/' + id)
    else:
        raise ValueError('Can not create directory, please change the id:')
    para = {'B0':B0,
            'n0':n0,
            'Np':Np,
            'pitch_angle_degree':pitch_angle_degree,
            't_total_num':t_total_num,
            'record_num':record_num,
            'dt_num':dt_num,
            'w_res_num':w_res_num,
            'w_lc_num':w_lc_num,
            'w_uc_num':w_uc_num,
            'w_width_num':w_width_num,
            'nw':nw,
            'dz_num':dz_num,
            'Bw':Bw,
            'wave_distribution':wave_distribution}

    with open(folder+ '/'  + id + '/' + 'input.txt', 'w') as f:
        for key, value in para.items():
            
            if key == 'B0':
                f.write('# background magnetic field strength (Gauss)\n')
            if key == 'n0':
                f.write('# density in cm^-3\n')
            if key == 'Np':
                f.write('# number of charged particles\n')
            if key == 'pitch_angle_degree':
                f.write('# initial pitchangle, in the range of [90, 180]. e.g. 160 -> 20, 170->10\n')
            if key == 't_total_num':
                f.write('# length of run in gyroperiod (relativistic)\n')
            if key == 'record_num':
                f.write('# how many step make one record\n')
            if key == 'dt_num':
                f.write('# timestep in gyroperiod\n')
            if key == 'w_res_num':
                f.write('# resonant frequency in gyrofrequency\n')
            if key == 'w_lc_num':
                f.write('# lower cutoff and upper cutoff\n')
            if key == 'w_width_num':
                f.write('# wave frequency width, this para only works when the wave distribution set as Guassian\n')
            if key =='nw':
                f.write('# The number of wave frequency\n')

            if key == 'dz_num':
                f.write('# z range in units of res wave length\n')
            if key =='Bw':
                f.write('# wave amplitude in Gauss (Bw,j = Bw / sqrt(nw))\n')
            if key =='wave_distribution':
                f.write('# wave distribution, support Constant and Guassian\n')
            
            f.write(f"{key}: {value}\n")

if __name__ == '__main__':
    folder = '.'
    id = 'test'
    # background magnetic field strength (Gauss)

    B0 = 0.0014

    # density in cm^-3
    n0 = 10

    # number of charged particles
    Np = 400

    # initial pitch angle
    pitch_angle_degree = 140 # to the background magnetic field!


    # length of run (in gyroperiod)
    t_total_num = 1000

    # record x (how many step make a record)
    record_num = 250

    # time step (in gyroperiod)
    dt_num = 0.002

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
    Bw = 1.4e-7

    # init mass and charge
    mass = cst.Me
    charge = cst.Charge * -1 # -1 is electron

    # wave direction
    direction = 1

    #wave_distribution = "Gaussian"
    wave_distribution = "Constant"
    save_input( folder,
               id,
               B0, 
               n0, 
               Np,
               pitch_angle_degree,
               t_total_num,
               record_num,
               dt_num,
               w_res_num,
               w_lc_num,
               w_uc_num,
               w_width_num,
               nw,
               dz_num,
               Bw,
               wave_distribution,
               )