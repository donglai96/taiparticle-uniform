import pickle
def save_input(id,
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
               p0,
               k0):

    dict_input = dict({'id':id,
               'B0':B0, 
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
               'wave_distribution':wave_distribution,
               'p0':p0,
               'k0':k0})
    with open(id + '/' + 'input','wb') as f:
        pickle.dump(dict_input, f)
    with open(id + '/' + 'input.txt','w') as f:
        f.write(str(dict))


    return dict_input
