"""
@author Donglai Ma
@email donglaima96@gmail.com
@create date 2023-01-25 17:10:34
@modify date 2023-01-25 17:10:34
@desc generate initial numpy array for the wave info
"""
import numpy as np
from scipy import integrate
import constants as cst
from taichiphysics import *
import random





class Waves_generate(object):
    """

    Because the operation of wave particle interaction is in taichi part,
    This code is just for initialization wave information

    """
    def __init__(self,direction,frequencies, B0, ne, Bw, w_width, distribution ="Gaussian"):
        self.ne = ne
        self.B0 = B0 # B0 is background and Bw is wave
        self.direction = direction # 1 or -1
        if direction == -1:
            print("Using anti-parallel propagation wave")
        if np.abs(direction)!= 1:
            raise ValueError("The direction has to be 1 (parallel) or -1(anitiparallel)")
        self.ws =frequencies
        self.w_width = w_width
        self.nw = self.ws.shape[0]
        print('total number of wave frequency is:',self.nw)

        if self.nw == 1:
            print('The input wave is monochromatic wave')
            print('The input frequency is :', self.ws)
            print('*************')
            print('Using the input magnetic field as field strength')
            self.Bwy = np.array([Bw])
        else:
            # using wave with a broadband
            if self.ws[-1] <= self.ws[0]:
                raise ValueError("The given frequencies should be from lower to higher!")
            
            Power_sum_square = 0
            self.Bwy = np.zeros(self.nw) 
            
            self.w_m = (self.ws[-1] + self.ws[0])/2
            
            # A Guassian distribution
            if distribution == "Gaussian":
                print("Using Guassian distribution of wave")
                for i in range(self.nw):
                    tmp = (self.ws[i] - self.w_m)/ self.w_width
                    self.Bwy[i] = np.exp(-tmp**2)
                    Power_sum_square += self.Bwy[i]
                Power_sum = np.sqrt(Power_sum_square)


                
                for i in range(self.nw):
                    self.Bwy[i] = np.sqrt(self.Bwy[i]) * Bw / Power_sum
            elif distribution == "Constant":
                print("Using constant distribution of wave")
                self.Bwy = Bw / np.sqrt(self.nw)
            else:
                raise ValueError("Unknown distribution")


    def generate_parallel_wave(self):
        # get wave phase, wave number, and the amplitude of E (calculated from B)
        self.wce = gyrofrequency(cst.Charge,cst.Me, self.B0)
        self.wpe = plasmafrequency(cst.Charge,cst.Me, self.ne)
        RR = 1 + self.wpe**2 / ((self.wce - self.ws) * self.ws) 
       
        self.k = self.ws * np.sqrt(RR)/cst.C 
        mu = cst.C * self.k / self.ws
        ExByp = self.direction / mu # Ex / By

        self.Ewx = self.Bwy * ExByp

        self.phi0= np.random.rand(self.nw) * 2 * np.pi



        # !! remember to get the random wave phase
        print("The information of wave:")
        print("gyrofrequency:", self.wce)
        print("plasmafrequency", self.wpe) 
        print("wave number k", self.k)
        print("wave frequency", self.ws)
        print("wave amplitude Bwy", self.Bwy)
        print("wave amplitude Ewx", self.Ewx)
        print("Initial phase", self.phi0)



            




        


    