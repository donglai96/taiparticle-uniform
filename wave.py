"""
@author Donglai Ma
@email donglaima96@gmail.com
@create date 2023-01-26 20:49:03
@modify date 2023-01-26 20:49:03
@desc taichi dataclass of wave
"""
import taichi as ti
from taichiphysics import * 
import constants as cst


@ti.dataclass
class Wave:
    w:ti.f64
    phi0:ti.f64
    Ewx: ti.f64
    Bwy: ti.f64
    k: ti.f64
    Bw:ti.types.vector(3,ti.f64)
    Ew:ti.types.vector(3,ti.f64)
    phi:ti.f64
    

    @ti.func
    def initialize(self, w, phi0, Ewx, Bwy, k):
        self.w = w
        self.phi0 = phi0
        self.Ewx = Ewx
        self.Bwy = Bwy
        self.k = k
        
    
    @ti.func
    def get_wavefield(self, r, t):
        #print('tttt',t)
        z = r[2]
        
        # get wave wphase
        phi = self.k * z - self.w * t + self.phi0
        By = self.Bwy
        Bx = By

        Ex = self.Ewx
        Ey = Ex

        self.Bw = [Bx * ti.cos(phi),
                   -By * ti.sin(phi),
                   0.0]
        self.Ew = [-Ex * ti.sin(phi),
                    -Ey * ti.cos(phi),
                    0.0]
        self.phi = phi
