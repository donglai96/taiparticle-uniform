"""
@author Donglai Ma
@email donglaima96@gmail.com
@create date 2023-01-25 16:35:44
@modify date 2023-01-25 16:35:44
@desc test particle simulation with taichi code
This file is for particle class
each particle is a struct in taichi
"""

import taichi as ti
import constants as cst
# from uniformfield import uniform_field
# from taichiphysics import *
# from wavepacket import Wave

@ti.dataclass
class Particle:
    r: ti.types.vector(3,ti.f64)
    p: ti.types.vector(3,ti.f64) #momentum
    #gamma: ti.f64
    m: ti.f64
    q: ti.f64
    t: ti.f64
    #deltat: ti.f64
    alpha: ti.f64 #pitch angle
    alpha0: ti.f64 #equator pitch angle
    phi: ti.f64 # record particle phase
    Ep: ti.types.vector(3,ti.f64) #record E
    Bp: ti.types.vector(3,ti.f64) #record B
    @ti.func
    def initParticles(self, mm, qq):
        self.m = mm
        self.q = qq
        self.t = 0
        self.phi = 0
        self.Ep = ti.Vector([0, 0, 0])
        self.Bp = ti.Vector([0,0,0])
        #self.deltat = deltat
    @ti.func
    def initPos(self, x, y, z):
        self.r = ti.Vector([x, y, z]) #z would be latitude
    @ti.func
    def initMomentum(self, px, py, pz):
        self.p = ti.Vector([px, py, pz])
    
    # @ti.func
    # def getgamma(self):
    #     self.gamma = ti.sqrt(1 + self.p.norm()**2 /(self.m**2 * cst.C**2))
    #     return self.gamma
    @ti.func
    def get_pitchangle(self):
        self.alpha = ti.acos(self.p[2]/self.p.norm())
    

    @ti.func
    def boris_push(self, dt, E, B):
        """Push the particles using Boris' method
        Update the velocity of particles
        An example for non-relativistic:
         https://www.particleincell.com/2011/vxb-rotation/
        p_n_minus = p_minus - qE*dt/2
        p_n_plus = p_plus + qE*dt/2
        (p_plus - p_minus) / dt = q/2 * (p_plus + p_minus)/(gamma * m0)
        ...
        Args:
            dt (_type_): _description_
            E (_type_): _description_
            B (_type_): _description_
        """

        p_minus = self.p + self.q * E * dt / 2.0
 
        gamma = ti.sqrt(1 + p_minus.norm()**2 /(self.m**2 * cst.C**2))
        #print('gamma',gamma)
        #print(self.q)
        t = self.q * dt * B / (2 * gamma* cst.C * self.m) # Gaussian unit
        
        p_p = p_minus + p_minus.cross(t)

        s = 2 * t /(1 + t.dot(t))

        p_plus = p_minus + p_p.cross(s)

        self.p = p_plus + self.q * E * dt / 2.0
        #print(self.p)
        

    @ti.func
    def leap_frog(self, dt, E, B):

        # The idead of lp is first half step, then boris and the leap frog
        # or 
        self.boris_push(dt, E, B)
        gammam = self.m * ti.sqrt(1 + self.p.norm()**2 /(self.m**2 * cst.C**2))
        
        v_xyz = self.p / gammam
        #print('xyz',v_xyz)
        

    
       
        self.r += dt * v_xyz 