# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 12:27:09 2020

@author: similarities
"""

import numpy as np
import matplotlib.pyplot as plt


class HoleBoringVelocity:

    def __init__(self, ne1, ne2, Z, A, a0):
        self.me_mp_ratio = 1 / 1836
        self.Z = Z
        self.ne1 = ne1
        self.ne2 = ne2
        self.A = A
        self.nc_ne = np.linspace(self.ne1, self.ne2, 100)
        self.vHB = np.zeros([len(self.nc_ne), 1])
        self.beta = float
        self.pulse_duration = 20 #in fs
        # print(self.nc_ne)
        # print(self.vHB)
        self.a0 = a0

    def calc_HB(self):
        for x in range(0, len(self.nc_ne)):
            self.vHB[x] = (self.Z / self.A) * (self.me_mp_ratio) * (1 / self.nc_ne[x])
            self.vHB[x] = ((self.vHB[x]) ** 0.5) * self.a0

        return self.vHB

    def intensity_to_a0(self, intensity):
        self.a0 = ((intensity * 1E19 / 1.38 * 1E-18) * (0.8 ** 2)) ** 0.5
        return self.a0
    
    def a0_to_intensity(self):
        return (self.a0 ** 2) * 1.38 * 1E-18/ ((0.8 **2)*1E19)
    
    
    def hole_boring_velocity_single_value(self):
        self.beta = (self.Z / self.A) * (self.me_mp_ratio) * (1 / self.nc_ne[x])
        self.beta = ((self.beta) ** 0.5) * self.a0
        return self.beta
    
    def change_intensity_by_beamwaist(self, peak_intensity, w0, wz):
        return peak_intensity*(w0 **2)/ (wz **2)
    
    def change_intensity_by_pulse_duration(self, peak_intensity, pulse_duration_old, pulse_duration_new):
        return peak_intensity *(pulse_duration_old/pulse_duration_new)
    
    def change_pulse_duration(self, new_pulse_duration):
        self.pulse_duration = new_pulse_duration
        return new_pulse_duration
    
    def denting_maximum_half_pulse_duration(self):
        denting_max = 1E9 * self.vHB * 3E8 * self.pulse_duration * 0.5 * 1E-15
        return denting_max

    def plot_denting(self):
        plt.figure(2)
        plt.plot(self.nc_ne, self.denting_maximum_half_pulse_duration(), label='Dmax for a0: ' + str(self.a0))
        plt.xlabel('ne/nc')
        plt.ylabel('Dmax nm')
        plt.yscale('log')
        plt.xlim(80, 200)
        plt.legend()

    def plot_HB(self):
        plt.figure(1)
        plt.plot(self.nc_ne, self.vHB, label='a0 = ' + str(self.a0))
        plt.xlabel("ne/nc")
        plt.ylabel("beta HB")
        plt.legend()

#def __init__(self, ne1, ne2, Z, A, a0):
a0_1 = HoleBoringVelocity(10, 200, 8, 10, 0.1)
a0_1.calc_HB()
a0_1.plot_HB()
a0_1.plot_denting()

a0_1 = HoleBoringVelocity(10, 200, 8, 10, 3)
a0_1.calc_HB()
a0_1.plot_HB()
a0_1.plot_denting()

a0_1 = HoleBoringVelocity(10, 200, 8, 10, 1)
a0_1.calc_HB()
a0_1.plot_HB()
a0_1.plot_denting()
cd
a0_1 = HoleBoringVelocity(10, 200, 8, 10, 6)
a0_1.calc_HB()
a0_1.plot_HB()
a0_1.plot_denting()


def a0():
    intensity = np.linspace(1, 10000, 1000)
    intensity[::] = intensity[::] * 1E-3
    a0 = ((intensity * 1E19 / 1.38 * 1E-18) * ((0.8 )** 2)) ** 0.5
    plt.figure(3)
    plt.plot(intensity, a0)
    plt.xlabel('intensity in 1E19')
    plt.ylabel('a0')
    plt.xscale('log')
    plt.legend()


a0()






plt.show()
