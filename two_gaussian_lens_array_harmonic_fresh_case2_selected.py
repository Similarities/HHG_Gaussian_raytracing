# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:17:34 2020

@author: similarities
"""
# Gaussian beam profile (focused by lens 1) which passes lens2 close to 
# focal position. The resulting beamwaist and new focal position 
# here only for inward denting and placed (last plots) before focus

import math
import numpy as np
import matplotlib.pyplot as plt


class GaussianBeamSecondLens:
    def __init__(self, w0, lambdaL, defocusing_range, denting_depth, harmonic_number):
        self.w0 = w0
        self.w0_fundamental = w0
        self.harmonic_number = harmonic_number
        self.lambdaL = lambdaL  # fundamental!!!
        self.z = np.arange(-defocusing_range, defocusing_range, 0.01)
        self.v_array = np.zeros(len(self.z))
        self.w_new = np.zeros(len(self.z))
        # print(self.z)
        self.denting_depth = denting_depth
        self.f_array = np.zeros(len(self.z))
        self.wz_array = np.zeros(len(self.z))
        self.wz_array_harmonic_number = np.zeros(len(self.z))
        self.denting_array = np.zeros(len(self.z))
        self.zr = self.rayleigh_length()
        print(self.zr)
        self.radius_array = np.zeros(len(self.z))
        self.beamdivergence_array = np.zeros(len(self.z))
        self.wz = float
        self.f = float
        self.q = self.q_initial()
        self.N_array = np.arange(1, 32, 1)
        self.focal_lens_of_wz()
        self.beam_waist_of_z_harmonic()

    def q_initial(self):
        self.q = (self.w0_fundamental ** 2) * math.pi / (self.lambdaL / self.harmonic_number)
        # print(self.harmonic_number, self.q)
        return self.q

    def rayleigh_length(self):
        # remains ry as fundamental
        ry = math.pi * (self.w0_fundamental ** 2) / (self.lambdaL)
        print('Rayleigh length: ', ry, 'for wo', self.w0_fundamental)
        return ry

    def beamwaist_of_z(self):
        self.wz_array[::] = self.w0_fundamental * (1 + (self.z[::] / self.zr) ** 2) ** 0.5
        return self.wz_array

    def beam_waist_of_z_harmonic(self):
        # definition of rayleigh length - which would scale with 1/harmonic_number
        # can be used is not used now
        zr = math.pi * self.w0 ** 2 / (self.lambdaL / self.harmonic_number)
        # print('calculation of w(z) for harmonic_number:', self.harmonic_number, 'with new Ry(harmonic_number)', zr)
        self.wz_array_harmonic_number[::] = self.w0 * (1 + (self.z[::] / zr) ** 2) ** 0.5
        # print('harmonic_number:', self.harmonic_number, 'Ry', zr, 'w0 inital', self.w0, 'for calculation w(z) before lens2')
        name = 'harmonic_number: ' + f'{self.harmonic_number}'
        plt.figure(2)
        plt.plot(self.z, self.wz_array_harmonic_number, label=name)
        plt.legend()
        return self.wz_array_harmonic_number

    def radius_of_lens_and_beamwaist_and_intensity(self):
        self.denting_array[::] = self.denting_depth * ((self.w0_fundamental) / (self.wz_array[::]))
        self.radius_array[::] = (4 * (self.denting_array[::] ** 2) + (self.wz_array[::] ** 2)) / (
                8 * self.denting_array[::])
        return self.radius_array

    def radius_of_lens_and_beamwaist(self):

        self.radius_array[::] = (4 * (self.denting_depth ** 2) + (self.wz_array[::] ** 2)) / (
                8 * self.denting_depth)
        return self.radius_array

    def radius_constant(self):
        radius_constant = (4 * (self.denting_depth ** 2) + (self.w0_fundamental ** 2)) / (
                8 * self.denting_depth)
        return radius_constant / 2

    def choose_focal_length_dependency(self, switch):
        if switch == 'w0_aperture_and_IL_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist_and_intensity() * 0.5


        elif switch == 'w0_aperture_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist() * 0.5


        else:
            self.f_array[::] = self.radius_constant() * 0.5

        return self.f_array

    def focal_lens_of_wz(self):

        name1 = 'focalL(w(z)), D0:' + str(self.denting_depth)
        plt.figure(1)
        plt.plot(self.z, self.f_array, label=name1)
        plt.xlabel = 'defocusing [mm]'
        plt.ylabel = 'focal length [mm]'
        plt.legend()
        return self.f_array

    def new_focal_position_single_value(self, index):
        i = index
        self.q_initial()
        AA = (self.q ** 2 / self.f_array[i]) - self.z[i] * (1 - self.z[i] / self.f_array[i])
        BB = (self.q ** 2 / self.f_array[i] ** 2) + (1 - self.z[i] / self.f_array[i]) ** 2
        # print('single value for harmonic_number:', self.harmonic_number, 'new focal position', AA/BB)
        return AA / BB

    def new_beam_waist_single_value(self, index):
        # print

        v_single = self.new_focal_position_single_value(index)
        new = ((1 - v_single / self.f_array[index]) ** 2) + (1 / self.q ** 2) \
              * (self.z[index] + v_single * (1 - (self.z[index] / self.f_array[index]))) ** 2
        new = self.w0 * (new ** 0.5)
        return new

    def new_divergence_from_w0_new(self, w_new):
        return self.lambdaL / (self.harmonic_number * math.pi * w_new)

    def create_v_array(self):

        for i in range(0, len(self.z)):
            self.v_array[i] = self.new_focal_position_single_value(i)
        name1 = 'f(w(z)) for Dmax: ' + str(self.denting_depth) + 'harmonic_number: ' + str(self.harmonic_number)
        name2 = 'f(w(z)) for Dmax: ' + str(self.denting_depth)
        plt.figure(4)
        plt.plot(self.z, self.v_array, label=name1)
        plt.xlabel = 'defocusing [mm]'
        plt.ylabel = 'new focal position [mm]'
        plt.legend()
        # plt.savefig("caseII_zN_over_N_50nm" + ".png", bbox_inches="tight", dpi=1000)
        plt.figure(3)
        plt.plot(self.z, self.f_array, label=name2)
        plt.ylabel = 'focal lens [mm]'
        plt.xlabel = 'defocusing [mm]'

        plt.legend()
        # plt.savefig("caseII_f_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        return self.v_array

    def new_beamwaist(self):

        # self.q_initial()
        # print('for harmonic_number: ', self.harmonic_number, 'we start with a beamwaist for lens 1 of about:', self.w0)
        self.create_v_array()
        self.w_new[::] = ((1 - self.v_array[::] / self.f_array[::]) ** 2) + (1 / self.q ** 2) * (
                self.z[::] + self.v_array[::] * (1 - (self.z[::] / self.f_array[::]))) ** 2
        self.w_new[::] = self.w0 * (self.w_new[::] ** 0.5)

        plt.figure(5)

        name = 'new beamwaist of f(w(z)), D0: ' + str(self.denting_depth) + 'harmonic_number: ' + str(
            self.harmonic_number)
        plt.plot(self.z, self.w_new, label=name)
        plt.xlabel = 'defocusing [mm]'
        plt.ylabel = 'w_new(z) [mm]'
        plt.legend()
        # plt.show()
        # plt.savefig("caseII_wN_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        return self.w_new

    def resulting_beam_divergence(self):
        self.new_beamwaist()
        self.beamdivergence_array[::] = (self.lambdaL / (self.harmonic_number * math.pi * self.w_new[::]))
        name = 'Theta Dmax:' + str(self.denting_depth) + 'harmonic_number: ' + str(self.harmonic_number)
        plt.figure(6)
        plt.plot(self.z, self.beamdivergence_array, label=name)
        plt.yscale('log')
        # plt.ylim(1,6)
        plt.xlabel = 'defocusing [mm]'
        plt.ylabel = '[rad]'
        plt.legend()
        # plt.savefig("caseII_Theta_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)

    def resulting_divergence_over_N(self, z):

        # intensity (position) dependent focal length lens 2

        index = list(zip(*np.where(self.z >= z)))
        index = index[0]
        # print(index, self.z[index])
        result_w0_N = np.zeros(len(self.N_array))
        result_div_N = np.zeros(len(self.N_array))

        for x in range(0, len(self.N_array)):
            self.harmonic_number = self.N_array[x]

            self.q_initial()

            result_w0_N[x] = self.new_beam_waist_single_value(index)
            print(result_w0_N[x], 'w0 new for z', self.z[x], 'N:', self.harmonic_number)

            result_div_N[x] = self.new_divergence_from_w0_new(result_w0_N[x])

        name1 = 'z: ' + str(z) + '[mm]' + 'w(0,N)'
        name2 = 'z: ' + str(z) + '[mm]' + 'div(z,N)'

        plt.figure(7)
        plt.plot(self.N_array, result_w0_N, label=name1)
        plt.xlabel = 'harmonic_number'
        plt.ylabel = 'w0(harmonic_number)* in [mm]'
        plt.legend()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.figure(8)
        plt.plot(self.N_array, result_div_N, label=name2, marker='.')
        plt.xlabel = 'harmonic_number'
        plt.ylabel = 'div in [rad]'
        plt.legend()
        plt.yscale('log')

    def plot_diffraction_limit(self):
        N_list = np.arange(1, 30, 1)
        N_diffraction_limit = np.zeros([29, 1])
        for x in range(0, 30 - 1):
            # halfangle
            N_list[x] = 1 + x
            N_diffraction_limit[x] = (60. / 1500.) / (1 + x)

        plt.figure(8)
        plt.scatter(N_list, N_diffraction_limit, marker="o", color="c", label="Theta_L/harmonic_number")

        plt.hlines(0.007, 0, 30, label="detector limit")
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# choose the dependency of the focal length function : 'w0_aperture_and_IL_dependent', 'w0_aperture_dependent' or False

Test = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 1)
Test.choose_focal_length_dependency('w0_aperture_dependent')
Test.resulting_beam_divergence()
Test.resulting_divergence_over_N(-2.)
Test.resulting_divergence_over_N(0.)
Test.resulting_divergence_over_N(-1.)
Test.resulting_divergence_over_N(1.)
Test.plot_diffraction_limit()

Test2 = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 12)
Test2.choose_focal_length_dependency('w0_aperture_dependent')
Test2.resulting_beam_divergence()
Test2 = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 22)
Test2.choose_focal_length_dependency('w0_aperture_dependent')
Test2.resulting_beam_divergence()
# plt.savefig("caseII_div_over_N_50nm" + ".png", bbox_inches="tight", dpi=1000)


# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()
