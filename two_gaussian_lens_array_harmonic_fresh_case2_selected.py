# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:17:34 2020

@author: similarities
"""
# Gaussian beam profile (focused by lens 1) which passes lens2 close to 
# focal position. The second lens is created by high intense laser target
# interaction, as the light pressure induces a target denting. The dented
# target corresponds to a curved mirror, here treated as lens2
# The resulting beamwaist and beam divergence is calculated in dependency
# of the second lens position relative to the focal position of lens 1.
# The implementation considers different cases with the following initial assumptions:
# case 1: High harmonic radiation follow laser divergence divided by the harmonic number
# in the place of origin (e.g. target out of focus / in focus). This case violates the
# the conservation of wave front curvature in the place of origin.
# case 2: High harmonic radiation
# with either a conversation of the fundamental beam waist w0 of the fundamental
# when the target is not dented but displaced from focal position (like a reflecting
# mirror without curvature). This case corresponds to the conversation of wave front
# curvature for the high harmonic radiation in the place of origin.

# The plasma lens (lens 2) can be chosen to be:
# constant for the whole defocusing range
# to scale with the aperture (laser spot size) on the target ~ 1/w(z)^2
# to scale with the aperture (laser spot size) on the target + a decrease
# of the denting that follows additionally with  ~1/w(z) - most realistic assumption.
# As a consequence, the focal lens of the plasma lens is varying in the last to case over
# the defocusing range very rapidly.

import math
import numpy as np
import matplotlib.pyplot as plt


class GaussianBeamSecondLens:
    def __init__(self, w0, lambdaL, defocusing_range, denting_depth, harmonic_number, case):

        self.w0_fundamental = w0
        self.w0 = w0
        self.harmonic_number = harmonic_number
        self.lambdaL = lambdaL  # fundamental!!!
        self.z = np.arange(-defocusing_range, defocusing_range, 0.01)
        self.wz_harmonic = np.zeros(len(self.z))
        self.v_array = np.zeros(len(self.z))
        self.w_new = np.zeros(len(self.z))
        # print(self.z)
        self.denting_depth = denting_depth
        self.f_array = np.zeros(len(self.z))
        self.wz_fundamental = np.zeros(len(self.z))
        self.denting_array = np.zeros(len(self.z))
        self.zr_fundamental = self.rayleigh_length()
        self.radius_array = np.zeros(len(self.z))
        self.beam_divergence = np.zeros(len(self.z))
        self.q = self.q_initial()
        self.description = str
        self.case = case
        self.select_case_inital_beamwaist()


    def select_case_inital_beamwaist(self):
        if self.case == 2:
            self.w0 = self.w0_fundamental
            self.description = 'case_II'

        elif self.case == 1:
            self.w0 = self.w0_fundamental / self.harmonic_number
            self.description = 'case_I'

        self.q_initial()
        self.beam_waist_of_z_harmonic()

        return self.description

    def q_initial(self):
        self.q = (self.w0 ** 2) * math.pi / (self.lambdaL / self.harmonic_number)
        return self.q

    def rayleigh_length(self):
        # remains Ry as fundamental
        ry = math.pi * (self.w0_fundamental ** 2) / (self.lambdaL)
        print('Rayleigh length: ', ry, 'for wo(fundamental)', self.w0_fundamental)
        return ry

    def diffraction_limit(self):
        focal_length_initial = 1500
        beam_diameter_initial = 60
        self.w0 = focal_length_initial * self.lambdaL / (self.harmonic_number * math.pi * beam_diameter_initial)
        return self.w0

    def beamwaist_of_z(self):
        self.wz_fundamental[::] = self.w0_fundamental * (1 + (self.z[::] / self.zr_fundamental) ** 2) ** 0.5
        print('initial beamwaist with 800nm:', self.w0_fundamental)
        return self.wz_fundamental

    def beam_waist_of_z_harmonic(self):
        # definition of rayleigh length - which would scale with 1/harmonic_number
        # can be used is not used now
        zr = math.pi * self.w0 ** 2 / (self.lambdaL / self.harmonic_number)
        self.wz_harmonic[::] = self.w0 * (1 + (self.z[::] / zr) ** 2) ** 0.5
        # print('harmonic_number:', self.harmonic_number, 'Ry', zr_fundamental, 'w0 inital', self.w0, 'for calculation w(z) before lens2')
        name = 'w(N,z) ' + f'{self.harmonic_number} ' + self.description
        self.plot_results(self.z, self.wz_harmonic, name, 'z mm', 'w(N,z) initial mm', 1, None)
        return self.wz_harmonic

    def radius_of_lens_and_beamwaist_and_intensity(self):
        self.denting_array[::] = self.denting_depth * ((self.w0_fundamental) / (self.wz_fundamental[::]))
        self.radius_array[::] = (4 * (self.denting_array[::] ** 2) + (self.wz_fundamental[::] ** 2)) / (
                8 * self.denting_array[::])
        return self.radius_array

    def radius_of_lens_and_beamwaist(self):
        print(self.denting_depth)
        self.radius_array[::] = (4 * pow(self.denting_depth, 2) + pow(self.wz_fundamental[::], 2)) \
                                / (8 * self.denting_depth)
        return self.radius_array

    def radius_constant(self):
        radius_constant = (4 * (self.denting_depth ** 2) + (self.w0_fundamental ** 2)) / (
                8 * self.denting_depth)
        return radius_constant

    def choose_focal_length_dependency(self, switch):
        if switch == 'w0_aperture_and_IL_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist_and_intensity() * 0.5
            print(np.min(self.f_array), 'minimum focal length apterure and intensity dependent')
            self.description = self.description + ' ~w(Z)^3'

        elif switch == 'w0_aperture_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist() * 0.5
            print(np.min(self.f_array), 'minimum focal length apterture w(z) dependent')
            self.description = self.description + ' ~w(z)^2'

        else:
            self.f_array[::] = self.radius_constant() * 0.5
            print('constant focal length', self.f_array[10])
            self.description = self.description + ' f constant'

        self.plot_results(self.z, self.f_array, self.description, 'z mm', 'focal length mm', 2, None)
        return self.f_array, self.description

    def new_focal_position_single_value(self, index):
        i = index
        AA = (self.q ** 2 / self.f_array[i]) - self.z[i] * (1 - self.z[i] / self.f_array[i])
        BB = (self.q ** 2 / self.f_array[i] ** 2) + (1 - self.z[i] / self.f_array[i]) ** 2
        return AA / BB

    def new_beam_waist_single_value(self, index):
        v_single = self.new_focal_position_single_value(index)
        new = ((1 - v_single / self.f_array[index]) ** 2) + (1 / self.q ** 2) * (
                self.z[index] + v_single * (1 - (self.z[index] / self.f_array[index]))) ** 2
        new = self.w0 * (new ** 0.5)
        return new

    def new_divergence_from_w0_new(self, w_new):
        return self.lambdaL / (self.harmonic_number * math.pi * w_new)

    def create_v_array(self):
        for i in range(0, len(self.z)):
            self.v_array[i] = self.new_focal_position_single_value(i)
        name1 = 'v(w(z)) for Dmax: ' + str(self.denting_depth) + 'N: ' + str(self.harmonic_number) + self.description
        self.plot_results(self.z, self.v_array, name1, 'z mm', 'new focal position mm', 5, None)
        # plt.savefig("caseII_f_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        return self.v_array

    def new_beamwaist(self):
        #print('for harmonic_number: ', self.harmonic_number, 'we start with a beamwaist for lens 1 of about:', self.w0)
        self.create_v_array()
        # print(self.f, 'insert in wnew')
        self.w_new[::] = ((1 - self.v_array[::] / self.f_array[::]) ** 2) + (1 / self.q ** 2) * (
                self.z[::] + self.v_array[::] * (1 - (self.z[::] / self.f_array[::]))) ** 2
        self.w_new[::] = self.w0 * (self.w_new[::] ** 0.5)
        name = 'w0_new for N:' + str(self.harmonic_number) + self.description
        self.plot_results(self.z, self.w_new, name, 'z mm', 'w(z)', 4, None)
        return self.w_new

    def resulting_beam_divergence(self):
        self.new_beamwaist()
        self.beam_divergence[::] = (self.lambdaL / (self.harmonic_number * math.pi * self.w_new[::]))
        name = 'Theta Dmax:' + str(self.denting_depth) + 'N: ' + str(self.harmonic_number) + self.description
        self.plot_results(self.z, self.beam_divergence, name, 'z mm', 'w(N,z)', 6, None)
        # plt.savefig("caseII_Theta_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)

    def resulting_divergence_over_N(self, z):
        # intensity (position) dependent focal length lens 2
        index = list(zip(*np.where(self.z >= z)))
        index = index[0]
        # print(index, self.z[index])
        harmonic_number_array = np.arange(1, 32, 1)
        result_w0_N = np.zeros(len(harmonic_number_array))
        result_div_N = np.zeros(len(harmonic_number_array))
        for x in range(0, len(harmonic_number_array)):
            self.harmonic_number = harmonic_number_array[x]
            self.q_initial()
            result_w0_N[x] = self.new_beam_waist_single_value(index)
            result_div_N[x] = self.new_divergence_from_w0_new(result_w0_N[x])
        name1 = self.description + 'z: ' + str(z) + '[mm] '
        self.plot_results(harmonic_number_array, result_w0_N, name1, 'N', 'w0(N) mm', 10, marker='.')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        self.plot_results(harmonic_number_array, result_div_N, name1, 'N', 'Theta(N) mm', 9, marker='.')
        # zip the 2 arrays to get the exact coordinates

    # index = list(zip(index[0])
    def plot_diffraction_limit(self):
        N_list = np.arange(1, 30, 1)
        N_diffraction_limit = np.zeros([29, 1])
        for x in range(0, 30 - 1):
            # halfangle
            N_list[x] = 1 + x
            N_diffraction_limit[x] = (60. / 1500.) / (1 + x)

        self.plot_results(N_list, N_diffraction_limit, 'Theta(L)/N', 'N', 'Theta rad', 9, 'o')

        # plt.savefig("20190123_divergence_mrad_halfangle_and_theo" +".png",  bbox_inches="tight", dpi = 1000)

    def plot_case_vincenti(self):
        N_list = np.arange(1, 30, 1)
        N_div_vincenti = np.zeros([29, 1])
        denting_part = (4 * math.pi * self.denting_depth) ** 2
        for x in range(0, 30 - 1):
            # halfangle
            N_list[x] = 1 + x
            N_div_vincenti[x] = (1 / (math.pi * self.w0_fundamental)) * (
                    denting_part + (self.lambdaL / (1 + x)) ** 2) ** 0.5
        self.plot_results(N_list, N_div_vincenti, 'Theta(L)/N vincenti model', 'N', 'Theta rad', 9, 'o')

    def plot_results(self, x, y, label, name_x, name_y, figure_number, marker):
        plt.figure(figure_number)
        plt.plot(x, y, label=label, marker=marker)
        plt.xlabel = name_x
        plt.ylabel = name_y
        plt.legend()
        plt.get_figlabels()

#GaussianBeamSecondLens(focal_radius[1/e in mm], wavelength_fundamental [mm], defocusing_range [mm], denting in [mm], harmonic_number [int], case selection [1 or 2])
# case_selection: select:1 for initial beamwaist scales with 1/harmonic number, select 2: for beamwaist of harmonic is just beamwaist fundamental (scales sourcesize 1/harmonic_number)
# choose the dependency of the focal length function : 'w0_aperture_and_IL_dependent', 'w0_aperture_dependent' or False
Test = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 12, 2)
Test.choose_focal_length_dependency('w0_aperture_and_IL_dependent')
Test.resulting_beam_divergence()
Test.resulting_divergence_over_N(-2.)
Test.resulting_divergence_over_N(0.)
Test.resulting_divergence_over_N(-1.)
Test.resulting_divergence_over_N(1.)
Test.plot_diffraction_limit()

Test2 = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 15, 2)
Test2.choose_focal_length_dependency('w0_aperture_and_IL_dependent')
Test2.resulting_beam_divergence()
Test2 = GaussianBeamSecondLens(0.012, 0.0008, 5, 0.00005, 32, 2)
Test2.choose_focal_length_dependency('w0_aperture_and_IL_dependent')
Test2.resulting_beam_divergence()
Test2.plot_case_vincenti()
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# plt.savefig("caseII_div_over_N_50nm" + ".png", bbox_inches="tight", dpi=1000)
plt.show()
