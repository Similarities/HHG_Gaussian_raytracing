# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:17:34 2020

@author: similarities
"""


# The plasma lens (lens 2) can be chosen to be:
# constant for the whole defocusing range
# to scale with the aperture (laser spot size) on the target ~ 1/w(z)^2
# to scale with the aperture (laser spot size) on the target + a decrease
# of the denting that follows additionally with  ~1/w(z) - most realistic assumption.
# As a consequence, the focal lens of the plasma lens is varying in the last to case over
# the defocusing range very rapidly.

#corrected jan 2022 for geometrical factor of two in beamwaist

import math
import numpy as np
import matplotlib.pyplot as plt


class GaussianBeamSecondLens:
    def __init__(self, w0, lambdaL, defocusing_range, denting_depth, color):

        self.w0_fundamental = w0
        self.w0 = w0
        self.lambdaL = lambdaL  # fundamental!!!
        self.z = np.arange(-defocusing_range, defocusing_range, 0.01)
        # print(self.z)
        self.denting_depth = denting_depth
        self.f_array = np.zeros(len(self.z))
        self.wz_fundamental = np.zeros(len(self.z))
        self.denting_array = np.zeros(len(self.z))
        self.zr_fundamental = self.rayleigh_length()
        self.radius_array = np.zeros(len(self.z))
        self.beam_divergence = np.zeros(len(self.z))
        self.description = "focal_length"
        self.color = color



    def rayleigh_length(self):
        # remains Ry as fundamental
        #M = Msquare, is M^2 already quadratic?
        M = 2.3
        ry = math.pi * (self.w0_fundamental ** 2) * M / (self.lambdaL)
        print('Rayleigh length: ', ry, 'for wo(fundamental)', self.w0_fundamental)
        return ry


    def beamwaist_of_z(self):
        self.wz_fundamental[::] = self.w0_fundamental * (1 + (self.z[::] / self.zr_fundamental) ** 2) ** 0.5
        print('initial beamwaist with 800nm:', self.w0_fundamental)
        return self.wz_fundamental


    def radius_of_lens_and_beamwaist_and_intensity(self):
        #assumes denting ~ a0 ~ (IL ** 0.5) ~ (w0**2/wz**2)**0.5
        self.denting_array[::] = self.denting_depth * ((self.w0_fundamental) / (self.wz_fundamental[::]))
        self.radius_array[::] = (4 * (self.denting_array[::] ** 2) + ((self.wz_fundamental[::]*2) ** 2)) / (
                8 * self.denting_array[::])
        return self.radius_array

    def radius_of_lens_and_beamwaist_and_intensity_and_chirp(self, factor):
        #assumes denting ~ a0 ~ (IL ** 0.5) ~ (w0**2/wz**2)**0.5 * (tauL/tChirped)**0.5
        # factor = tauL/tauChirped
        self.denting_array[::] = self.denting_depth * ((self.w0_fundamental) / (self.wz_fundamental[::]))*(factor ** 0.5)
        self.radius_array[::] = (4 * (self.denting_array[::] ** 2) + ((self.wz_fundamental[::]*2) ** 2)) / (
                8 * self.denting_array[::])
        return self.radius_array

    def radius_of_lens_and_beamwaist(self):
        print(self.denting_depth)
        self.radius_array[::] = (4 * pow(self.denting_depth, 2) + pow((self.wz_fundamental[::]*2), 2)) \
                                / (8 * self.denting_depth)
        return self.radius_array

    def radius_constant(self):
        radius_constant = (4 * (self.denting_depth ** 2) + ((self.wz_fundamental[::]*2) ** 2)) / (
                8 * self.denting_depth)
        return radius_constant

    def choose_focal_length_dependency(self, switch):
        if switch == 'w0_aperture_and_IL_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist_and_intensity() * 0.5
            print(np.min(self.f_array), 'minimum focal length apterure and intensity dependent')
            self.description = ' D(w(z), IL(wz)) Dmax: ' + str(self.denting_depth*1E6) + "nm"
            marker = "dashed"


        elif switch == 'w0_aperture_dependent':
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist() * 0.5
            print(np.min(self.f_array), 'minimum focal length apterture w(z) dependent')
            self.description =  ' D(w(z)) Dmax: ' + str(self.denting_depth *1E6) + "nm"
            marker = "dotted"


        elif switch == "w0_apertur_and_IL_and_chirp_dependent":
            self.beamwaist_of_z()
            self.f_array[::] = self.radius_of_lens_and_beamwaist_and_intensity_and_chirp(25/120) * 0.5
            print(np.min(self.f_array), 'minimum focal length apterture w(z), IL(w(z)) and tauL/tauChirp dependent')
            self.description = ' D(wz, I(wz), tauChirp)  Dmax: ' +str(self.denting_depth * 1E6) + "nm"
            marker = "solid"



        else:
            self.f_array[::] = self.radius_constant() * 0.5
            print('constant focal length', self.f_array[10])
            self.description = ' f constant Dmax: ' + str(self.denting_depth *1E6) +"nm"
            marker = None


        self.plot_results(self.z, self.f_array, self.description, 'z mm', 'focal length mm', 2, marker)
        return self.f_array, self.description



    def plot_results(self, x, y, label, name_x, name_y, figure_number, marker):
        plt.figure(figure_number)
        plt.plot(x, y, label=label, linestyle= marker, color = self.color)
        plt.xlabel(name_x)
        plt.ylabel(name_y)
        plt.legend()
        #plt.get_figlabels()


#GaussianBeamSecondLens(focal_radius[1/e in mm], wavelength_fundamental [mm], defocusing_range [mm], denting in [mm], harmonic_number [int], case selection [1 or 2])
# case_selection: select:1 for ini7ial beamwaist scales with 1/harmonic number, select 2: for beamwaist of harmonic is just beamwaist fundamental (scales sourcesize 1/harmonic_number)
# choose the dependency of the focal length function : 'w0_aperture_and_IL_dependent', 'w0_aperture_dependent' or False

#w0 [mm], lambda[mm], a0, DentingDepth[mm].

Test = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.0005, "c")
Test.choose_focal_length_dependency('w0_apertur_and_IL_and_chirp_dependent')


Test2 = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.0005, "c")
Test2.choose_focal_length_dependency('w0_aperture_and_IL_dependent')

Test2 = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.0005, "c")
Test2.choose_focal_length_dependency('w0_aperture_dependent')

Test = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.00050, "b")
Test.choose_focal_length_dependency('w0_apertur_and_IL_and_chirp_dependent')


Test2 = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.0010, "b")
Test2.choose_focal_length_dependency('w0_aperture_and_IL_dependent')

Test2 = GaussianBeamSecondLens(0.007, 0.0008, 10, 0.00100, "b")
Test2.choose_focal_length_dependency('w0_aperture_dependent')

plt.yscale("log")
plt.legend()
plt.xlim(-4,4)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig("plasma_lens_units_a0_10_D0.5um_1um" + ".png", bbox_inches="tight", dpi=1000)
plt.show()
