import math
import numpy as np
import matplotlib.pyplot as plt


class GaussianSecondLensSingleValue:
    def __init__(self, w0, lambdaL, defocusing_range, denting_depth):

        print('xxxxxxxx run with denting_depth: ' + str(denting_depth) + 'xxxxxxx')
        self.w0 = w0
        self.w0_fundamental = w0
        self.harmonic_number = 1
        self.lambdaL = lambdaL  # fundamental!!!
        self.z = defocusing_range
        # self.v = 0
        # self.w_new = 0
        self.zr = self.rayleigh_length()
        self.denting_depth = denting_depth
        self.wz = self.beam_waist_of_z_value(self.z)
        #print(self.wz)
        self.f = self.choose_f_dependency(True)
        #print(self.f, 'focal length for z:', self.z)
        self.q = self.q_initial()
        self.harmonic_number_array = np.arange(1, 32, 1)
        self.test_beam_waist()
        # self.beam_waist_of_z_harmonic()

    def q_initial(self):
        self.q = (self.w0 ** 2) * math.pi / (self.lambdaL / self.harmonic_number)
        return self.q

    def rayleigh_length(self):
        # remains Ry as fundamental
        ry = math.pi * (self.w0_fundamental ** 2) / (self.lambdaL)
        print('Rayleigh length: ', ry, 'for wo', self.w0_fundamental)
        return ry

    def beam_waist_of_z_value(self, z):
        self.wz = self.w0_fundamental* ((1 + (z / self.zr) ** 2) ** 0.5)
        #print(self.wz, z, 'wz and z fundamental')
        return self.wz

    def test_beam_waist(self):
        test = np.linspace(-5,5,1000)
        wz = np.zeros([len(test)])

        wz[::] = self.beam_waist_of_z_value(test[::])

        plt.figure(1)
        plt.plot(test,wz)


    def beam_waist_of_z_harmonic(self):
        # definition of rayleigh length - which would scale with 1/harmonic_number
        # can be used is not used now
        zr = math.pi * self.w0 ** 2 / (self.lambdaL / self.harmonic_number)
        wz_harmonic = self.w0 * (1 + (self.z / zr) ** 2) ** 0.5
        return wz_harmonic

    def radius_of_lens(self):
        return (4 * (self.denting_depth ** 2) + self.w0_fundamental ** 2) / (8 * self.denting_depth)

    def radius_of_lens_and_beamwaist(self):
        denting_z = self.denting_depth * (self.w0_fundamental / self.wz)
        radius = (4 * (denting_z ** 2) + (self.wz ** 2)) / (8 * denting_z)
        return radius

    def focal_lens_constant(self):
        return self.radius_of_lens() / 2

    def focal_lens_of_wz(self):
        # print( self.radius_of_lens_and_beamwaist(), 'radius for z:', self.z, self.wz, 'wz')
        self.f = self.radius_of_lens_and_beamwaist() / 2
        name1 = 'focalL(w(z)), Dmax:' + str(self.denting_depth)
        return self.f

    def choose_f_dependency(self, switch):
        if switch == True:
            self.f = self.focal_lens_of_wz()
            print(self.f, 'wz dependent focal lens')
        else:
            self.f = self.focal_lens_constant()
            print(self.f, 'constant focal length')

        return self.f

    def new_focal_position_constant_focal_length(self):
        # this often named v
        self.f = self.focal_lens_constant()
        AA = (self.q ** 2 / self.f) - self.z * (1 - self.z / self.f)
        BB = (self.q ** 2 / self.f ** 2) + (1 - self.z / self.f) ** 2
        return AA / BB

    def new_focal_position_wz(self):
        self.q_initial()
        AA = (self.q ** 2 / self.f) - self.z * (1 - self.z / self.f)
        BB = (self.q ** 2 / self.f ** 2) + (1 - self.z / self.f) ** 2
        # print(self.harmonic_number, 'new v', AA / BB)
        # print('single value for harmonic_number:', self.harmonic_number, 'new focal position', AA/BB)
        return AA / BB

    def new_beam_waist_single_value(self):
        # print('initial beamwaist', self.w0, 'for single value and harmonic_number:', self.harmonic_number)
        v_single = self.new_focal_position_wz()
        new_wz_harmonic = ((1 - v_single / self.f) ** 2) + (1 / self.q ** 2) * (
                self.z + v_single * (1 - (self.z / self.f))) ** 2
        new_wz_harmonic = self.w0 * (new_wz_harmonic ** 0.5)
        # print(new_wz_harmonic, v_single)
        return new_wz_harmonic

    def new_divergence_from_w0_new(self, w_new):

        return self.lambdaL / (self.harmonic_number * math.pi * w_new)

    def switch_sign(self, var):
        return -var

    def resulting_divergence_over_harmonic_number(self):
        result_w0_new = np.zeros([len(self.harmonic_number_array)])

        result_div_harmonic_number = np.zeros([len(self.harmonic_number_array)])

        for x in range(0, len(self.harmonic_number_array)):
            self.harmonic_number = self.harmonic_number_array[x]
            # self.q_initial()
            # print(self.f, 'focal length')
            result_w0_new[x] = self.new_beam_waist_single_value()
            result_div_harmonic_number[x] = self.new_divergence_from_w0_new(result_w0_new[x])

        name1 = 'z: ' + str(self.z) + '[mm]  lens+'

        plt.figure(10)
        plt.plot(self.harmonic_number_array, result_w0_new, label=name1 + 'w(harmonic_number)')
        plt.xlabel = 'harmonic_number'
        plt.ylabel = 'w0(harmonic_number)* in [mm]'
        plt.legend()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.figure(9)
        plt.plot(self.harmonic_number_array, result_div_harmonic_number, label=name1 + 'div detctor', marker='.')
        plt.xlabel = 'harmonic_number'
        plt.ylabel = 'div in [rad]'
        plt.legend()
        # plt.yscale ('log')
        # plt.ylim(0.001,0.006)
        # zip the 2 arrays to get the exact coordinates

    # index = list(zip(index[0])

    def plot_diffraction_limit(self):
        N_list = np.arange(1, 30, 1)
        N_diffraction_limit = np.zeros([29, 1])
        for x in range(0, 30 - 1):
            # halfangle
            N_list[x] = 1 + x
            N_diffraction_limit[x] = (60. / 1500.) / (1 + x)

        plt.figure(9)
        plt.scatter(N_list, N_diffraction_limit, marker="o", color="c", label="Theta_L/harmonic_number")
        plt.hlines(0.007, 0, 30, label="detector limit")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        # plt.savefig("20190123_divergence_mrad_halfangle_and_theo" +".png",  bbox_inches="tight", dpi = 1000)


Test = GaussianSecondLensSingleValue(0.012, 0.0008, 0, 0.00005)
Test.resulting_divergence_over_harmonic_number()
Test = GaussianSecondLensSingleValue(0.012, 0.0008, -1, 0.00005)
Test.resulting_divergence_over_harmonic_number()
Test = GaussianSecondLensSingleValue(0.012, 0.0008, -2, 0.00005)
Test.resulting_divergence_over_harmonic_number()
Test = GaussianSecondLensSingleValue(0.012, 0.0008, -3, 0.00005)
Test.resulting_divergence_over_harmonic_number()
Test.plot_diffraction_limit()



# plt.savefig("caseII_div_over_N_50nm" + ".png", bbox_inches="tight", dpi=1000)

plt.legend()
plt.xlim(1, 28)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

# Test3.resulting_divergence_over_N(2.5)
