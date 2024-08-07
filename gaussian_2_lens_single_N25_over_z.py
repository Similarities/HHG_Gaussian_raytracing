import math
import numpy as np
import matplotlib.pyplot as plt

#git not upstream!
# corrected geometric radius of mirror curvature equation


class GaussianSecondLensSingleValue:
    def __init__(self, w0, lambdaL, denting_depth, color, marker):

        print('xxxxxxxx run with denting_depth: ' + str(denting_depth) + 'xxxxxxx')
        self.w0 = w0
        print("input angle in pi", np.cos(35/180*np.pi))
        print("increase of w0 by angle", 1/(np.cos(35/180*np.pi)))
        self.w0_fundamental = w0
        self.harmonic_number = 25
        self.lambdaL = lambdaL  # fundamental!!!
        #note: self.z is relative to fundamental focal position. z<0 calculates beamwaist (etc.) before focus
        self.z = np.linspace(-4, 4, 500)

        self.focal_lenght_wz_tauL = np.zeros(len(self.z))
        self.focal_lenght_wz = np.zeros(len(self.z))
        self.focal_lenght = np.zeros(len(self.z))

        self.zr = self.rayleigh_length()
        self.denting_depth = denting_depth
        self.wz = self.beam_waist_of_z_value(self.z)

        # focal length can be caluculated via different approaches:
        # " " = constant (determined by denting depth, geometrical and independent of defocusing)
        # "wzIL" denting depth is function of a0~ (IL(wz))**0.5 (compressed pulse)
        # "wzILtauChirp" denting depth is a function of a0 ~ (IL(wz, tauL))**0.5
        # -> pulse duration correction: defined in self.radius_of_lens_and_beamwaist_and_intensity_and_chirp()
        self.f = self.choose_f_dependency("wzILtauChirp")
        #print(self.f, 'focal length for z:', self.z)
        self.q = self.q_initial()
        self.harmonic_number_array = np.arange(1, 32, 1)
        self.test_beam_waist()
        # self.beam_waist_of_z_harmonic()

    def q_initial(self):
        self.q = (self.w0 ** 2) * math.pi / (self.lambdaL / self.harmonic_number)
        return self.q

    def q_initial_vincenti(self):
        self.q = (self.beam_waist_of_z_harmonic_scaling_wN(0) ** 2) * math.pi / (self.lambdaL / self.harmonic_number)
        return self.q


    def rayleigh_length(self):
        # remains Ry as fundamental
        M2 = 2.2 # from experimental values
        ry = math.pi * (self.w0_fundamental ** 2)*M2/ (self.lambdaL)
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
        plt.plot(test,wz, label = "w(z,N) with sourcesize = w0 ")
        plt.legend()
        #plt.show()


    def beam_waist_of_z_harmonic(self):
        # definition of rayleigh length - which would scale with 1/harmonic_number
        # can be used is not used now
        zr = math.pi * self.w0 ** 2 / (self.lambdaL / self.harmonic_number)
        wz_harmonic = self.w0 * (1 + (self.z / zr) ** 2) ** 0.5
        return wz_harmonic

    def test_beam_waist2(self):
        plt.title("beamwaist w(z,N)_N="+str(self.harmonic_number))
        test = np.linspace(-5,5,1000)
        wz = np.zeros([len(test)])
        wz[::] = self.beam_waist_of_z_harmonic_scaling_wN(test[::])
        plt.figure(1)
        plt.plot(test,wz, label = "w(z) with source size dependency ")
        plt.legend()
        #plt.show()

    def beam_waist_of_z_harmonic_scaling_wN(self,z):
        # definition of rayleigh length - which would scale with 1/harmonic_number
        # can be used is not used now
        #former wN = self.w0/self.harmonic_number
        wN= self.w0*(0.72-9*self.harmonic_number*10**-3)
        zr = math.pi * (wN ** 2) / (self.lambdaL / self.harmonic_number)
        wz_harmonic = wN * (1 + (z / zr) ** 2) ** 0.5

        return wz_harmonic

    def radius_of_lens(self):
        return (4 * (self.denting_depth ** 2) + (self.w0_fundamental/np.cos(35/(np.pi*180))*2) ** 2) / (8 * self.denting_depth)

    def radius_of_lens_and_beamwaist(self):
        denting_z = self.denting_depth * (self.w0_fundamental / self.wz)
        radius = (4 * (denting_z ** 2) + ((self.wz/np.cos(35/(np.pi*180))*2) ** 2)) / (8 * denting_z)
        return radius

    def radius_of_lens_and_beamwaist_and_intensity_and_chirp(self):
        #assumes denting ~ a0 ~ (IL ** 0.5) ~ (w0**2/wz**2)**0.5
        # factor = (tauL/tauChirped)
        factor = (25/ 120) ** 0.5
        #print("intensity reduction by chirp:", factor**2)
        denting_z = (self.denting_depth * ((self.w0_fundamental) / ((self.wz))) *factor)
        radius = (4 * (denting_z ** 2) + ((self.wz/np.cos(35/(np.pi*180))*2)** 2)) / (
                8 * denting_z)
        return radius

    def focal_lens_constant(self):
        return self.radius_of_lens() / 2

    def focal_lens_of_wz(self):
        # print( self.radius_of_lens_and_beamwaist(), 'radius for z:', self.z, self.wz, 'wz')
        self.f = self.radius_of_lens_and_beamwaist() / 2
        name1 = 'f (w(z), IL(wz)), Dmax:' + str(self.denting_depth)
        return self.f

    def focal_lens_of_wz_tauL(self):
        # print( self.radius_of_lens_and_beamwaist(), 'radius for z:', self.z, self.wz, 'wz')
        self.f = self.radius_of_lens_and_beamwaist_and_intensity_and_chirp() / 2
        name1 = 'f (w(z), IL(wz), tauL), Dmax:' + str(self.denting_depth)
        return self.f


    def choose_f_dependency(self, switch):
        if switch == "wzIL":
            self.f = self.focal_lens_of_wz()
            print(self.f, 'wz dependent focal lens')

        elif switch == "wzILtauChirp":
            self.f = self.focal_lens_of_wz_tauL()
            print(self.f, 'wz dependent focal lens, IL dependent with chirp')

        else:
            self.f = self.focal_lens_constant()
            print(self.f, 'constant focal length')

        return self.f

    def new_focal_position_constant_focal_length(self, index_z):
        # this often named v
        self.f = self.focal_lens_constant()
        AA = (self.q ** 2 / self.f) - self.z * (1 - self.z[index_z] / self.f)
        BB = (self.q ** 2 / self.f ** 2) + (1 - self.z[index_z] / self.f) ** 2
        return AA / BB


    # second model that is Vincenti like follows here
    def new_focal_position_wz_N_dependence_of_q(self, index_z):
        AA = (self.q ** 2 / self.f) - self.z[index_z] * (1 - self.z[index_z] / self.f)
        BB = (self.q ** 2 / self.f ** 2) + (1 - self.z[index_z] / self.f) ** 2
        # print(self.harmonic_number, 'new v', AA / BB)
        #print('single value for harmonic_number:', self.harmonic_number, 'new focal position', AA/BB)
        return AA / BB

    #this relies on q(N) as well
    def new_beam_waist_single_value(self, index_z):
        # print('initial beamwaist', self.w0, 'for single value and harmonic_number:', self.harmonic_number)
        v_single = self.new_focal_position_wz_N_dependence_of_q(index_z)
        new_wz_harmonic = ((1 - v_single / self.f) ** 2) + (1 / self.q ** 2) * (
                self.z[index_z] + v_single * (1 - (self.z[index_z] / self.f))) ** 2
        new_wz_harmonic = self.w0 * (new_wz_harmonic ** 0.5)
        #print("w0(N)", new_wz_harmonic, "focal position:" ,v_single)
        return new_wz_harmonic

    def new_beam_waist_single_value_vincenti(self, index_z):

        # print('initial beamwaist', self.w0, 'for single value and harmonic_number:', self.harmonic_number)
        #wz_harmonic = self.w0/self.harmonic_number
        wz_harmonic = self.beam_waist_of_z_harmonic_scaling_wN(0)
        self.q = self.q_initial_vincenti()
        v_single = self.new_focal_position_wz_N_dependence_of_q(index_z)
        new_wz_harmonic = ((1 - v_single / self.f) ** 2) + (1 / self.q ** 2) * (
                self.z[index_z] + v_single * (1 - (self.z[index_z] / self.f))) ** 2
        new_wz_harmonic = wz_harmonic  * (new_wz_harmonic ** 0.5)
        #print("w0(N)", new_wz_harmonic, "focal position:" ,v_single)
        return new_wz_harmonic

    def new_divergence_from_w0_new(self, w_new):
        return self.lambdaL / (self.harmonic_number * math.pi * w_new)

    def switch_sign(self, var):
        return -var

    def focal_lens_wz_tauL_z_scan(self):
        for x in range(0,len(self.z)):
            self.beam_waist_of_z_value(self.z[x])
            #print(self.q)
            #print(self.wz, "beamwaist", "for N:", self.harmonic_number, "at position z:", self.z[x])

            self.focal_lenght_wz_tauL[x]=self.focal_lens_of_wz_tauL()
            self.focal_lenght_wz[x] = self.focal_lens_constant()
            self.focal_lenght[x] = self.focal_lens_constant()

        plt.figure(100)
        plt.title("focal lens over z")
        plt.plot(self.z, self.focal_lenght_wz_tauL, label ="focal lens of z")
        plt.plot(self.z, self.focal_lenght, label ="focal lens const")
        plt.legend()

        return self.focal_lenght_wz_tauL, self.focal_lenght_wz, self.focal_lenght
    


    def save_data(self,array_x, array_y):
        result = np.hstack((array_x, array_y))
        print('saved data')
        np.savetxt("Fig3_div_overN_theo" + '_' + "20240531" + ".txt", result, delimiter=' ',
                   header='string', comments='',
                   fmt='%s')

# needs to be change over z with constant harmonic number
    def resulting_divergence_over_harmonic_number_Vincenti(self):

        result_w0_new = np.zeros([len(self.z)])
        result_div_harmonic_number = np.zeros([len(self.z)])
        #self.q = self.q_initial()/(self.harmonic_number **2)
        self.q = self.q_initial_vincenti()
        self.test_beam_waist2()
        #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx q xxxxxxxxx", self.q)
        self.focal_lens_wz_tauL_z_scan()

        for x in range(0, len(self.z)):

            self.f = self.focal_lenght_wz_tauL[x]
            # print(self.f, 'focal length')
            # print(self.q)
            #print(self.z[x],"z position in mm")
            result_w0_new[x] = self.new_beam_waist_single_value_vincenti(x)

            result_div_harmonic_number[x] = self.new_divergence_from_w0_new(result_w0_new[x])
            #print(result_div_harmonic_number[x], "divergence at this position of harmonic")
            
        self.save_data(self.z, result_div_harmonic_number)

        plt.figure(10)
        plt.title("N"+str(self.harmonic_number))
        plt.plot(self.z, result_w0_new/self.w0, label = "w0(N)/w0(L)")
        plt.xlabel( 'z in mm relative to focal position')
        plt.ylabel ( 'w0(N,z)/w0 [mm]')
        plt.legend()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.figure(9)
        plt.title("N"+str(self.harmonic_number)+" divergence")
        plt.plot(self.z, result_div_harmonic_number/(60/1500), label="div over N_vincenti_like")
        plt.xlabel ('z in mm relative to focal position')
        plt.ylabel ('Theta(N)/Theta(L) ')
        #plt.xlim(16.5, 28)
        #plt.ylim(1,10)
        #plt.ylim(1, 10)
        plt.title("chirped pulse")
        plt.yscale("log")
        plt.legend()

    def resulting_divergence_over_harmonic_number_julia(self):

        result_w0_new = np.zeros([len(self.z)])
        result_div_harmonic_number = np.zeros([len(self.z)])
        self.q = self.q_initial()
        self.focal_lens_wz_tauL_z_scan()
        for i in range(0, len(self.z)):

        # difference to the above model is: the function self.focal_length
        #
            self.f = self.focal_lenght_wz_tauL[i]
            print(self.f, 'focal length')
            print(self.z[i],"z position in mm")
            result_w0_new[i] = self.new_beam_waist_single_value(i)

            result_div_harmonic_number[i] = self.new_divergence_from_w0_new(result_w0_new[i])
            #print(result_div_harmonic_number[i], "divergence at this position of harmonic")

        plt.figure(9)
        plt.title("N"+str(self.harmonic_number)+" divergence")
        plt.plot(self.z, result_div_harmonic_number/(60/1500), label = "w(N)(z)=w(L)(z)")
        plt.xlabel ('z in mm relative to focal position')
        plt.ylabel ('Theta(N)/Theta(L)')
        #plt.xlim(16.5, 28)
        #plt.ylim(1,10)
        #plt.ylim(1, 10)
        plt.title("chirped pulse")
        plt.yscale("log")
        plt.legend()




    # index = list(zip(index[0])





#Denting has to be given in mm
# experimental beamwaist for fundamental: in mm
w0 = 0.011
D0 = 0.0001
Test = GaussianSecondLensSingleValue(w0, 0.0008, D0, 'tab:grey', marker = ".")
Test.resulting_divergence_over_harmonic_number_julia()
Test.resulting_divergence_over_harmonic_number_Vincenti()



#plt.xlabel('harmonic number N')
#plt.ylabel('divergence half angle 1/e [mrad]')

plt.figure(9)

plt.legend(bbox_to_anchor=(1.05, 1), loc=4, borderaxespad=0.)
plt.savefig("test_" + ".png", bbox_inches="tight", dpi=1000)
plt.show()

# Test3.resulting_divergence_over_N(2.5)
