# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:17:34 2020

@author: similiarities
"""
# Gaussian beam profile (focused by lens 1) which passes lens2 close to 
# focal position. The resulting beamwaist and new focal position 
# here only for inward denting and placed (last plots) before focus

import math
import numpy as np
import matplotlib.pyplot as plt


class Gaussian_second_lens:
    def __init__(self, w0, lambdaL, defocusing_range, Denting, N):
        
        print('xxxxxxxx run with Denting: ' +str(Denting) + 'xxxxxxx')
        self.w0 = w0
        self.w0_fundamental = w0
        self.N = N
        self.lambdaL = lambdaL #fundamental!!! 
        self.z = np.arange(-defocusing_range, defocusing_range, 0.01)
        self.v_array = np.zeros(len(self.z))
        self.w_new = np.zeros(len(self.z))
        #print(self.z)
        self.Denting = Denting
        self.f_array = np.zeros(len(self.z))
        self.wz_array = np.zeros(len(self.z))
        self.wz_array_N = np.zeros(len(self.z))
        self.denting_array = np.zeros(len(self.z))
        self.zr = self.Rayleigh_length()
        print(self.zr)
        self.Radius_array = np.zeros(len(self.z))
        self.beamdivergence_array = np.zeros(len(self.z))
        self.wz = float
        self.f = float
        self.q = self.q_initial()
        self.N_array = np.arange(1, 32, 1)
        self.focal_lens_of_wz()
        self.beam_waist_of_z_harmonic()
        
        

    def q_initial(self):

        self.q = (self.w0**2)*math.pi/(self.lambdaL/self.N)
        #print(self.N, self.q)
        return self.q
    
    def Rayleigh_length(self): 
        #remains Ry as fundamental 
        Ry = math.pi*(self.w0_fundamental**2)/(self.lambdaL)
        print('Rayleigh length: ' , Ry, 'for wo', self.w0_fundamental )
        return Ry
    

    
    def diffraction_limit(self):
        self.w0 = 1500*self.lambdaL*2/(self.N*math.pi*60)
        print('diff limit initial focal length for N:', self.N, 'is:', self.w0)
        self.q_initial()
        return self.w0
        
    
        

    def beamwaist_of_z(self):
        self.wz_array[::] = self.w0_fundamental *(1+(self.z[::]/self.zr)**2)**0.5
        print('initial beamwaist with 800nm:', self.w0_fundamental)
        
        plt.figure(11)
        plt.plot(self.z, self.wz_array, label = 'fundamental')
        plt.legend()
        

        return self.wz_array
    
    
    def beam_waist_of_z_harmonic(self):
        # definition of rayleigh length - which would scale with 1/N
        # can be used is not used now
        zr = math.pi*self.w0**2/(self.lambdaL/self.N)
        print('calculation of w(z) for N:', self.N, 'with new Ry(N)', zr)

        self.wz_array_N[::] = self.w0*(1+(self.z[::]/zr)**2)**0.5
        
        print('N:', self.N, 'Ry', zr, 'w0 inital', self.w0,'for calculation w(z) before lens2')

        name = 'N: ' + f'{self.N}'
        
        plt.figure(11)
        plt.plot(self.z, self.wz_array_N, label = name)
        plt.legend()
        

        return self.wz_array_N
        
        



    def Radius_of_lens(self):
        return (4*(self.Denting**2)+self.w0_fundamental**2)/(8*self.Denting)
    
    def Radius_of_lens_and_beamwaist(self):
        
        self.denting_array[::] = self.Denting * ((self.w0_fundamental)/(self.wz_array[::]))
        self.Radius_array[::] = (4*(self.denting_array[::]**2) + (self.wz_array[::]**2))/(8*self.denting_array[::]) 
 
        return self.Radius_array
    

    def focal_lens(self, Radius):
        return Radius/2
    
    def focal_lens_of_wz(self):
        
        self.beamwaist_of_z()
        self.Radius_of_lens_and_beamwaist()
        
        self.f_array[::] = self.Radius_array[::]/2
        
        name1 = 'focalL(w(z)), Dmax:' +str(self.Denting)
        
        plt.figure(3)
        plt.plot(self.z, self.f_array, label = name1)
        plt.xlabel='defocusing [mm]'
        plt.ylabel='focal length [mm]'
        plt.legend()
        return self.f_array

    def calculate_new_focal_position_constant_focal_length(self, z):
        # this often named v 
        self.f = self.focal_lens(self.Radius_of_lens())

        AA = (self.q**2/self.f)-z*(1-z/self.f)
        BB = (self.q**2/self.f**2) + (1-z/self.f)**2
        return AA/BB

    def calculate_new_focal_position_beamwaist_dependent(self, index):
        i = index
        
        AA = (self.q**2/self.f_array[i])-self.z[i]*(1-self.z[i]/self.f_array[i])
        BB = (self.q**2/self.f_array[i]**2) + (1-self.z[i]/self.f_array[i])**2
        #print('single value for N:', self.N, 'new focal position', AA/BB)
        return AA/BB
    
    def calulate_beam_waist_single_value(self, index):
       # print('initial beamwaist', self.w0, 'for single value and N:', self.N)
        
        v_single = self.calculate_new_focal_position_beamwaist_dependent(index)
        new = ((1 - v_single/self.f_array[index])**2) + (1 /self.q **2) *(self.z[index]+v_single *(1-(self.z[index]/self.f_array[index])))**2
        new =  self.w0*(new**0.5)
        return new
    
    def calculate_div_from_w0_new(self, w_new):
        return (self.lambdaL/(self.N *math.pi *w_new))
        
    

    
    def create_v_array_constant_focal_lens(self):
        
        self.f_array[::] = self.f
        
        for i in range(0, len(self.z)):
           self.v_array[i] = self.calculate_new_focal_position_constant_focal_length(self.z[i])
        name = 'f contst.: ' + str(round(self.f, 3))
        plt.figure(1)
        plt.plot(self.z, self.v_array, label = name)
        plt.xlabel='defocusing [mm]'
        plt.ylabel='new focal position [mm]'
        plt.title('new focal position')
        plt.legend()
        
        plt.figure(2)
        name2 = 'f: ' + str(round(self.f,3))
        plt.hlines(self.f, self.z[0], self.z[-1], label =  name2)
        
        plt.xlabel='defocusing [mm]'
        plt.ylabel='focal length [mm]'
        plt.legend()
        return self.v_array
    
    def create_v_array_z_dependent_focal_lens(self):

        
        
        for i in range(0, len(self.z)):

            self.v_array[i] = self.calculate_new_focal_position_beamwaist_dependent(i)

        name1 = 'f(w(z)) for Dmax: ' +str(self.Denting) + 'N: ' + str(self.N)
        name2 = 'f(w(z)) for Dmax: ' +str(self.Denting)
        
           
        plt.figure(1)
        plt.plot(self.z, self.v_array, label = name1)
        plt.xlabel='defocusing [mm]'
        plt.ylabel='new focal position [mm]'
        plt.legend()
        plt.savefig("caseII_zN_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        
        plt.figure(2)
        plt.plot(self.z, self.f_array, label = name2)
        plt.ylabel='focal lens [mm]'
        plt.xlabel='defocusing [mm]'
        
        plt.legend()
        #plt.savefig("caseII_f_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        return self.v_array
    
    
    def calculate_new_beamwaist(self):
        print('for N: ', self.N, 'we start with a beamwaist for lens 1 of about:' , self.w0)
        self.create_v_array_z_dependent_focal_lens()
        
            #print(self.f, 'insert in wnew')
        self.w_new[::] = ((1-self.v_array[::]/self.f_array[::])**2)+(1/self.q**2)*(self.z[::]+self.v_array[::]*(1-(self.z[::]/self.f_array[::])))**2
        self.w_new[::] =  self.w0*(self.w_new[::]**0.5)
        
        plt.figure(4)
        name = 'new beamwaist of f(w(z)), D0: ' +str(self.Denting) + 'N: ' + str(self.N)
        plt.plot(self.z, self.w_new, label = name)
        plt.xlabel='defocusing [mm]'
        plt.ylabel='w_new(z) [mm]'
        plt.legend()
        #plt.savefig("caseII_wN_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
        return self.w_new
    
    
    def switch_sign(self, var):
        return -var
    
    def resulting_beam_divergence(self):
        
        self.beamdivergence_array[::] = (self.lambdaL/(self.N *math.pi*self.w_new[::]))
        
        name = 'Theta Dmax:' +str(self.Denting) + 'N: ' + str(self.N)
        
        plt.figure(6)
        plt.plot(self.z, self.beamdivergence_array, label = name)
        plt.yscale('log')
        #plt.ylim(1,6)
        plt.xlabel='defocusing [mm]'
        plt.ylabel='[rad]'
        plt.legend()
        #plt.savefig("caseII_Theta_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)
     
        
        
        
        
    def resulting_divergence_over_N(self, z):
        
        #intensity (position) dependent focal length lens 2
 
        index = list(zip(*np.where(self.z >= z)))
        index = index[0]
        #print(index, self.z[index])
        result_w0_N = np.zeros(len(self.N_array))
        result_div_N = np.zeros(len(self.N_array))
        
        for x in range(0,len(self.N_array)):
            
            
            
            self.N = self.N_array[x]
            self.q_initial()
            

            result_w0_N[x] = self.calulate_beam_waist_single_value(index)

            result_div_N[x] =self.calculate_div_from_w0_new(result_w0_N[x])
            
        name1 = 'z: ' + str(z) +'[mm]  lens+' 
        
            
        plt.figure(10)
        plt.plot(self.N_array, result_w0_N, label = name1+'w(N)')
        plt.xlabel='N'
        plt.ylabel='w0(N)* in [mm]'
        plt.legend()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        
        plt.figure(9)
        plt.plot(self.N_array, result_div_N, label = name1 + 'div detctor', marker = '.')
        plt.xlabel = 'N'
        plt.ylabel='div in [rad]'
        plt.legend()
        #plt.yscale ('log')
        #plt.ylim(0.001,0.006)


        

         	
        # zip the 2 arrays to get the exact coordinates
       # index = list(zip(index[0])
    def plot_diffraction_limit(self):
        N_list =np.arange(1, 30, 1)
        N_diffraction_limit = np.zeros([29,1])
        for x in range(0, 30-1):

                #halfangle
            N_list[x] = 1 + x
            N_diffraction_limit[x] = (60./1500.)/(1+x)

                
        plt.figure(9)
        plt.scatter(N_list, N_diffraction_limit, marker = "o", color ="c", label = "Theta_L/N")
        plt.hlines(0.007,0,30,label="detector limit")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

            #plt.savefig("20190123_divergence_mrad_halfangle_and_theo" +".png",  bbox_inches="tight", dpi = 1000)

    
    
        
    
    
    
        
        


Test = Gaussian_second_lens(0.012, 0.0008, 5, 0.00005, 1)
#Test.beamwaist_of_z()
#Test.Radius_of_lens_and_beamwaist()
#Test.create_v_array_z_dependent_focal_lens()
#Test.q_harmonic()
Test.calculate_new_beamwaist()
Test.resulting_beam_divergence()

Test2 = Gaussian_second_lens(0.012, 0.0008, 5, 0.00005, 25)
#Test.q_harmonic()

Test2.calculate_new_beamwaist()
Test2.resulting_beam_divergence()


Test3 = Gaussian_second_lens(0.012, 0.0008, 5, 0.00005, 2)
#Test.q_harmonic()
Test3.calculate_new_beamwaist()
Test3.resulting_beam_divergence()

Test3.resulting_divergence_over_N(0.0)
Test3.resulting_divergence_over_N(-0.5)
#Test3.resulting_divergence_over_N(2.0)
Test3.resulting_divergence_over_N(-1.0)
Test3.resulting_divergence_over_N(-2.0)
Test3.resulting_divergence_over_N(-2.5)
Test3.plot_diffraction_limit()
plt.savefig("caseII_div_over_N_50nm" +".png",  bbox_inches="tight", dpi = 1000)




plt.legend()
plt.xlim(1,28)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()



#Test3.resulting_divergence_over_N(2.5)

