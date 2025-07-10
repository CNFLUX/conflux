#Author: Anosh Irani

#CEvNS cross section model
import numpy as np
import scipy as sp
import scipy.constants as const
import math
import matplotlib.pyplot as plt

#Fermi Coupling constant, converted to MeV
G_F, G_F_unit, G_F_unc = const.physical_constants["Fermi coupling constant"]
G_F = G_F * (1e-3)**2


#Weak mixing angle
theta_w, theta_w_unit, theta_w_unc = const.physical_constants["weak mixing angle"]


#Atomic mass unit, to convert nucleons to MeV
amu, amu_unit, amu_unc = const.physical_constants["atomic mass constant energy equivalent in MeV"]

#MeV^-2 to cm^2 conversion factor 1 MeV^-2 = x cm^2
MeV_to_cm = (197.3**2) / (1e26)



class CEvNS():
    """
    A class to calculate the CEvNS spectrum on a specific target
   
    """
    def __init__(self, Z, N, E_r, E_v):
        self.N = N
        self.Z = Z
        self.M = (N + Z) * amu
        self.E_r = E_r
        self.E_v = E_v
        self.min_E_v()
        
        
   # minimum E_v = (T + sqrt(T**2 + 2*M*T))/2
    def min_E_v(self):
        self.E_v_min = []
        for i in self.E_r:
            self.E_v_min.append(float(0.5 * (i + math.sqrt((i**2) + (2 * self.M * i)))))

    # max E_r = 2E_v**2 / (M + 2E_v)
    def max_recoil(self):
        self.E_r_max = []
        for i in self.E_v:
            self.E_r_max.append(2 * i**2 / (self.M + 2*i))
            

    def Huber_sevens(self):
        
        self.diff_xsec = np.zeros((len(self.E_v), len(self.E_r)))
        self.min_E_v()
        
        constant = G_F**2 / (4 * np.pi)
        constant *= self.N**2 * self.M
        for i in range(len(self.E_v)):
            for j in range(len(self.E_r)):
                if(self.E_v[i] < self.E_v_min[j]):
                    self.diff_xsec[i][j] = 0
                    continue
                
                term = (self.M * self.E_r[j]) / (2 * self.E_v[i]**2)
                
                total = constant * (1 - term)

                if (math.isnan(total)):
                    print("isNan")
                    self.diff_xsec[i][j] = 0
                    continue
                elif(total < 0):
                    print("Is negative")
                    self.diff_xsec[i][j] = 0
                    continue                   
                else:
                    self.diff_xsec[i][j] =  total
            
    #https://arxiv.org/pdf/2406.16081
    #Method to Calculate the differential CEvNS x-section, and put it into a
    #Matrix of values.
    def diff_sevens(self, mag_moment = False):
        
        #Nuclear Mass in MeV
        M = (self.N + self.Z) * amu
        
        #Weak charge
        Q_w = (0.5 - 2 * theta_w) * self.Z - 0.5 * self.N

        #constant term
        const = ((G_F**2 * M) / np.pi) * (Q_w**2)


        self.diff_xsec = np.zeros((len(self.E_v), len(self.E_r)))
        self.min_E_v()
        for i in range(len(self.E_v)):
            for j in range(len(self.E_r)):
                #Each recoil energy will have an associated minimum 
                #Neutrino energy that can induce that recoil. If 
                #Neutrino energy is less than that minimum neutrino energy,
                #Set x-section to 0 for that E_v and E_r value. 
                if(self.E_v[i] < self.E_v_min[j]):
                    self.diff_xsec[i][j] = 0
                    continue


                first_term = (M * self.E_r[j]) / (2 * self.E_v[i]**2)
                second_term = self.E_r[j] / self.E_v[i]
                third_term = 0.5 * (self.E_r[j] / self.E_v[i])**2
                total = 1 - first_term - second_term + third_term
                if (math.isnan(total)):
                    print("isNan")
                    self.diff_xsec[i][j] = 0
                    continue
                elif(total < 0):
                    print("Is negative")
                    self.diff_xsec[i][j] = 0
                    continue                   
                else:
                    self.diff_xsec[i][j] = const * total

    def integrate_xsec_r(self, min_recoil_threshold):
        integrate = np.zeros(len(self.E_v))

        
        if (min_recoil_threshold < self.E_r[0] 
            or min_recoil_threshold > self.E_r[len(self.E_r) - 1]):
            print("Your threshold is out of bounds of the recoil spectrum")
            print("Please input a valid recoil threshold energy")
            return integrate
        
        self.max_recoil()
                
        recoil_bin_width = self.E_r[2] - self.E_r[1]
        index_min = np.where(np.isclose(self.E_r, min_recoil_threshold))[0][0]
        
        for i in range(len(self.E_v)):
            integrate[i] = sum(self.diff_xsec[i][index_min:]) * recoil_bin_width
        
        integrate *= MeV_to_cm
        
        return integrate
    
    def plot(self,figname, threshold = False):
        binwidth_E_v = self.E_v[1] - self.E_v[0]
        binwidth_E_r = self.E_r[1] - self.E_r[0]
        lower_e_v = self.E_v[0]
        upper_e_v = self.E_v[len(self.E_v) - 1]
        lower_e_r = self.E_r[0] * 1e6
        upper_e_r = self.E_r[len(self.E_r) - 1] * 1e6

        fig, ax = plt.subplots()
        im = ax.imshow(self.diff_xsec, origin = "lower", cmap = "Wistia", aspect="auto", extent=[lower_e_r, upper_e_r, lower_e_v , upper_e_v])


        index = int(1.8 / binwidth_E_v)
        if(threshold):

            for j in range(len(self.E_r)):
                if self.diff_xsec[int(index)][j] == 0:
                    ax.axvline(x = j, color = 'black', linestyle="--")
                    ax.axhline(y = self.E_v[index], color = "black", linestyle="--")
                    break

        fig.colorbar(im)
        ax.set_ylabel("E_v (MeV)")
        ax.set_xlabel("E_r (eV)")
        fig.savefig(str(figname) + ".png")
        fig.clf()




