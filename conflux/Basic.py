import numpy as np
import csv

class Spectrum:
    '''
    A general class of a spectrum to initialize, add, and modify the object
    '''
    def __init__(self, binwidths=0.1, spectRange=[0.0, 20.0]):
        '''
        initializer
        '''
        self.binwidths=binwidths
        self.bins = int(spectRange[1]/binwidths)
        self.xbins = np.arange(*spectRange, binwidths)
        self.spectrum = np.zeros(self.bins)
        self.uncertainty = np.zeros(self.bins)
        self.integral = 0
        
    def Add(self, targetSpec, W=1, sigma_W = 0):
        '''
        An operator to add (or subtract) a spectrum
        '''
        assert(targetSpec.xbins == self.xbins) # ensure the identical binning
        targetSpec.Norm(W, sigma_W)
        self.spectrum += targetSpec.spectrum
        
    def Norm(self, W, sigma_W = 0):
        '''
        An operator to scale the spectrum by a weighting factor
        '''
        uncertainty = self.uncertainty
        self.uncertainty = W*self.spectrum*np.sqrt((sigma_W/W)**2+(uncertainty*self.spectrum)**2)
        self.spectrum *= W
    
    def Integral(self):
        '''
        The absolute integral of the spectrum
        '''
        self.integral = self.spectrum.sum()*self.binwidths
        return self.integral
        
    def SaveToFile(self, filename):
        '''
        To save the spectrum to a csv file
        '''
        with open(filename, 'w', newline='') as outputfile:
            colNames = ['E', 'count', 'error']
            writer = csv.DictWriter(outputfile, fieldnames=colNames)
            writer.writeheader()
            for i in range(self.bins):
                E = self.xbins[i]
                count = self.spectrum[i]
                error = self.uncertainty[i]
                writer.writerow({'E':E, 'count':count, 'error':error})


class Summed:
    '''
    A general class for summed objects, such as summed spectra
    '''
    def __init__(self):
        self.ID = 0
        self.branches = {}
        self.spectra = {}
        self.uncertainty = {}
    
    def Clear(self):
        """
            Clears out all associated dictionaries inside the SumEngine

            Parameters:
                None
            Returns:
                None
        """
        self.branches = {}
        self.spectra = {}
        self.uncertainty = {}
    
    def AddBranch(self, branch):
        '''
        Adding a branch element to the summing object
        '''
        self.branches[branch.ID] = branch
    
    def SumSpectra(self):
        '''
        Summing units together
        '''
        spectrum = Spectrum()
        for ZAI in self.branches:
            spectrum.Add(self.spectra[ZAI], self.branches[ZAI])
