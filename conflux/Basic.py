import numpy as np
import csv

class Spectrum:
    """A general class of a spectrum to initialize, add, and modify the object."""
    
    def __init__(self, xbins=np.arange(0, 20, 0.1)):
        self.xbins = xbins
        """The particle energy (MeV), or x values of the spectrum"""
        self.nbin = len(xbins)
        """Size of the spectrum array"""
        self.spectrum = np.zeros(self.nbin)
        """Contents of the spectrum"""
        self.uncertainty = np.zeros(self.nbin)
        """Spectrum content uncertainties"""
        self.integral = 0
        """Integral of spectrum"""
        
    def Add(self, targetSpec, W=1, sigma_W = 0):
        '''
        An operator to add (or subtract) a spectrum
        
        Parameters:
            targetSpec (list): A list containing the spectrum we want to add to the current spectrum
            W (float): The weight of the added spectrum
            sigma_W (float): The uncertainty in the weight of the added spectrum
        Returns:
            None        
        '''
        assert(targetSpec.xbins == self.xbins) # ensure the identical binning
        targetSpec.Norm(W, sigma_W) #Rescale the added spectrum to the current spectrum
        self.spectrum += targetSpec.spectrum #Add the target spectrum to the current spectrum
        
    def Norm(self, W, sigma_W = 0):
        '''
        An operator to scale the spectrum by a weighting factor

        Parameters:
            W (float): The weight of the spectrum
            sigma_W (float): The uncertainty in the weight
        Returns:
            None
        '''
        uncertainty = self.uncertainty
        #Calculate the relative uncertainty of the spectrum
        relativeUnc = np.sqrt((sigma_W/W)**2+(uncertainty*self.spectrum)**2)
        #calculate the normed uncertainty        
        self.uncertainty = W*self.spectrum*relativeUnc
        #Scale the spectrum by our Normalization weight
        self.spectrum *= W
    
    def Integral(self):
        '''
        The absolute integral of the spectrum
        
        Parameters:
            None
        Returns:
            integral (float): The integral of the spectrum
        '''
        #Calculate the integral of the spectrum
        self.integral = np.dot(self.spectrum.sum(), self.xbins)
        return self.integral
        
    def SaveToFile(self, filename):
        '''
        To save the spectrum to a csv file

        Parameters:
            filename (String): The name of the file you want to save the spectrum to.
        Returns:
            None
        '''

        #Save the spectrum information in the file in the format [Energy, spectrum, error]
        with open(filename, 'w', newline='') as outputfile:
            colNames = ['E', 'count', 'error']
            writer = csv.DictWriter(outputfile, fieldnames=colNames)
            writer.writeheader()
            for i in range(self.nbin):
                E = self.xbins[i]
                count = self.spectrum[i]
                error = self.uncertainty[i]
                writer.writerow({'E':E, 'count':count, 'error':error})


class Summed:
    '''
    A general class for summed objects, such as summed spectra

    ...

    Parameters
    ----------

    ID : int
        The isotopic ID of a particular spectrum in the summed spectrum
    branches : dictionary
        A dictionary that contains the branch IDs of all spectra in this summed object
    spectra : dictionary
        A dictionary containing all spectral information for each branch
    uncertainty : dictionary
        A dictionary containing all uncertainty information for each branch
    
    Methods
    -------
    Clear():
        Clears out all associated dictionaries inside the summed spectra
    AddBranch(branch):
        Adds a branch object to the summed specctra
    SumSpectra():
        Sums the spectra of all branches
    '''
    def __init__(self):
        self.ID = 0
        self.branches = {}
        self.spectra = {}
        self.uncertainty = {}
    
    def Clear(self):
        """
            Clears out all associated dictionaries inside the summed spectra

            Parameters:
                None
            Returns:
                None
        """

        #Pretty self explanatory what these do
        self.branches = {}
        self.spectra = {}
        self.uncertainty = {}
    
    def AddBranch(self, branch):
        '''
            Adding a branch element to the summing object

            Parameters:
                branch (int): The ID of the branch you want to add to the summed spectrum
            Returns:
                None
        '''
        self.branches[branch.ID] = branch
    
    def SumSpectra(self):
        '''
            Summing units together

            Paremeters:
                None
            Returns:
                None
        '''
        
        spectrum = Spectrum()
        #walk through all branch IDs in the branch dictionary
        for ZAI in self.branches:
            #Add all the individual branch spectra to the total spectra
            spectrum.Add(self.spectra[ZAI], self.branches[ZAI])
