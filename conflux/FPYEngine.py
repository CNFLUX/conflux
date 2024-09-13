# fission product yield engine:
# input: dictionary of fission fractions, with {'ZA', fission_fraction}
# output: dictionary of fission yield isotopes, with {'ZA', isotope_fraction}

# universal modules
import sys
import os
import csv
import copy
import numpy as np
from os import listdir
import xml
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import pkg_resources
from tqdm import tqdm

from conflux.config import CONFLUX_DB
from conflux.Basic import *


# this class saves nuclide info of fission products
class FPNuclide:
    """
    Class to handle all the Fission Nuclide Information.

    ...

    Attributes
    ----------
    Z : (int)
        Atomic number of your nuclide
    A : (int)
        atomic mass of your nuclide
    N : (int)
        Neutron count in your nuclide
    isomer : (int)
        The isomer number of your nuclide
    FPZAI : (dictionary)
        Fission Product dictionary for this nuclide. Is denoted with it's ZAI number
    y : (float)
        The yield of the specific nuclide
    cov : (dictionary)
        The covariance matrix for the specific nuclide
    corr : (dictionary)
        The correlation matrix for the specific nuclide
    yerr : (float)
        The yield error for this specific nuclide

    Methods
    -------
    Contribute(fraction, d_fraction=0):
        Adds the fission yield of a specific branch to the overall fission yield
    AddCovariance(newNuclide):
        Add the covariance matrices of different nuclides together
    """
    def __init__(self, FPZAI, y, yerr):
        self.Z = int(FPZAI/10000)
        self.A = int((FPZAI-self.Z*10000)/10)
        self.N = self.A-self.Z
        self.isomer = int(FPZAI-self.Z*10000-self.A*10)
        self.FPZAI = FPZAI
        self.y = y
        self.cov = {}
        self.corr = {}
        self.yerr = yerr

    # Method that add fission yield of this branch to total fission yield.
    # When use_corr == True, assume correlation matrix is loaded, calculate a
    # new covariance matrix.
    def Contribute(self, fraction, d_fraction=0, use_corr = False):
        """
           Adds the fission yield of this branch to the total fission yield

            Parameters:
                fraction (float) : The fractional contribution of this nuclide to the overall yield
                d_fraction (float) : the uncertainty in the fractional contribution of this nuclide
            Returns:
                None
        """

        #
        # if use_corr:
        #     assert self.corr, "Correlation matrix was not loaded!"
        #     self.cov = {}
        #     for element in self.corr:
        #         self.cov

        # Adding the fission fraction uncertainty and fission yield uncertainty together
        self.yerr = self.y*fraction*np.sqrt((self.yerr/self.y)**2 + (d_fraction/fraction**2))
        # force yeild uncertainty to equal 0, when yeild is zero
        self.yerr = np.nan_to_num(self.yerr, nan=0.0)

        # Also scale the covariance matrix.
        for key in self.cov:
            if key == self.FPZAI:
                self.cov[key] = self.y*fraction*np.sqrt((self.cov[key]/self.y)**2 + (d_fraction/fraction**2))
            else:
                self.cov[key] *= fraction

        self.y *= fraction  # multiply fission yield with fission fraction

    # Method to add covariance matrices together
    def AddCovariance(self, newNuclide):
        """
            Adds covariance matrices together

            Parameters:
                newNuclide (FPNuclide) : The nuclide whose covariance matrix you want to add to this one.
            Returns:
                None
        """
        # print("Added FPY covaraince matrix")
        for key in newNuclide.cov:
            if key not in self.cov:
                self.cov[key] = newNuclide.cov[key]
            else:
                self.cov[key] += newNuclide.cov[key]

            #if key == self.FPZAI: print(key, self.cov[key])

# Class that counts fission products of a specified fission isotope
class FissionIstp(Spectrum, Summed):
    """
    Class to handle all the Fission Nuclide Information.

    ...

    Attributes
    ----------
    Z : (int)
        Atomic number of your nuclide
    A : (int)
        Atomic mass of your nuclide
    FPYlist : (dictionary)
        Dictionary of cumulative or independent fission yields {"FPZAI", FPNuclide}
    DBTitle : (dictionary)
        Dictionary of the two fission databases included in CONFLUX (can be downloaded using the file in the aux folder)

    Methods
    -------
    LoadFissionDB(DB='ENDF'):
        Function to load data from a specific fission database. DB has format /path/to/DB. Defaults to ENDF
    LoadCovariance(DB='ENDF', percent=True):
        Load the covariance matrices associated with the fission isotopes. DB has format /path/to/DB. Defaults to ENDF
    """

    def __init__(self, Z, A, Ei, DB='ENDF', IFPY=False):
        self.Z = Z
        self.A = A
        self.Ei = Ei
        self.IFPY=IFPY
        self.FPYlist = {}
        self.DBtitle = {'ENDF':'nfy', 'JEFF':'nfpy'}

        self.LoadFissionDB(Ei=self.Ei, DB=DB)

    # method that load xml database of FPY and save nuclide info in dictionaries.
    def LoadFissionDB(self, Ei=None, DB='ENDF'):
        """
           Loads the fission products from a given database into the Independant and cumulative fission dictionaries in the model.

            Parameters:
                DB (String) : The path to a user inputed fission database. Has the form "/path/to/database
            Returns:
                None
        """

        if Ei != None: self.Ei = Ei
        DBname = DB
        # Check if the user gave a valid Database path
        # if there is no specified fissionDB, look for the default one
        if DBname == None or not os.path.exists(DBname):
            DBpath = CONFLUX_DB+"/fissionDB/"+DB+"/"
            if DBname != None and DBname not in list(self.DBtitle.keys()):
                print('Custom DB: '+ DBname + ' NOT found!')
            print('Reading default FPY DB from folder: '+DBpath+'...')
            fileList = listdir(DBpath) #Get the list of files in the Database directory
            istpfound = False
            for filename in fileList: #iterate through the list of files in the directory
                namecache = filename.split('.')
                if namecache[-1] != 'xml': #if the file is not an xml file, continue
                    continue
                if (self.DBtitle[DB] not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                    continue # Alternatively, if the isotope is not in the list, continue
                istpfound = True # Else, assert that the isotope is in the DB
                break
            if (not istpfound):
                print(f"WARNING!!! Fisson DB {DBpath} {filename} not found")# assert error if isotope not found in DB

                return

            DBname = DBpath+filename # this is the isotope that we found in the list.
            assert(DBname)

        print('Reading FPY DB: '+DBname+'...')
        tree = ET.parse(DBname)
        root = tree.getroot()

        Ei_exist = False
        Ei_list = []
        for HEAD in root:
            MT = HEAD.attrib['MT']
            # "MT" indicates decay data type, designated by the ENDF-6 format
            # In the case of fission products, there are cumulative fission
            # fission products (CFP) and indpendent fission products (IFP)
            if (MT == 'CFP') == self.IFPY:
                continue

            for LIST in HEAD: # begin to read fission modes of each fissile isotope
                Ei = float(LIST.attrib['Ei'])   # Induced neutron energy
                Ei /= 1e6                       # convert to MeV
                if (Ei < 0.01): Ei = 0          # simplify thermal neutron E
                Ei_list.append(Ei)
                if Ei != self.Ei: continue
                Ei_exist = True
                branch = int(LIST.attrib['NFPi'])
                assert(branch == len(LIST))
                nuclidelist = {}

                for CONT in LIST:
                    FPZA = int(CONT.attrib['ZA'])
                    isomeric = float(CONT.attrib['FPS'])
                    Y = float(CONT.attrib['Y'])
                    DY = float(CONT.attrib['DY'])
                    #print(MT, Ei, FPZA, Y)
                    FPZAI = FPZA*10+isomeric
                    nuclide = FPNuclide(FPZAI, y=Y, yerr=DY)
                    nuclidelist[FPZAI] = nuclide
                    self.FPYlist = nuclidelist

                FPYsize = len(nuclidelist)

                print(f"FPY list loaded. {FPYsize} fission products were found in "+MT+f" DB with Ei of {Ei} MeV ")

        assert Ei_exist, 'Isotope in '+DBname+' has only fission data with Ei = '+str(Ei_list)+' MeV! (input Ei = ' +str(self.Ei)+')'

    # Method to read the prepackaged covariance csv file
    # This function has to be called after loading the fission DB for neutrino
    # flux calcuation.
    def LoadCovariance(self, DB = 'ENDF', percent=True):
        """
           Loads the covariance matrices from a given database into the Independant and cumulative fission dictionaries in the model.

            Parameters:
                DB (String) : The path to a user inputed covariance file. Has the form "/path/to/database
                percent (boolean) : determines whether the covariances are relative or not.
            Returns:
                None
        """

        DBpath = DB
        #Figure out where the DB and files are in a similar method we loaded the
        #Fission Database.
        if DBpath in ['JEFF', 'ENDF'] or not os.path.exists(DBpath):
            DBpath = CONFLUX_DB+"/fissionDB/"+DB+"/"
            print("Reading covariance matrices in: "+DBpath+"...")
        fileList = listdir(DBpath)
        assert(DBpath)
        istpfound = False
        #By default, set Neutron energy to 0 (thermal), 0.5 (fast), or 14 (relativistic) << default DB is ENDF
        e_neutron = {'T': 0, 'F': 0.5, 'H': 14}
        #If JEFF is chosen, adjust the neutron energies.
        if DB =='JEFF':
            e_neutron = {'T': 0, 'F': 0.4, 'H': 14}
        filesfound = []
        #Check for database files, and append them into a files list.
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'csv':
                continue
            # if ("normed_cov" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
            #     continue
            if ("cov" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            print("Filename is ", filename)
            filesfound.append(filename)

        # determine the whether the covariance are relative or absolute
        rate = 1e4 if percent else 1
        for mode in e_neutron:
            for filename in filesfound:
                if (mode not in filename) or self.Ei !=e_neutron[mode]:
                    continue#If our matrix is not defined at a given neutron energy, skip.

                #Read through the file, pull out the relevant FPZAI number of this isotope
                DBname = DBpath+filename
                print('Reading covariance data: '+DBpath+filename+'...')
                with open(DBname) as inputfile:
                    reader = csv.DictReader(inputfile, dialect='excel', delimiter=',')

                    for row in reader:
                        row_id = int(row[''])
                        z = int(row_id/10000)
                        i = int((row_id-z*10000)/1000)
                        a = int(row_id-z*10000-i*1000)
                        fpzai = z*10000+a*10+i
                        if fpzai not in self.FPYlist:
                            continue
                        #If the isotope is not in the Cumulative fission dictionary
                        #At the given energy, skip
                        for corrzai in self.FPYlist:
                            col_id = int(corrzai)
                            z = int(col_id/10000)
                            a = int((col_id-z*10000)/10)
                            i = int(col_id-z*10000-a*10)
                            key = z*10000+a+i*1000
                            keystr = ' '+str(key)

                            # if key is not found in the covaraince matrix, set
                            # value to zero
                            if keystr in row:
                                self.FPYlist[fpzai].cov[corrzai] = float(row[keystr])/rate
                            else:
                                self.FPYlist[fpzai].cov[corrzai] = 0.0

                    for nuclide in self.FPYlist:
                    # if element is not found in the covariance matrix, add the element into the matrix and set the diagonal element to 'yerr'
                        if nuclide not in self.FPYlist[nuclide].cov:
                            self.FPYlist[nuclide].cov = {otherNuclide: 0 for otherNuclide in self.FPYlist}
                            self.FPYlist[nuclide].cov[nuclide] = self.FPYlist[nuclide].yerr**2

        if not filesfound:
            print(f"WARNING!!! FPY covariance matrix {DBpath} {filename} not found")

    # Method to read the prepackaged correlation csv file
    # This function has to be called after loading the fission DB for neutrino
    # flux calcuation.
    #Most of the comments will be the same for both this and for the LoadCovariance function
    def LoadCorrelation(self, DB='ENDF'):
        """
           Loads the correlation matrices from a given database into the Independant and cumulative fission dictionaries in the model.

            Parameters:
                DB (String) : The path to a user inputed correlation file. Has the form "/path/to/database
            Returns:
                None
        """

        DBpath = DB  #Figure out where the DB and files are in a similar method we loaded the
        #Fission Database.
        if DBpath in ['JEFF', 'ENDF'] or not os.path.exists(DBpath):
            DBpath = CONFLUX_DB+"/fissionDB/"+DB+"/"
            print("Reading correlation matrices in: "+DBpath+"...")
        fileList = listdir(DBpath)
        assert DBpath, DBpath+" is not found"
        istpfound = False
        #By default, set Neutron energy to 0 (thermal), 0.5 (fast), or 14 (relativistic) << default DB is ENDF
        e_neutron = {'T': 0, 'F': 0.5, 'H': 14}
        if DB=='JEFF':  #If JEFF is chosen, adjust the neutron energies.
            e_neutron = {'T': 0, 'F': 0.4, 'H': 14}
        filesfound = []

        #Check for database files, and append them into a files list.
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'csv':
                continue #If the file you are looking at is not a csv covariance file, skip
            if ("corr" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            filesfound.append(filename) # assert error if isotope not found in DB

        for mode in e_neutron:
            for filename in filesfound:
                if (mode not in filename) or self.Ei !=e_neutron[mode]:
                    continue

                #Read through the file, pull out the relevant FPZAI number of this isotope
                DBname = DBpath+filename
                print('Reading correlation data: '+DBpath+filename+'...')
                with open(DBname) as inputfile:
                    reader = csv.DictReader(inputfile, dialect='excel', delimiter=',')
                    for row in reader:
                        row_id = int(row[''])
                        z = int(row_id/10000)
                        i = int((row_id-z*10000)/1000)
                        a = int(row_id-z*10000-i*1000)
                        fpzai = z*10000+a*10+i
                        if fpzai not in self.FPYlist:
                            continue #If the isotope is not in the Cumulative fission dictionary
                            #At the given energy, skip
                        for corrzai in self.FPYlist:

                            col_id = int(corrzai)
                            z = int(col_id/10000)
                            a = int((col_id-z*10000)/10)
                            i = int(col_id-z*10000-a*10)
                            key = z*10000+a+i*1000
                            keystr = ' '+str(key)
                            # if key is not found in the correlation matrix, set
                            # value to zero
                            if keystr in row:
                                self.FPYlist[fpzai].corr[corrzai] = float(row[keystr])
                            else:
                                self.FPYlist[fpzai].corr[corrzai] = 0.0

                    # if element not found, set diagonal element to 1.0
                    for nuclide in self.FPYlist:
                        if nuclide not in self.FPYlist[nuclide].corr:
                            self.FPYlist[nuclide].corr = {otherNuclide: 0 for otherNuclide in self.FPYlist}
                            self.FPYlist[nuclide].corr[nuclide] = 1.0

                self.CalcCovariance() #Once you've loaded the correlation database, calculate the
                    #subsequent covariance matrix from those correlations.
        if not filesfound:
            print(f"WARNING!!! FPY correlation matrix {DBpath} {filename} not found")  # assert error if isotope not found in DB

    # Method to calculate a covariance matrix from the correlation information
    # of each fission product at a given neutron energy.
    def CalcCovariance(self):
        """
            Calculates the covariance matrix from inputted correlation information
            Returns:
                None
        """

        #Iterate over all fission nuclides inside the Cumulative fission dictionary
        for i in self.FPYlist:
            yerri = self.FPYlist[i].yerr #Load the yield error of the nuclide i
            corr = self.FPYlist[i].corr  #Load the correlation of the nuclide i
            for j in self.FPYlist:
                yerrj = self.FPYlist[j].yerr #Load the yield error of nuclide j
                if j not in corr:
                    corr[j] = 0 if j!=i else 1
                self.FPYlist[i].cov[j]=yerri*corr[j]*yerrj #carry out the covariance calculation
                #(error of i) * (correlation of j) * (error of j)

    def CalcBetaSpectra(self, betaSpectraDB, processMissing=False, ifp_begin = 0, ifp_end = 0, modelunc = True, silent = False):

        self.xbins = betaSpectraDB.xbins
        nbins = len(self.xbins)
        betaSpectraList = {}
        betaUncertainty = {}
        #initialize total spectrum & uncertainty
        self.spectrum = np.zeros(nbins)
        self.uncertainty = np.zeros(nbins)

        #Initialize model and Yield Uncertainties, a list of missing branches to the total
        #Contribution, the missing yield, and the total yield.
        self.modelUnc = np.zeros(nbins)
        self.yieldUnc = np.zeros(nbins)

        # Find common isotopes among beta decaying isotopes and fission products
        self.betaFPYlist = set(self.FPYlist.keys()).intersection(betaSpectraDB.istplist.keys())
        self.missing_list = set(self.FPYlist.keys()) - set(betaSpectraDB.istplist.keys())

        #Iterate through every fission product in the Reactor Model.
        for FPZAI in tqdm(self.betaFPYlist, desc="Summing beta/neutrino spectra for "+str(len(self.betaFPYlist))+ " fission products", disable=silent):
            # get the yield of the fission products
            thisyield = self.FPYlist[FPZAI].y
            yielderr = self.FPYlist[FPZAI].yerr

            # for IFP calculation, adjust the decay rate of the target isotope
            # by the rate of isotope that are decayed in the time window
            if (ifp_begin < ifp_end):
                adjustment = betaSpectraDB.istplist[FPZAI].decay_time_adjust(ifp_begin, ifp_end)
                thisyield *= adjustment
                yielderr *= adjustment

            #Pull the beta spectrum and the beta uncertainties, add the product of the Beta Spectra and yield
            #To the total spectrum
            betaSpectraList[FPZAI] = betaSpectraDB.istplist[FPZAI].spectrum
            betaUncertainty[FPZAI] = betaSpectraDB.istplist[FPZAI].uncertainty
            self.spectrum += betaSpectraList[FPZAI]*thisyield

            '''
            obsolete code
            '''
            self.yieldUnc += betaSpectraList[FPZAI]*yielderr
            self.modelUnc += betaUncertainty[FPZAI]*thisyield

        # Uncertainty calculation
        #Have to make a 2D Lattice of every single combination of fission products i_j
        for i in tqdm(self.betaFPYlist, desc="Calculating uncertainties of "+str(len(self.betaFPYlist))+ " fission products"):
            for j in self.betaFPYlist:
                adjustmenti = 1
                adjustmentj = 1
                #If we're doing an IFP calculation, calculate the adjustments need to be made to
                #Both fission products.
                if (ifp_begin < ifp_end):
                    adjustmenti = betaSpectraDB.istplist[i].decay_time_adjust(ifp_begin, ifp_end)
                    adjustmentj = betaSpectraDB.istplist[j].decay_time_adjust(ifp_begin, ifp_end)

                #Pull the yields, yeild errors, the beta Spectra, and beta Spectral uncertainty
                #For both fission products (i and j)
                yi = self.FPYlist[i].y*adjustmenti
                yerri = self.FPYlist[i].yerr*adjustmenti
                fi = betaSpectraList[i]*adjustmenti
                ferri = betaUncertainty[i]*adjustmenti

                yj = self.FPYlist[j].y*adjustmentj
                # yerrj = self.FPYlist[j].yerr*adjustmentj
                fj = betaSpectraList[j]*adjustmentj
                ferrj = betaUncertainty[j]*adjustmentj

                #Carry out the uncertainty calculation
                # if covariance matrix were not loaded, make cov diagonal variance
                cov_ij = 0
                if j not in self.FPYlist[i].cov:
                    cov_ij = yerri*yerri if i == j else 0
                else:
                    cov_ij = self.FPYlist[i].cov[j]

                sigmay_ij = fi*cov_ij*fj

                self.uncertainty += sigmay_ij

        # if allowed, add beta model uncertainty to the result
        if modelunc:
            self.uncertainty += self.modelUnc**2

        self.uncertainty = np.sqrt(self.uncertainty)
# class that builds a fission reactor model with dictionary of all FPYs from all
# reactor compositions.
class FissionModel:
    '''
    Attributes
    ----------
    FPYlist : (dictionary)
        A dictionary of fission product yields. contains FissionIstps
    W : (int)
        The weight of this fission model

    Methods
    -------
    AddContribution(isotope, Ei, fraction, d_fraction=0, IFP=False):
        Add all the daughter fission products associated with the given isotope into the FPYlist
    AddIstp(Z, A, fraction, isomer=0, d_frac=0.0):
        Add a specific fission product isotope into the FPYlist
    CustomCovariance(DBname, percent = False, rel = False):
        Method to load a customized covariance matrix into the model
    GetNuclide(ZAI):
        return nuclide ZAI < (ZAI = Atomic number, Atomic mass, Isomer number)
    SaveToFile(filename):
        save the Fission Product yields into a CSV file.
    DrawBranches(figname):
        Draw the branch fractions of associated with this fission model.
    '''

    def __init__(self, W = 1.0):
        self.FPYlist = {}
        self.W = 1.0

    # method that accumulates FPYs of fission isotopes in the list of FPY
    def AddContribution(self, isotope, fraction, d_frac=0.0):
        """
           Add the FPYs of a fission isotope into the the list of FPYs in the model

            Parameters:
                isotope (FissionIstp) : The fission isotope whose products you want to add to the model.
                Ei (float) : The neutron energy that is causing the fissions to occur (0.0, 0.4/0.5, 14)
                fraction (float) : The fractional contribution that this isotope has on the overall model
                d_frac (float) : The uncertainty in the contribution
                IFP (boolean) : determines whether to include the independant fission products in the model

            Returns:
                None
        """

        #Set FPYlist to the Cumulative fission product list of the inputed isotope
        FPYs = copy.deepcopy(isotope.FPYlist)

        for FPZAI in FPYs:
            if FPYs[FPZAI].y == 0: continue #If the yield for the specific fission product is 0, skip it (ZAI)
            FPYs[FPZAI].Contribute(self.W*fraction, d_frac) # add the contribution of a fission product
            #and scale it by the fractional contribution of the parent isotope
            if FPZAI not in self.FPYlist: #Check to see if the fission product is in the FPYlist
                self.FPYlist[FPZAI] = FPYs[FPZAI] #if it is not, add the fission product to our FPYlist
            else:
                #If it is in the list, add the yield and error from the inputted fission product
                #To the fission product in the model.
                self.FPYlist[FPZAI].y += FPYs[FPZAI].y
                self.FPYlist[FPZAI].yerr += FPYs[FPZAI].yerr
                #Also add the covariance from the inputted fission product to the fission product in the model
                self.FPYlist[FPZAI].AddCovariance(FPYs[FPZAI]) # summing the covariance matrix

    # method that accumulates beta-decaying isotopes into the list of FPY
    def AddIstp(self, Z, A, fraction, isomer=0, d_frac = 0.0):
        """
           Adds a specific beta decaying isotope to the list of FPYs

            Parameters:
                Z (int) : The isotopes atomic number
                A (int) : The isotopes atomic weight
                fraction (float) : the fractional contribution this isotope has on the model
                isomer (int) : The isomeric number of this isotope
                d_frac (float) : The uncertainty in the fractional contribution
            Returns:
                None
        """

        nuclide = FPNuclide(Z*10000+A*10+isomer, fraction, d_frac) #Create a nuclide with the given
        #Isotopic information
        FPZAI = nuclide.FPZAI # Generate the FPZAI number associated with the inputted isotope
        if FPZAI not in self.FPYlist:
            self.FPYlist[FPZAI] = nuclide #If the nuclide is not in our model, add it to the model
        else:
            #If it is in the model, add the contribution and error of this isotope to the models copy of this isotope
            self.FPYlist[FPZAI].y + nuclide.y
            self.FPYlist[FPZAI].yerr + nuclide.yerr

    def CustomCovariance(self, DBname, percent = False, rel = False):
        """
           Loads a user defined covariance matrix into the model. The loading of the covariance matrix into the model
           is very similar to the FissionIstp method LoadCovarianceDB

            Parameters:
                DBname (String) : The path to the user defined covariance csv file. has the format "/path/to/file"
                percent (boolean) : Determines whether the covarainces are relative or absolute
                rel (boolean) : Determines whether the product of the yields of each pair of isotopes is a value, or 1
            Returns:
                None
        """

        # determine the whether the covariance are relative or absolute
        rate = 1e4 if percent else 1
        with open(DBname) as inputfile:
            reader = csv.DictReader(inputfile, dialect='excel', delimiter=',')

            for row in reader:
                row_id = int(row[''])
                z = int(row_id/10000)
                i = int(row_id%10)
                a = int((row_id-z*10000-i*1000)/10)
                fpzai = z*10000+a*10+i
                # keystr = str(key)
                if fpzai not in self.FPYlist:
                    continue

                for corrzai in self.FPYlist:
                    col_id = int(corrzai)
                    # if col_id is not found in the covaraince matrix, set
                    # value to zero
                    y1 = self.FPYlist[fpzai].y
                    y2 = self.FPYlist[corrzai].y
                    y_prod = (y1*y2)**rel
                    if str(col_id) in row:
                        self.FPYlist[fpzai].cov[corrzai] = y_prod*float(row[str(corrzai)])/rate
                    else:
                        self.FPYlist[fpzai].cov[corrzai] = 0.0

            # if element not found, set diagonal element to 'yerr'
            for nuclide in self.FPYlist:
                if nuclide not in self.FPYlist[nuclide].cov:
                    self.FPYlist[nuclide].cov = {otherNuclide: 0 for otherNuclide in self.FPYlist}
                    self.FPYlist[nuclide].cov[nuclide] = self.FPYlist[nuclide].yerr**2

    # Get the nuclide information based on ZAI
    def GetNuclide(self, ZAI):
        """
            Returns the nuclide information based on its' ZAI number

            Parameters:
                ZAI (int) : The ZAI number of the nuclide
            Returns:
                self.FPYlist[ZAI] (FissionNuclide) : The FissionNuclide associated with that ZAI number

        """
        return self.FPYlist[ZAI]

    # Function that save FPYs in a csv file
    def SaveToFile(self, filename):
        """
            A function that saves all FPYs into a CSV file for external use

            Parameters:
                filename (String) : The filename of the csv file you want to save
            Returns:
                None

        """
        #create a csv file and open it.
        with open(filename, 'w', newline='') as outputfile:
            colNames = ['Z', 'A', 'I', 'Y', 'Yerr'] #Create the column names
            writer = csv.DictWriter(outputfile, fieldnames=colNames)
            writer.writeheader()
            #iterate through all the nuclides in the FPYlist and populate the csv file with the associated information
            for FPZAI in self.FPYlist:
                nuclide = self.FPYlist[FPZAI]
                writer.writerow({'Z':nuclide.Z, 'A':nuclide.A, 'I':nuclide.isomer, 'Y':nuclide.y, 'Yerr':nuclide.yerr})

    # Draw a histogram of branch fractions
    def DrawBranches(self, figname):
        """
            Draw a histogram of the branch fractions for the model

            Parameters:
                figname (String) : The name of the image that will get generated
            Returns:
                None
        """
        print("Drawing branches...")
        fig, ax = plt.subplots()
        alist = np.arange(50, 180, 1)
        branchlist = np.zeros(len(alist))
        errlist = np.zeros(len(alist))

        for FPZAI in self.FPYlist:
            A = self.FPYlist[FPZAI].A
            branchlist[A-50] += self.FPYlist[FPZAI].y
            errlist[A-50] += self.FPYlist[FPZAI].yerr

        ax.errorbar(alist, branchlist, yerr=errlist)
        ax.set(xlabel='A', ylabel='fraction', title='Branch fractions')
        fig.savefig(figname)

if __name__ == "__main__":
    from conflux.BetaEngine import BetaEngine
    from conflux.SumEngine import SumEngine
    xbins = np.arange(0, 20, 0.1)

    U235 = FissionIstp(92, 235, Ei = 0.5, DB='ENDF', IFPY=True)
    U235.LoadFissionDB(Ei = 0.5)
    U235.LoadCorrelation(DB='ENDF')

    betaSpectraDB = BetaEngine(xbins=xbins)
    #betaSpectraDB = BetaEngine(newlist)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    U235.CalcBetaSpectra(betaSpectraDB)
    # Pu239.CalcBetaSpectra(betaSpectraDB)

    newsum = SumEngine(xbins = xbins)
    newsum.AddFissionIstp(U235, "U235", 1, 0)
    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    ax.set(xlabel='E (MeV)', ylabel='neutrino count')
    ax.errorbar(newsum.xbins, newsum.spectrum, yerr=newsum.uncertainty, label="test spectrum")
    ax.legend()
    plt.show()
