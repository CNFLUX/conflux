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
class FissionIstp:
    """
    Class to handle all the Fission Nuclide Information.

    ...

    Attributes
    ----------
    Z : (int)
        Atomic number of your nuclide
    A : (int)
        Atomic mass of your nuclide
    CFPY : (dictionary)
        Dictionary of cumulative fission yields {"FPZAI", FPNuclide}
    IFPY : (dictionary)
        Dictionary of independent fission yields {"FPZAI", FPNuclide}
    DBTitle : (dictionary)
        Dictionary of the two fission databases included in CONFLUX (can be downloaded using the file in the aux folder)
    
    Methods
    -------
    LoadFissionDB(customDB = None, defaultDB='ENDF'):
        Function to load data from a specific fission database. customDB has format /path/to/DB. Defaults to ENDF
    LoadCovariance(customDB = None, defaultDB='ENDF', percent=True):
        Load the covariance matrices associated with the fission isotopes. customDB has format /path/to/DB. Defaults to ENDF
    """

    def __init__(self, Z, A):
        self.Z = Z      
        self.A = A      
        self.CFPY = {}  
        self.IFPY = {} 
        self.DBtitle = {'ENDF':'nfy', 'JEFF':'nfpy'}

    # method that load xml database of FPY and save nuclide info in dictionaries.
    def LoadFissionDB(self, customDB = None, defaultDB='ENDF'):
        DBname = customDB
        if DBname == None or not os.path.exists(DBname): #Check if the user gave a valid Database path
            DBpath = os.environ["CONFLUX_DB"]+"/fissionDB/"+defaultDB+"/"
            if DBname != None:
                print('Custom DB: '+ DBname + ' NOT found!')
            print('Reading default FPY DB from folder: '+DBpath+'...')
            fileList = listdir(DBpath) #Get the list of files in the Database directory
            istpfound = False
            for filename in fileList: #iterate through the list of files in the directory
                namecache = filename.split('.')
                if namecache[-1] != 'xml': #if the file is not an xml file, continue
                    continue
                if (self.DBtitle[defaultDB] not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                    continue #Alternatively, if the isotope is not in the list, continue
                istpfound = True #Else, assert that the isotope is in the DB
                break
            assert(istpfound) # assert error if isotope not found in DB

            DBname = DBpath+filename #this is the isotope that we found in the list.
            assert(DBname)

        print('Reading FPY DB: '+DBname+'...')
        tree = ET.parse(DBname)
        root = tree.getroot()

        for HEAD in root:
            MT = HEAD.attrib['MT'] #I'm assuming this is measurement type, whether it is 
            #Individual fission products or cumulative fission products.

            for LIST in HEAD: #work through this particular fission isotope, pull
                #Out the relevant branch information. 
                Ei = float(LIST.attrib['Ei'])
                Ei /= 1e6
                if (Ei < 0.01): Ei = 0.
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

                #Sort out this isotope and all of its' branches into either 
                #Cumulative or Independant yields. 
                if (MT == 'CFP'):
                    self.CFPY[Ei] = nuclidelist
                if (MT == 'IFP'):
                    self.IFPY[Ei] = nuclidelist

                print('FPY list loaded. '+ str(len(nuclidelist))+" fission products were found with Ei of "+str(Ei)+" MeV ")

    # Method to read the prepackaged covariance csv file
    # This function has to be called after loading the fission DB for neutrino
    # flux calcuation.
    def LoadCovariance(self, customDB = None, defaultDB='ENDF', percent=True):
        DBpath = customDB
        #Figure out where the DB and files are in a similar method we loaded the 
        #Fission Database.
        if DBpath == None or not os.path.exists(DBpath):
            DBpath = os.environ["CONFLUX_DB"]+"/fissionDB/"+defaultDB+"/"
            print("Reading covariance matrices in: "+DBpath+"...")
        fileList = listdir(DBpath)
        assert(DBpath)
        istpfound = False
        #By default, set Neutron energy to 0 (thermal), 0.5 (fast), or 14 (relativistic) << default DB is ENDF
        e_neutron = {'T': 0, 'F': 0.5, 'H': 14}
        #If JEFF is chosen, adjust the neutron energies. 
        if defaultDB=='JEFF':
            e_neutron = {'T': 0, 'F': 0.4, 'H': 14}
        filesfound = []
        #Check for database files, and append them into a files list.
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'csv':
                continue
            if ("normed_cov" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            filesfound.append(filename)

        # determine the whether the covariance are relative or absolute
        rate = 1e4 if percent else 1
        Ei = 0
        for mode in e_neutron:
            for filename in filesfound:
                if mode in filename:
                    Ei = e_neutron[mode]
                    #Check to see if the specific neutron energy
                    #Is in the file we are looking at
                else:
                    continue

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
                        if fpzai not in self.CFPY[Ei]:
                            continue 
                        #If the isotope is not in the Cumulative fission dictionary    
                        #At the given energy, skip
                        for corrzai in self.CFPY[Ei]:
                            col_id = int(corrzai)
                            z = int(col_id/10000)
                            a = int((col_id-z*10000)/10)
                            i = int(col_id-z*10000-a*10)
                            key = z*10000+a+i*1000
                            keystr = ' '+str(key)
                            # if key is not found in the covaraince matrix, set
                            # value to zero

                            if keystr in row:
                                self.CFPY[Ei][fpzai].cov[corrzai] = float(row[keystr])/rate
                            else:
                                self.CFPY[Ei][fpzai].cov[corrzai] = 0.0

                    # if element not found, set diagonal element to 'yerr'
                    for nuclide in self.CFPY[Ei]:

                        if nuclide not in self.CFPY[Ei][nuclide].cov:
                            self.CFPY[Ei][nuclide].cov = {otherNuclide: 0 for otherNuclide in self.CFPY[Ei]}
                            self.CFPY[Ei][nuclide].cov[nuclide] = self.CFPY[Ei][nuclide].yerr**2

        assert(filesfound) # assert error if isotope not found in DB

    # Method to read the prepackaged correlation csv file
    # This function has to be called after loading the fission DB for neutrino
    # flux calcuation.

    #Most of the comments will be the same for both this and for the LoadCovariance function
    def LoadCorrelation(self, customDB = None, defaultDB='ENDF'):
        DBpath = customDB  #Figure out where the DB and files are in a similar method we loaded the 
        #Fission Database. 
        if DBpath == None or not os.path.exists(DBpath):
            DBpath = os.environ["CONFLUX_DB"]+"/fissionDB/"+defaultDB+"/"
            print("Reading correlation matrices in: "+DBpath+"...")
        fileList = listdir(DBpath)
        assert(DBpath)
        istpfound = False
        #By default, set Neutron energy to 0 (thermal), 0.5 (fast), or 14 (relativistic) << default DB is ENDF
        e_neutron = {'T': 0, 'F': 0.5, 'H': 14}
        if defaultDB=='JEFF':  #If JEFF is chosen, adjust the neutron energies. 
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

        Ei = 0
        for mode in e_neutron:
            for filename in filesfound:
                if mode in filename:
                    Ei = e_neutron[mode]
                     #Check to see if the specific neutron energy
                    #Is in the file we are looking at
                else:
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
                        if fpzai not in self.CFPY[Ei]:
                            continue #If the isotope is not in the Cumulative fission dictionary    
                            #At the given energy, skip
                        for corrzai in self.CFPY[Ei]:
                            
                            col_id = int(corrzai)
                            z = int(col_id/10000)
                            a = int((col_id-z*10000)/10)
                            i = int(col_id-z*10000-a*10)
                            key = z*10000+a+i*1000
                            keystr = ' '+str(key)
                            # if key is not found in the correlation matrix, set
                            # value to zero
                            if keystr in row:
                                self.CFPY[Ei][fpzai].corr[corrzai] = float(row[keystr])
                            else:
                                self.CFPY[Ei][fpzai].corr[corrzai] = 0.0

                    # if element not found, set diagonal element to 1.0
                    for nuclide in self.CFPY[Ei]:
                        if nuclide not in self.CFPY[Ei][nuclide].corr:
                            self.CFPY[Ei][nuclide].corr = {otherNuclide: 0 for otherNuclide in self.CFPY[Ei]}
                            self.CFPY[Ei][nuclide].corr[nuclide] = 1.0

                self.CalcCovariance(Ei) #Once you've loaded the correlation database, calculate the
                    #subsequent covariance matrix from those correlations. 
        assert(filesfound) # assert error if isotope not found in DB

    # Method to calculate a covariance matrix from the correlation information
    # of each fission product at a given neutron energy. 
    def CalcCovariance(self, Ei):
        #Iterate over all fission nuclides inside the Cumulative fission dictionary
        for i in self.CFPY[Ei]:
            yerri = self.CFPY[Ei][i].yerr #Load the yield error of the nuclide i
            corr = self.CFPY[Ei][i].corr #Load the correlation of the nuclide i
            for j in self.CFPY[Ei]:
                yerrj = self.CFPY[Ei][j].yerr #Load the yield error of nuclide j
                self.CFPY[Ei][i].cov[j]=yerri*corr[j]*yerrj #carry out the covariance calculation 
                #(error of i) * (correlation of j) * (error of j)


# class that builds a fission reactor model with dictionary of all FPYs from all
# reactor compositions.
class FissionModel:
    def __init__(self, W = 1.0):
        self.FPYlist = {}
        self.W = 1.0

    # method that accumulates FPYs of fission isotopes in the list of FPY
    def AddContribution(self, isotope, Ei, fraction, d_frac=0.0, IFP=False):
        assert Ei in isotope.CFPY, 'Isotope '+str(isotope.A)+' has no such fission type with Ei = '+str(Ei)+' MeV!'
        FPYlist = copy.deepcopy(isotope.CFPY[Ei])
        if IFP == True:
            FPYLIST = copy.deepcopy(isotope.IFPY[Ei])
        for FPZAI in FPYlist:
            if FPYlist[FPZAI].y == 0: continue
            FPYlist[FPZAI].Contribute(self.W*fraction, d_frac) # add the contribution of a fission product, also modify uncertainties
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = FPYlist[FPZAI]
            else:
                self.FPYlist[FPZAI].y += FPYlist[FPZAI].y
                self.FPYlist[FPZAI].yerr += FPYlist[FPZAI].yerr
                self.FPYlist[FPZAI].AddCovariance(FPYlist[FPZAI]) # summing the covariance matrix

    # method that accumulates beta-decaying isotopes into the list of FPY
    def AddIstp(self, Z, A, fraction, isomer=0, d_frac = 0.0):
        nuclide = FPNuclide(Z*10000+A*10+isomer, fraction, d_frac)
        FPZAI = nuclide.FPZAI
        if FPZAI not in self.FPYlist:
            self.FPYlist[FPZAI] = nuclide
        else:
            self.FPYlist[FPZAI].y + nuclide.y
            self.FPYlist[FPZAI].yerr + nuclide.yerr

    def CustomCovariance(self, DBname, percent = False, rel = False):
        """
        Method to load a customized covariance matrix
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
        return self.FPYlist[ZAI]

    # Function that save FPYs in a csv file
    def SaveToFile(self, filename):
        with open(filename, 'w', newline='') as outputfile:
            colNames = ['Z', 'A', 'I', 'Y', 'Yerr']
            writer = csv.DictWriter(outputfile, fieldnames=colNames)
            writer.writeheader()
            for FPZAI in self.FPYlist:
                nuclide = self.FPYlist[FPZAI]
                writer.writerow({'Z':nuclide.Z, 'A':nuclide.A, 'I':nuclide.isomer, 'Y':nuclide.y, 'Yerr':nuclide.yerr})

    # Draw a histogram of branch fractions
    def DrawBranches(self, figname):
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

# U235 = FissionIstp(92,235)
# U235.LoadFissionDB()
# U235.LoadCorrelation()
#
# # U238 = FissionIstp(92,238)
# # U238.LoadFissionDB()
# # U238.LoadCorrelation(Ei = 0.5)
# # for nuclide in U235.CFPY[0.]:
# #     print(nuclide, U235.CFPY[0.][nuclide].cov)
#
# model = FissionModel()
# model.AddContribution(isotope=U235, Ei=0., fraction=1)
