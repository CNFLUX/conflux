# fission product yield engine:
# input: dictionary of fission fractions, with {'ZA', fission_fraction}
# output: dictionary of fission yield isotopes, with {'ZA', isotope_fraction}

# universal modules
import sys
import csv
import numpy as np
from os import listdir
import xml
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

# this class saves nuclide info of fission products
class FPNuclide:
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

        # Also normalize the covariance matrix.
        for key in self.cov:
            if key == self.FPZAI:
                self.cov[key] = self.y*fraction*np.sqrt((self.cov[key]/self.y)**2 + (d_fraction/fraction**2))
            else:
                self.cov[key] *= fraction

        self.y *= fraction  # multiply fission yield with fission fraction

    # Method to add covariance matrices together
    def AddCovariance(self, newNuclide):
        print("ADDED COV MAT")
        for key in newNuclide.cov:
            if key not in self.cov:
                self.cov[key] = newNuclide[key]
            else:
                self.cov[key] += newNuclide[key]
            #if key == self.FPZAI: print(key, self.cov[key])

# Class that counts fission products of a specified fission isotope
class FissionIstp:
    def __init__(self, Z, A):
        self.Z = Z      # fission isotope Z
        self.A = A      # fission isotope A
        self.CFPY = {}  # dictionary of cumulative fission yields {"FPZAI", FPNuclide}
        self.IFPY = {}  # dictionary of independent fission yields {"FPZAI", FPNuclide}

    # method that load xml database of FPY and save nuclide info in dictionaries.
    def LoadDB(self, DBpath = './fissionDB/'):
        print('Reading FPY DB from folder: '+DBpath+'...')
        fileList = listdir(DBpath)

        istpfound = False
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'xml':
                continue
            if (str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            istpfound = True
            break
        assert(istpfound) # assert error if isotope not found in DB

        DBname = DBpath+filename
        print('Reading FPY DB: '+DBpath+filename+'...')
        tree = ET.parse(DBname)
        root = tree.getroot()

        for HEAD in root:
            MT = HEAD.attrib['MT']

            for LIST in HEAD:
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

                if (MT == 'CFP'):
                    self.CFPY[Ei] = nuclidelist
                if (MT == 'IFP'):
                    self.IFPY[Ei] = nuclidelist


        print('FPY list loaded.')

    # Method to read the prepackaged covariance csv file
    def LoadCovariance(self, DBpath = './fissionDB/'):
        print("Reading covariance matrix: "+DBpath+"...")
        fileList = listdir(DBpath)
        istpfound = False
        e_neutron = {'Thermal': 0, 'Fast': 0.5, 'High': 14}
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'csv':
                continue
            if ("cov" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            istpfound = True

            Ei = 0
            for mode in e_neutron:
                if mode in namecache[0]:
                    Ei = e_neutron[mode]

            DBname = DBpath+filename
            print('Reading FPY DB: '+DBpath+filename+'...')
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
                    for key in row:
                        if key == '':
                            continue
                        else:
                            col_id = int(key)
                            z = int(col_id/10000)
                            i = int((col_id-z*10000)/1000)
                            a = int(col_id-z*10000-i*1000)
                            corrzai = z*10000+a*10+i
                            self.CFPY[Ei][fpzai].cov[corrzai] = float(row[key])

                # if element not found, set diagonal element to 'yerr'
                for nuclide in self.CFPY[Ei]:
                    if nuclide not in self.CFPY[Ei][nuclide].cov:
                        self.CFPY[Ei][nuclide].cov = {key: 0 for key in row}
                        self.CFPY[Ei][nuclide].cov[nuclide] = self.CFPY[Ei][nuclide].yerr

        assert(istpfound) # assert error if isotope not found in DB

    # Method to read the prepackaged correlation csv file
    def LoadCorrelation(self, DBpath = './fissionDB/'):
        print("Reading correlation matrix from: "+DBpath+"...")
        fileList = listdir(DBpath)
        istpfound = False
        e_neutron = {'Thermal': 0, 'Fast': 0.5, 'High': 14}
        for filename in fileList:
            namecache = filename.split('.')
            if namecache[-1] != 'csv':
                continue
            if ("corr" not in namecache[0] or str(self.Z) not in namecache[0] or str(self.A) not in namecache[0]):
                continue
            istpfound = True

            Ei = 0
            for mode in e_neutron:
                if mode in namecache[0]:
                    Ei = e_neutron[mode]

            DBname = DBpath+filename
            print('Reading FPY DB: '+DBpath+filename+'...')
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
                    for key in row:
                        if key == '':
                            continue
                        else:
                            col_id = int(key)
                            z = int(col_id/10000)
                            i = int((col_id-z*10000)/1000)
                            a = int(col_id-z*10000-i*1000)
                            colzai = z*10000+a*10+i
                            self.CFPY[Ei][fpzai].corr[colzai] = float(row[key])

                # if element not found, set diagonal element to 1.0
                for nuclide in self.CFPY[Ei]:
                    if nuclide not in self.CFPY[Ei][nuclide].cov:
                        self.CFPY[Ei][nuclide].corr = {key: 0 for key in row}
                        self.CFPY[Ei][nuclide].corr[nuclide] = self.CFPY[Ei][nuclide].yerr

        assert(istpfound) # assert error if isotope not found in DB

# class that builds a fission reactor model with dictionary of all FPYs from all reactor compositions.
class FissionModel:
    def __init__(self, W = 1.0):
        self.FPYlist = {}
        self.W = 1.0

    # method that accumulates FPYs of fission isotopes in the list of FPY
    def AddContribution(self, isotope, Ei, fraction, d_frac=0.0):
        assert Ei in isotope.CFPY, 'Isotope '+str(isotope.A)+' has no such fission type with Ei = '+str(Ei)+' MeV!'
        CFPYlist = isotope.CFPY[Ei]
        for FPZAI in CFPYlist:
            if CFPYlist[FPZAI].y == 0: continue
            CFPYlist[FPZAI].Contribute(self.W*fraction, d_frac) # add the contribution of a fission product, also modify uncertainties
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = CFPYlist[FPZAI]
            else:
                self.FPYlist[FPZAI].y += CFPYlist[FPZAI].y
                self.FPYlist[FPZAI].yerr += CFPYlist[FPZAI].yerr
                self.FPYlist[FPZAI].AddCovariance(CFPYlist[FPZAI]) # summing the covariance matrix

    # method that accumulates beta-decaying isotopes into the list of FPY
    def AddIstp(self, Z, A, fraction, isomer=0, d_frac = 0.0):
        nuclide = FPNuclide(Z*10000+A*10+isomer, fraction, d_frac)
        FPZAI = nuclide.FPZAI
        if FPZAI not in self.FPYlist:
            self.FPYlist[FPZAI] = nuclide
        else:
            self.FPYlist[FPZAI].y + nuclide.y
            self.FPYlist[FPZAI].yerr + nuclide.yerr

    # Get the nuclide information based on ZAI
    def GetNuclide(self, ZAI):
        return self.FPYlist[ZAI]

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

# test on how to define fission isotopes and add them to the reactor model
if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadDB()
    U235.LoadCovariance()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadDB()

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=1, d_frac=0.0)
    #model.AddContribution(isotope=Pu239, Ei = 0, fraction=0.4, d_frac=0.05)
    #model.AddIstp(390960)
    # for FPZAI in model.FPYlist:
    #     print('nuclide: ', FPZAI, 'y: ', model.FPYlist[FPZAI].y, 'yerr: ', model.FPYlist[FPZAI].yerr )
    model.DrawBranches("frac.png")
