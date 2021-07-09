# fission product yield engine:
# input: dictionary of fission fractions, with {'ZA', fission_fraction}
# output: dictionary of fission yield isotopes, with {'ZA', isotope_fraction}

# universal modules
import sys
import numpy as np
from os import listdir
import xml
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

class FPNuclide:
    def __init__(self, ZA, isomer, y, yerr):
        self.Z = int(ZA/1000)
        self.A = int(ZA%1000)
        self.N = self.A-self.Z
        self.isomer = isomer
        self.y = y
        self.yerr = yerr

    def Contribute(self, fraction, d_fraction=0):
        self.y *= fraction
        self.yerr = self.y*np.sqrt((self.yerr/self.y)**2 + (d_fraction/fraction**2))

class FissionIstp:
    def __init__(self, Z, A):
        self.Z = Z
        self.A = A
        self.CFPY = {}
        self.IFPY = {}

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
                    nuclidelist = []

                    for CONT in LIST:
                        FPZA = int(CONT.attrib['ZA'])
                        isomeric = float(CONT.attrib['FPS'])
                        Y = float(CONT.attrib['Y'])
                        DY = float(CONT.attrib['DY'])
                        #print(MT, Ei, FPZA, Y)
                        nuclide = FPNuclide(ZA=FPZA, isomer=isomeric, y=Y, yerr=DY)
                        nuclidelist.append(nuclide)

                    if (MT == 'CFP'):
                        self.CFPY[Ei] = nuclidelist
                    if (MT == 'IFP'):
                        self.IFPY[Ei] = nuclidelist

        assert(istpfound) # assert error if isotope not found in DB
        print('FPY list loaded.')

class FissionModel:
    def __init__(self, W = 1.0):
        self.FPYlist = {}
        self.W = 1.0

    def AddContribution(self, isotope, Ei, fraction, d_frac=0.0):
        if Ei not in isotope.CFPY:
            print('Isotope '+str(isotope.A)+' has no such fission type with Ei = '+str(Ei)+' MeV!')
            return

        for nuclide in isotope.CFPY[Ei]:
            if nuclide.y == 0: continue
            FPZAI = int(nuclide.Z*10000+nuclide.A*10+nuclide.isomer)
            nuclide.Contribute(self.W*fraction, d_frac)
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = nuclide
            else:
                print (self.FPYlist[FPZAI].y, nuclide.y)
                self.FPYlist[FPZAI].y + nuclide.y
                self.FPYlist[FPZAI].yerr + nuclide.yerr

    def AddNFIstp(self, Z, A, fraction, isomer=0, d_frac = 0.0):
        nuclide = FPNuclide(Z*1000+A, isomer, fraction, d_frac)
        FPZAI = int(nuclide.Z*10000+nuclide.A*10+nuclide.isomer)
        if FPZAI not in self.FPYlist:
            self.FPYlist[FPZAI] = nuclide
        else:
            self.FPYlist[FPZAI].y + nuclide.y
            self.FPYlist[FPZAI].yerr + nuclide.yerr

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

if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadDB()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadDB()

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=0.6, d_frac=0.05)
    model.AddContribution(isotope=Pu239, Ei = 0, fraction=0.4, d_frac=0.05)
    for FPZAI in model.FPYlist:
        print('nuclide: ', FPZAI, 'y: ', model.FPYlist[FPZAI].y, 'yerr: ', model.FPYlist[FPZAI].yerr )
    model.DrawBranches("frac.png")
