# universal modules
import sys
from xml.dom import minidom
from os import listdir
import csv
import numpy as np

# global method to determine if string contains float
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# global method to generate a dictionary of element and Z
def element_to_Z():
    zdict = {}
    with open("./Z_to_element.csv") as csvinput:
        csvreader = csv.reader(csvinput, delimiter=',')
        for row in csvreader:
            if row[0].isdigit:
                zdict[row[1]] = row[0]
    return zdict

# global dictionary of element and Z
elementdict = element_to_Z()

# global method to quantify simplified uncertainty indication
def transUncert(mean, unc):
    decimal = mean.find(".")
    exponent = mean.find("E")
    magnitude = len(mean[decimal:(exponent-1)]) if exponent>0 else len(mean[decimal:(exponent)])
    correction = int(mean[exponent+1:]) if exponent>0 else 0
    return float(unc)*pow(10, -magnitude+correction)

# global method to convert spin to numerical veriable
def convert_J(chars):

    chars = chars.replace('(', '')
    chars = chars.replace(')', '')
    J = 0
    pi = 1
    if not any(c.isdigit() for c in chars):
        if '-' in chars:
            pi = -1
        return [pi, J]

    it = 0
    for Js in chars.split(','):
        J = np.zeros(len(chars.split(',')))
        pi = np.zeros(len(chars.split(',')))
        num, denom = Js.split('/')
        if '-' in denom:
            denom = denom.replace('-','')
            pi[it] = (-1)
        elif '+' in denom:
            denom = denom.replace('+','')
            pi[it] = (1)
        else:
            if '-' in chars: pi[it] = -1
            else: pi[it] = 1

        J[it] = (float(num)/float(denom))
        it +=1
    return [pi, J]

# prepare the xml data file structrue
class XMLedit:
    def __init__(self):
        self.root = minidom.Document()
        self.DB = self.root.createElement('betaDB')
        self.root.appendChild(self.DB)

    # beta-decay isotopes
    def createIsotope(self, isotopeID, Q = '0.0', HL = '0.0'):
        self.isotope = self.root.createElement('isotope')
        self.isotope.setAttribute('isotope', str(isotopeID))
        self.isotope.setAttribute('Q', str(Q))
        self.isotope.setAttribute('HL', str(HL))
        self.DB.appendChild(self.isotope)

    # decay branches of the isotope
    def editBranch(self, fraction, end_point_E, forbideness, sigma_frac='0.0', sigma_E0='0.0'):
        branch = self.root.createElement('branch')
        branch.setAttribute('fraction', str(fraction))
        branch.setAttribute('sigma_frac', str(sigma_frac))
        branch.setAttribute('end_point_E', str(end_point_E))
        branch.setAttribute('sigma_E0', str(sigma_E0))
        branch.setAttribute('forbideness', str(forbideness))
        self.isotope.appendChild(branch)

    def saveXML(self, save_path_file):
        xml_str = self.root.toprettyxml(indent ="\t")

        with open(save_path_file, "w") as f:
            f.write(xml_str)

# class to save parent isotope information
class ParentIstp:
    def __init__(self, line):
        self.A = int(line[:3].strip())
        element = line[3:5].capitalize().strip()
        self.Z = int(elementdict[element])
        self.level = float(line[9:19].strip())/1e3 #convert to MeV
        self.HL = line[39:49]
        self.Emax = float(line[64:74].strip())/1e3
        self.d_Emax = transUncert(line[64:74].strip(), line[74:76].strip())/1e3
        self.pi, self.J = convert_J(line[21:39])

# class to save information of the level before decay
class DecayLevel:
    def __init__(self, line):
        self.level = float(line[9:19].strip())/1e3 if is_number(line[9:19]) else 0 #convert to MeV
        self.d_level = transUncert(line[9:19].strip(), line[19:21].strip())/1e3 if is_number(line[19:21]) else 0
        self.HL = line[39:49]
        self.pi, self.J = convert_J(line[21:39])

# class to save information of decay branches
class DecayBranch:
    def __init__(self, line):
        self.E0 = float(line[9:19].strip())/1e3 if is_number(line[9:19]) else 0
        self.sigma_E0 = transUncert(line[9:19].strip(), line[19:21].strip())/1e3 if is_number(line[19:21]) else 0
        self.fraction = float(line[21:29].strip())/100 if is_number(line[21:29]) else 1
        self.sigma_frac = transUncert(line[21:29].strip(), line[29:31].strip())/100 if is_number(line[29:31]) else 0
        self.forbidden = line[77:79].strip() if not line[77:79].isspace() else "0"

def ENSDFbeta(inputfile):
    xmloutput = XMLedit()
    lastline = ""
    betabool = False
    beginlevel = 0.
    isomer = 0
    for line in inputfile.readlines():

        if line.isspace():
            if (betabool == True): print("EMPTY LINE")
            lastline = ""
            betabool = False

        # find parent isotope of beta decay (ignore other types of decay)
        MT = line[5:8]      # Check datatype
        betatag = "B- DECAY"
        if MT == "   " and betatag in line:
            lastline = line
            betabool = True

        if betabool == True:
            # save information of the parent isotope
            if MT == "  P":
                decayparent = ParentIstp(line)
                if decayparent.level > beginlevel:
                    beginlevel = decayparent.level
                    isomer+=1
                else:
                    beginlevel = 0
                    isomer = 0
                ZAI = int(decayparent.Z*1e4 + decayparent.A*10 + isomer)
                # save the parent isotope information
                xmloutput.createIsotope(ZAI, decayparent.Emax, decayparent.HL)
                lastline = line
                print(lastline)
            elif MT == "  L":
                lastline = line
                print(lastline)
            elif MT == "  B":
                decaybranch = DecayBranch(line)
                decaydaughter = DecayLevel(lastline)
                if decaybranch.E0 == 0:
                    decaybranch.E0 = decayparent.Emax - decaydaughter.level
                    decaybranch.sigma_E0 = decayparent.d_Emax + decaydaughter.d_level
                # calculate spin difference
                decaybranch.forbidden = [decayparent.pi*decaydaughter.pi*(decaydaughter.J != 0), abs(decayparent.J - decaydaughter.J)*(decaydaughter.J != 0)]
                print(decayparent.pi,decaydaughter.pi )
                print(abs(decayparent.J - decaydaughter.J))
                print(decaybranch.forbidden)
                forbid_str = ""
                for i in range(len(decaybranch.forbidden[1])):
                    dspin_str = str(int(decaybranch.forbidden[1][i]))
                    if decaybranch.forbidden[0][i] < 0:
                        dspin_str = "-"+dspin_str
                    if i == 0:
                        forbid_str = forbid_str+dspin_str
                    else:
                        forbid_str = forbid_str+","+dspin_str
                # print(forbid_str)
                # save the decay branch information
                xmloutput.editBranch(str(decaybranch.fraction), str(decaybranch.E0), forbid_str, str(decaybranch.sigma_frac), str(decaybranch.sigma_E0))
                lastline = line
                print(lastline)

    xmloutput.saveXML("ENSDFtest.xml")

    return 0


if __name__ == "__main__":
    inputfile = open("/Users/zhang39/Data/ENSDF/ensdf.003", "r")
    ENSDFbeta(inputfile)
