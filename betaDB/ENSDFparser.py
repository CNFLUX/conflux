# universal modules
import sys
from xml.dom import minidom
from os import listdir
import fortranformat as ff
import csv

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
        self.level = float(line[9:19].strip())
        self.HL = line[39:49]
        self.Emax = float(line[64:74].strip())
        self.d_Emax = transUncert(line[64:74].strip(), line[74:76].strip())

# class to save information of the level before decay
class DecayLevel:
    def __init__(self, line):
        self.level = float(line[9:19].strip()) if is_number(line[9:19]) else 0
        self.d_level = transUncert(line[9:19].strip(), line[19:21].strip()) if is_number(line[19:21]) else 0
        self.HL = line[39:49]

# class to save information of decay branches
class DecayBranch:
    def __init__(self, line):
        self.E0 = float(line[9:19].strip()) if is_number(line[9:19]) else 0
        self.sigma_E0 = transUncert(line[9:19].strip(), line[19:21].strip()) if is_number(line[19:21]) else 0
        self.fraction = float(line[21:29].strip()) if is_number(line[21:29]) else 1
        self.sigma_frac = transUncert(line[21:29].strip(), line[29:31].strip()) if is_number(line[29:31]) else 0
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
                if decaybranch.E0 == 0:
                    decaydaughter = DecayLevel(lastline)
                    decaybranch.E0 = decayparent.Emax - decaydaughter.level
                    decaybranch.sigma_E0 = decayparent.d_Emax + decaydaughter.d_level
                # save the decay branch information
                xmloutput.editBranch(str(decaybranch.fraction), str(decaybranch.E0), str(decaybranch.forbidden), str(decaybranch.sigma_frac), str(decaybranch.sigma_E0))
                lastline = line
                print(lastline)

    xmloutput.saveXML("ENSDFtest.xml")

    return 0


if __name__ == "__main__":
    inputfile = open("/Users/zhang39/Data/ENSDF/ensdf.137", "r")
    ENSDFbeta(inputfile)
