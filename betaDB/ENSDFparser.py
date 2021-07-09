# universal modules
import sys
from xml.dom import minidom
from os import listdir
import fortranformat as ff

# prepare the xml data file structrue
class XMLedit:
    def __init__(self):
        self.root = minidom.Document()
        self.DB = self.root.createElement('betaDB')
        self.root.appendChild(self.DB)

    def createIsotope(self, isotopeID):
        self.isotope = self.root.createElement('isotope')
        self.isotope.setAttribute('isotope', (isotopeID))
        self.DB.appendChild(self.isotope)

    def editBranch(self, fraction, end_point_E, forbideness, sigma_frac='0.0', sigma_E0='0.0'):
        branch = self.root.createElement('branch')
        branch.setAttribute('fraction', (fraction))
        branch.setAttribute('sigma_frac', (sigma_frac))
        branch.setAttribute('end_point_E', (end_point_E))
        branch.setAttribute('sigma_E0', (sigma_E0))
        branch.setAttribute('forbideness', (forbideness))
        self.isotope.appendChild(branch)

    def saveXML(self, save_path_file):
        xml_str = self.root.toprettyxml(indent ="\t")

        with open(save_path_file, "w") as f:
            f.write(xml_str)

inputfile = open("/Users/zhang39/Data/ENSDF/ensdf.137", "r")
for line in inputfile.readlines():
    A = (line[:3].strip())
    element = line[3:5]
    MT = line[5:8]

    if MT == "  B": print(int(A), element, MT)
