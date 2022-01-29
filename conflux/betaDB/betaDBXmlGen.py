import sys
import csv
from xml.dom import minidom

# prepare the xml data file structrue
class XMLedit:
    def __init__(self):
        self.root = minidom.Document()
        self.DB = self.root.createElement('betaDB')
        self.root.appendChild(self.DB)

    # beta-decay isotopes
    def createIsotope(self, isotopeID, Q = '0.0', HL = '0.0'):
        self.isotope = self.root.createElement('isotope')
        self.isotope.setAttribute('isotope', (isotopeID))
        self.DB.appendChild(self.isotope)

    # decay branches of the isotope
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

# reading the betaDB supplied by P. Huber
def betaDBReader(filename):
    with open(filename) as inputfile:
        reader = csv.reader(filter(lambda row: row[0]!='#', inputfile), dialect='excel', delimiter='\t')
        xmloutput = XMLedit()
        prev_isotopeID = 0
        for row in reader:
            if not row:
                continue

            isotopeID = row[0]
            fraction = row[1]
            end_point_E = row[2]
            forbideness = row[3]

            if (isotopeID!=prev_isotopeID):
                xmloutput.createIsotope(isotopeID)
                xmloutput.editBranch(fraction, end_point_E, forbideness)
                prev_isotopeID = isotopeID

            else:
                xmloutput.editBranch(fraction, end_point_E, forbideness)

	xmloutput.saveXML("betaDB.xml")

    return 0

if __name__ == "__main__":
    filename = sys.argv[1]
    betaDBReader(filename)
