#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 10:51:51 2021

@author: zhang39
"""

# universal modules
import sys
from xml.dom import minidom
import fortranformat as ff
from os import listdir, environ
import csv
import pkg_resources

# global method to generate a dictionary of element and Z
def Z_to_element():
    zdict = {}
    listname = environ['CONFLUX_DB']+'/betaDB/Z_to_element.csv'
    with open(listname) as csvinput:
        csvreader = csv.DictReader(csvinput, dialect='excel', delimiter=',')
        for row in csvreader:
            Z = int(row['Z'])
            zdict[Z] = row['Element']
            
    return zdict
    
# global dictionary of element and Z
elementdict = Z_to_element()
    
# In[0]: Prepare the class to write the output xml file
class XMLedit:
    def __init__(self, DBname):
        self.root = minidom.Document()
        self.DB = self.root.createElement(DBname)
        self.root.appendChild(self.DB)
        self.outputName = DBname

    # HEAD of the fission products data, contains:
    # ZA - indicate Z, A number of the fission isotope (ZA = Z*1000+A)
    # AWR - the isotope's multiple of neutron mass;
    # LISO - the energy state of the isotope
    #   (mostly thermal, fast and 14MeV)
    # MT - the ENDF data branch
    #   (454 for independent productions, 459 for cumulative productions)
    def createIsotope(self, istpName, isotopeID, Q = '0.0', HL = '0.0'):
        self.isotope = self.root.createElement('isotope')
        self.isotope.setAttribute('name', str(istpName))
        self.isotope.setAttribute('isotope', str(isotopeID))
        self.isotope.setAttribute('Q', str(Q))
        self.isotope.setAttribute('HL', str(HL))
        self.DB.appendChild(self.isotope)

    # LIST is the title of FPY list that corresponds to:
    # Ei - the fission triggering energy (mostly thermal, fast and 14MeV)
    # Ii - the interpolation scheme
    # NFPi - the total amount of fission products
    def editBranch(self, fraction, end_point_E, dJpi, sigma_frac='0.0',
                    sigma_E0='0.0'):
        branch = self.root.createElement('branch')
        branch.setAttribute('fraction', str(fraction))
        branch.setAttribute('sigma_frac', str(sigma_frac))
        branch.setAttribute('end_point_E', str(end_point_E))
        branch.setAttribute('sigma_E0', str(sigma_E0))
        branch.setAttribute('dJpi', str(dJpi))
        self.isotope.appendChild(branch)

    def saveXML(self, save_path_file):
        xml_str = self.root.toprettyxml(indent ="\t")

        with open(save_path_file, "w") as f:
            f.write(xml_str)

# In[1]: the method that reads ENDF fission production yield (FPY) data base
def ENDF6BetaDBreader(filename, xmloutput):
    record_reader = ff.FortranRecordReader('(A66,I4,I2,I3,I5)')
    DBname = filename.split('/')
    DBname_short = DBname[-1].split('.')
    

    with open(filename) as inputfile:

        linenum = 0
        ZA = -1
        name = ''
        Ei = []
        Ii = []
        NNi = []
        NFPi = []
        Q = 0
        NC = 0
        NDK = 0
        RTYP = -1

        datacache = []
        i = 0
        CONT = []
        for line in inputfile:

            TEXT, MAT, MF, MT, NS = record_reader.read(line)
            # independent fission products
            if (MT == 457):
                data_reader = ff.FortranRecordReader('(6E11.0)')
                data = data_reader.read(TEXT)
                # print(data, linenum)

                # reading the HEAD of independent fission products
                datacache += (data)
                
        ZA = int(datacache[0])
        AWR = float(datacache[1])
        LIS = int(datacache[2])
        LISO = int(datacache[3])    # isomer
        NST = int(datacache[4])
        NSP = int(datacache[5])
        Z = int(ZA/1000)
        A = int(ZA%1000)
        ZAI = ZA*10+LISO
        name = elementdict[Z]+str(A)
        
        HL = float(datacache[6])
        dHL = float(datacache[7])
        NC = int(datacache[10])
        E_ = []     #average decay product energy
        dE_ = []    #average decay product energy
                
        for i in range(12, 12+NC):
            newdata = datacache[i]
            if i%2 == 0:
                E_.append(newdata)
            else:
                dE_.append(newdata)
        
        newbegin = 12+NC
        spin = float(datacache[newbegin])
        parity = int(datacache[newbegin+1])
        NDK = int(datacache[newbegin+5])
        
        newbegin += 6
        decaymode = dict.fromkeys(['RTYP','RFS','Q','dQ','BR', 'dBR',])
        for i in range(NDK):
            decaymode['Rtype']=float(datacache[newbegin+i*newbegin])
            if decaymode['Rtype'] != None and int(decaymode['Rtype']%1)==0 and int(decaymode['Rtype']%1)>0:
                decaymode['Q']=float(datacache[newbegin+i*newbegin+2])
                decaymode['dQ']=float(datacache[newbegin+i*newbegin+3])
                break
        
        if not decaymode['Q']:
            return
            
        Q = decaymode['Q']/1e6
        dQ = decaymode['dQ']/1e6

        # print(name, ZAI, Q, HL)
        
        xmloutput.createIsotope(name, isotopeID=str(ZAI), Q=str(Q), HL=str(HL))
        
        newbegin += NDK*6
        
        currentline = int(newbegin/6)
        linelim = int(len(datacache)/6)
        
        for i in range(currentline, linelim):
            if float(datacache[i*6]) == 0 and float(datacache[i*6+1]) == 1:
                print(datacache[i*6:i*6+6])
                branches = int(datacache[i*6+5])
                print(datacache[i*6+6:i*6+12])
                newbegin = i*6+12
                for j in range(branches):
                    # print(datacache[newbegin+j*12:newbegin+j*12+12])
                    end_point_E = float(datacache[newbegin+j*12])/1e6
                    sigma_E0 = float(datacache[newbegin+j*12+1])1e6
                    NT = datacache[newbegin+j*12+4]
                    if NT != 6:
                        print("NT", NT)
                    forbidden = int(datacache[newbegin+j*12+7])
                    fraction = float(datacache[newbegin+j*12+8])
                    sigma_frac = float(datacache[newbegin+j*12+9])
                    
                    xmloutput.editBranch(str(fraction), str(end_point_E),
                                        str(forbidden), sigma_frac=str(sigma_frac),
                                        sigma_E0=str(sigma_frac))
                                            
    xmloutput.saveXML("ENDF_betaDB.xml")



if __name__ == "__main__":
    dirName = sys.argv[1]
    fileList = listdir(dirName)
    #print(fileList)
    xmloutput = XMLedit(dirName)
    for filename in fileList:
        if filename.split('.')[-1] == 'endf':
            ENDF6BetaDBreader(dirName+filename, xmloutput)
