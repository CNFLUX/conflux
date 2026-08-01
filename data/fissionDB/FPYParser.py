#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified Fission Product Yield (FPY) Parser for ENDF-6 format
Works for both ENDF and JEFF databases

Created on Sat Apr 17 10:51:51 2021
Modified: 2026-07-29 - Unified ENDF and JEFF parsers

@author: Xianyi Zhang, LLNL
"""

# universal modules
import sys
from xml.dom import minidom
import fortranformat as ff
from os import listdir
from pathlib import Path

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
    # LE - the number of fission triggering energies
    #   (mostly thermal, fast and 14MeV)
    # MT - the ENDF data branch
    #   (454 for independent productions, 459 for cumulative productions)
    def createHEAD(self, ZA, AWR, LE, MT):
        self.HEAD = self.root.createElement('HEAD')
        self.HEAD.setAttribute('FissionZA', str(ZA))
        self.HEAD.setAttribute('AWR', str(AWR))
        self.HEAD.setAttribute('LE', str(LE))
        self.HEAD.setAttribute('MT', str(MT))
        self.DB.appendChild(self.HEAD)

    # LIST is the title of FPY list that corresponds to:
    # Ei - the fission triggering energy (mostly thermal, fast and 14MeV)
    # Ii - the interpolation scheme
    # NFPi - the total amount of fission products
    def createLIST(self, Ei, I, NFP):
        self.LIST = self.root.createElement('LIST')
        self.LIST.setAttribute('Ei', str(Ei))
        self.LIST.setAttribute('Ii', str(I))
        self.LIST.setAttribute('NFPi', str(NFP))
        self.HEAD.appendChild(self.LIST)

    # CONT is the contents of FPY list under different energy, where
    # ZA - Z, A number of the fission product
    # FPS - isomeric state designato
    # Y - relative fission yield
    # DY - uncertainty of fission yield
    def editCONT(self, data):
        CONT = self.root.createElement('CONT')
        CONT.setAttribute('ZA', str(data[0]))
        CONT.setAttribute('FPS', str(data[1]))
        CONT.setAttribute('Y', str(data[2]))
        CONT.setAttribute('DY', str(data[3]))
        self.LIST.appendChild(CONT)

    def saveXML(self):
        xml_str = self.root.toprettyxml(indent ="\t")
        with open(self.outputName+".xml", "w") as f:
            f.write(xml_str)

# In[1]: the method that reads ENDF fission production yield (FPY) data base
def ENDF6FPYreader(filename):
    """
    Parse ENDF-6 format fission product yield file

    Args:
        filename: Path to ENDF-6 format file (.endf or .dat)

    Returns:
        XMLedit object with parsed data
    """
    record_reader = ff.FortranRecordReader('(A66,I4,I2,I3,I5)')
    DBname = filename.split('/')
    DBname_short = DBname[-1].split('.')
    xmloutput = XMLedit(DBname_short[0])

    with open(filename) as inputfile:

        linenum = 0
        ZA = -1
        LEplus1 = -1
        Ei = []
        Ii = []
        NNi = []
        NFPi = []

        i = 0
        CONT = []
        for line in inputfile:

            TEXT, MAT, MF, MT, NS = record_reader.read(line)
            # independent fission products
            if (MT == 454):
                data_reader = ff.FortranRecordReader('(6E11.0)')
                data = data_reader.read(TEXT)

                # reading the HEAD of independent fission products
                if (linenum == 0):
                    print(data)
                    ZA = int(data[0])
                    AWR = float(data[1])
                    LEplus1 = int(data[2])
                    Ei = [0]*LEplus1
                    Ii = [0]*LEplus1
                    NNi = [0]*LEplus1
                    NFPi = [0]*LEplus1
                    xmloutput.createHEAD(ZA, AWR, LEplus1-1, 454)
                    linenum += 1

                else:

                    # reading CONT []
                    # Add bounds check to prevent index errors
                    if (i < len(NFPi) and len(CONT) <= NFPi[i] and  NFPi[i]>0):
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                            datacache = [int(data[4]), float(data[5]), 0, 0]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        else:
                            if len(CONT[-1]) != 4:
                                CONT[-1][2] = float(data[0])
                                CONT[-1][3] = float(data[1])
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                            else:
                                if len(CONT[-1]) != 4:
                                    CONT[-1][2] = float(data[2])
                                    CONT[-1][3] = float(data[3])
                                    xmloutput.editCONT(CONT[-1])
                                datacache = [int(data[4]), float(data[5]), 0, 0]
                                if datacache[0] != 0:
                                    CONT.append(datacache)
                                else:
                                    if len(CONT[-1]) != 4:
                                        CONT[-1][2] = float(data[4])
                                        CONT[-1][3] = float(data[5])
                                        xmloutput.editCONT(CONT[-1])
                        if (len(CONT[-1]) != 4): continue

                    # Add bounds check here too
                    if (i < len(NFPi) and len(CONT) == NFPi[i] > 0 ):
                        print(len(CONT), CONT[-1])
                        CONT=[]
                        i+=1
                        continue

                    # reading the LIST []
                    if (i < LEplus1 and i < len(NFPi) and len(CONT) == NFPi[i] == 0):
                        Ei[i] = (data[0])
                        Ii[i] = int(data[2])
                        NNi[i] = int(data[4])
                        NFPi[i] = int(data[5])
                        xmloutput.createLIST(Ei[i], Ii[i], NFPi[i])
                        print(Ei[i], Ii[i], NNi[i], NFPi[i])
                        continue

            # cumulative fission products
            if (MT == 459):
                data_reader = ff.FortranRecordReader('(6E11.0)')
                data = data_reader.read(TEXT)

                # reading the HEAD of cumulative fission products
                if (linenum == 0):
                    print(data)
                    ZA = int(data[0])
                    AWR = float(data[1])
                    LEplus1 = int(data[2])
                    Ei = [0]*LEplus1
                    Ii = [0]*LEplus1
                    NNi = [0]*LEplus1
                    NFPi = [0]*LEplus1
                    xmloutput.createHEAD(ZA, AWR, LEplus1-1, 459)
                    linenum += 1

                else:
                    # reading CONT
                    # Add bounds check to prevent index errors
                    if (i < len(NFPi) and len(CONT) <= NFPi[i] and  NFPi[i]>0):
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                            datacache = [int(data[4]), float(data[5]), 0, 0]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        else:
                            if len(CONT[-1]) != 4:
                                CONT[-1][2] = float(data[0])
                                CONT[-1][3] = float(data[1])
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                            else:
                                if len(CONT[-1]) != 4:
                                    CONT[-1][2] = float(data[2])
                                    CONT[-1][3] = float(data[3])
                                    xmloutput.editCONT(CONT[-1])
                                datacache = [int(data[4]), float(data[5]), 0, 0]
                                if datacache[0] != 0:
                                    CONT.append(datacache)
                                else:
                                    if len(CONT[-1]) != 4:
                                        CONT[-1][2] = float(data[4])
                                        CONT[-1][3] = float(data[5])
                                        xmloutput.editCONT(CONT[-1])
                        if (len(CONT[-1]) != 4): continue

                    # Add bounds check here too
                    if (i < len(NFPi) and len(CONT) == NFPi[i] > 0 ):
                        print(len(CONT), CONT[-1])
                        CONT=[]
                        i+=1
                        continue

                    # reading the LIST
                    if (i < LEplus1 and i < len(NFPi) and len(CONT) == NFPi[i] == 0):
                        Ei[i] = (data[0])
                        Ii[i] = int(data[2])
                        NNi[i] = int(data[4])
                        NFPi[i] = int(data[5])
                        xmloutput.createLIST(Ei[i], Ii[i], NFPi[i])
                        print(Ei[i], Ii[i], NNi[i], NFPi[i])
                        continue
    xmloutput.saveXML()
    return xmloutput

def parse_directory(dirName, file_extensions=None):
    """
    Parse all FPY files in a directory

    Args:
        dirName: Directory containing FPY files
        file_extensions: List of file extensions to parse (default: ['.endf', '.dat'])
                        Use ['.endf'] for ENDF, ['.dat'] for JEFF, or both

    Returns:
        None (writes XML files to current directory)
    """
    if file_extensions is None:
        file_extensions = ['.endf', '.dat']  # Support both by default

    dirPath = Path(dirName)
    if not dirPath.exists():
        print(f"Error: Directory not found: {dirName}")
        return

    fileList = [f for f in dirPath.iterdir() if f.is_file()]

    for filepath in fileList:
        if filepath.suffix in file_extensions:
            print(f"Parsing: {filepath.name}")
            try:
                ENDF6FPYreader(str(filepath))
            except Exception as e:
                print(f"Error parsing {filepath.name}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python FPYParser.py <directory> [file_extension]")
        print("  directory: Path to directory containing FPY files")
        print("  file_extension: Optional, e.g., '.endf' or '.dat' (default: both)")
        sys.exit(1)

    dirName = sys.argv[1]

    # Determine file extension from command line or auto-detect
    if len(sys.argv) >= 3:
        ext = sys.argv[2]
        if not ext.startswith('.'):
            ext = '.' + ext
        file_extensions = [ext]
    else:
        # Auto-detect: if directory contains 'JEFF' or 'JENDL', use .dat, else .endf
        if 'JEFF' in dirName.upper() or 'JENDL' in dirName.upper():
            file_extensions = ['.dat']
            print("Detected JEFF/JENDL database - using .dat extension")
        else:
            file_extensions = ['.endf']
            print("Detected ENDF database - using .endf extension")

    parse_directory(dirName, file_extensions)
