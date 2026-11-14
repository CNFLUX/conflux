#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 10:51:51 2021

@author: zhang39
"""

# universal modules
import argparse
from pathlib import Path
import re
import tempfile
from xml.dom import minidom
import fortranformat as ff

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

    def saveXML(self, output_dir=None):
        for head in self.DB.childNodes:
            if len(head.childNodes) != int(head.getAttribute('LE')):
                raise ValueError('Parsed energy-table count does not match HEAD LE')
            for table in head.childNodes:
                if len(table.childNodes) != int(table.getAttribute('NFPi')):
                    raise ValueError('Parsed product count does not match LIST NFPi')
        xml_str = self.root.toprettyxml(indent ="\t")
        destination = Path(output_dir or '.') / (self.outputName + '.xml')
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w") as f:
            f.write(xml_str)
        return destination

# In[1]: the method that reads ENDF fission production yield (FPY) data base
def ENDF6FPYreader(filename, output_dir=None):
    record_reader = ff.FortranRecordReader('(A66,I4,I2,I3,I5)')
    data_reader = ff.FortranRecordReader('(6E11.0)')
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

                    xmloutput.createHEAD(ZA, AWR, LEplus1, 'IFP')

                else:

                        # reading CONT []
                    # Add bounds check to prevent index errors
                    if (i < len(NFPi) and len(CONT) <= NFPi[i] and  NFPi[i]>0):
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 4:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 2:
                            CONT[-1].append(data[0])
                            CONT[-1].append(data[1])
                            xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
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
                        xmloutput.createLIST(Ei[i],Ii[i],NFPi[i])

                linenum+=1

            # data type spliter
            if (linenum !=0 and MT == 0):
                linenum = 0

            # reading the HEAD of cumulative fission products
            if (MT == 459):
                data = data_reader.read(TEXT)
                if (linenum == 0):
                    print(data)
                    i = 0
                    CONT = []
                    ZA = int(data[0])
                    AWR = float(data[1])
                    LEplus1 = int(data[2])
                    Ei = [0]*LEplus1
                    Ii = [0]*LEplus1
                    NNi = [0]*LEplus1
                    NFPi = [0]*LEplus1

                    xmloutput.createHEAD(ZA, AWR, LEplus1, 'CFP')

                else:
                    # reading CONT
                    # Add bounds check to prevent index errors
                    if (i < len(NFPi) and len(CONT) <= NFPi[i] and  NFPi[i]>0):
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 4:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 2:
                            CONT[-1].append(data[0])
                            CONT[-1].append(data[1])
                            xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
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
                        xmloutput.createLIST(Ei[i],Ii[i],NFPi[i])

                linenum += 1

    return xmloutput.saveXML(output_dir)


def parse_fission_yield_library(filename, output_dir=None):
    """Split an ENDF library by MAT and use the existing per-material reader.

    Accepts combined JEFF .txt libraries as well as single-material .dat files.
    Output names follow CONFLUX's nfpy_<MAT>_<Z>-<element>-<A>[M] convention.
    Both MF=8/MT=454 (IFP) and MT=459 (CFP) are retained by the reader.
    """
    materials = {}
    for line in Path(filename).read_text().splitlines(keepends=True):
        if not line.strip():
            continue
        if len(line.rstrip('\r\n')) < 75:
            raise ValueError('ENDF record shorter than 75 columns')
        mat = int(line[66:70])
        if mat > 0:
            materials.setdefault(mat, []).append(line)

    outputs = []
    with tempfile.TemporaryDirectory(prefix='conflux-fpy-') as temporary:
        for mat, lines in materials.items():
            if not any(int(line[70:72]) == 8 and int(line[72:75]) in (454, 459)
                       for line in lines):
                continue
            description = [line for line in lines
                           if int(line[70:72]) == 1 and int(line[72:75]) == 451]
            if len(description) < 5:
                raise ValueError(f'MAT {mat}: missing MF=1/MT=451 isotope metadata')
            isotope = re.match(r'\s*(\d+)-\s*([A-Za-z]+)\s*-\s*(\d+)', description[4][:66])
            if isotope is None:
                raise ValueError(f'MAT {mat}: unrecognized isotope description')
            z, element, a = isotope.groups()
            state = int(description[1][33:44])  # LISO, not LIS or LFI
            suffix = '' if state == 0 else ('M' if state == 1 else f'M{state}')
            basename = f'nfpy_{mat}_{int(z)}-{element}-{int(a)}{suffix}'
            source = Path(temporary) / (basename + '.dat')
            source.write_text(''.join(lines))
            outputs.append(ENDF6FPYreader(str(source), output_dir=output_dir))
    if not outputs:
        raise ValueError(f'No MF=8 fission-yield materials found in {filename}')
    return outputs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse JEFF ENDF-6 fission yields into CONFLUX XML.')
    parser.add_argument('input', type=Path, help='Combined library file or directory of .dat files')
    parser.add_argument('--output-dir', type=Path, default=Path('.'), help='XML destination (default: current directory)')
    args = parser.parse_args()
    sources = sorted(args.input.glob('*.dat')) if args.input.is_dir() else [args.input]
    outputs = []
    for source in sources:
        outputs.extend(parse_fission_yield_library(source, output_dir=args.output_dir))
    print(f'Saved {len(outputs)} XML files to {args.output_dir.resolve()}')
