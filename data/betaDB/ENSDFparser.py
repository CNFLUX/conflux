#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Xianyi Zhang, LLNL <zhang39@llnl.gov>
# License: MIT License
#
# Description:
#   Generate CONFLUX xml format nuclear structure dataset for B- DECAY ONLY
#
# Usage:
#   python ENSDFparser.py <ENSDF folder path>

# # universal modules
import sys
from xml.dom import minidom
from os import listdir, environ
import csv
import numpy as np
import re
import pkg_resources
from tqdm import tqdm

# define time units to convert all half-lives in seconds
tu = {'NS':1e-9, 'MS':1e-3, 'S':1}
tu['M'] = tu['S']*60
tu['H'] = tu['M']*60
tu['D'] = tu['H']*24
tu['Y'] = tu['D']*365
print(tu)

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
    # listname = pkg_resources.resource_filename('conflux', 'betaDB/Z_to_element.csv')
    with open(environ["CONFLUX_DB"]+"/betaDB/Z_to_element.csv") as csvinput:
        csvreader = csv.reader(csvinput, delimiter=',')
        for row in csvreader:
            if row[0].isdigit:
                zdict[row[1]] = row[0]
    return zdict

# check if a string contain any element in a list
def contains_any(string, elements):
    for element in elements:
        if element in string:
            return True
    return False

# global dictionary of element and Z
elementdict = element_to_Z()

# global method to quantify simplified uncertainty indication
def transUncert(mean, unc):
    if (not unc.isdigit()):
        return 0
    decimal = mean.find(".")
    exponent = mean.find("E")
    magnitude = len(mean[decimal:(exponent-1)]) if exponent>0 else len(mean[decimal:(exponent)])
    correction = int(mean[exponent+1:]) if exponent>0 else 0
    return float(unc)*pow(10, -magnitude+correction)

# global method to convert spin to numerical veriable
def convert_J(chars):

    chars = chars.replace('(', '')
    chars = chars.replace(')', '')
    chars = chars.replace('[', '')
    chars = chars.replace(']', '')

    tags = ['LT', 'GT', 'LE', 'GE', 'AP', 'CA', 'SY', '|>', '<|', '< ', ' >']
    tagged = False
    whichtag = ''
    for tag in tags:
        if tag in chars:
            tagged = True
            whichtag += tag
    if tagged:
        pi = 1
        if '-' in chars:
            pi -= 1
            chars = chars.replace('-','')
            chars = chars.replace('+','')
        else:
            pi += 1
            chars = chars.replace('+','')

        chars = chars.replace(whichtag, "")
        print(whichtag)
        if '/' in chars:
            num, denom = chars.split('/')
            num = num.strip()
            denom = denom.strip()
            J = (float(num)/float(denom))
        else:
            chars = chars.strip()
            J = (float(chars))

        if whichtag in ['LT', '<']:
            J -= 1
        elif whichtag in ['GT', '>']:
            J += 1
        return [np.array([pi]),np.array([J]) ]

    links = ['TO', 'AND', ':', "&"]
    linked = False
    whichlink = ''
    for link in links:
        if link in chars:
            linked = True
            whichlink += link
    if linked:
        pi = 1
        if '-' in chars:
            pi -= 1
            chars = chars.replace('-','')
            chars = chars.replace('+','')
        else:
            pi += 1
            chars = chars.replace('+','')

        Js = chars.split(whichlink)
        Jlist = np.zeros(len(Js))
        it = 0
        for Jstr in Js:
            if '/' in Jstr:
                num, denom = Jstr.split('/')
                num = num.strip()
                denom = denom.strip()

                Jlist[it] = (float(num)/float(denom))
            else:
                Jstr = Jstr.strip()
                Jlist[it] = (float(Jstr))
            it +=1

        J = np.arange(float(Jlist[0]), float(Jlist[1]), 1)

        pi = np.zeros(len(J))+pi

        return [pi, J]

    J = np.zeros(len(chars.split(',')))
    pi = np.zeros(len(chars.split(',')))
    if not any(c.isdigit() for c in chars):
        if '-' in chars:
            pi += -1
        return [pi, np.array([1e3])]

    it = 0
    for Js in chars.split(','):
        # for n+1/2 spin angular momentum
        if '/' in chars:
            num, denom = Js.split('/')
            if '-' in denom:
                denom = denom.replace('-','')
                pi[it] = (-1)
            elif '+' in denom:
                denom = denom.replace('+','')
                pi[it] = (1)
            else:
                pi[it] = -1 if ('-' in Js) else 1
            J[it] = (float(num)/float(denom))
        # for n spin angular momentum
        else:
            pi[it] = -1 if ('-' in Js) else 1
            num = Js.replace('-', '') if ('-' in Js) else Js.replace('+', '')
            J[it] = num
        it +=1
    return [pi, J]

# prepare the xml data file structrue
class XMLedit:
    def __init__(self):
        self.root = minidom.Document()
        self.DB = self.root.createElement('betaDB')
        self.root.appendChild(self.DB)

    # beta-decay isotopes
    def createIsotope(self, istpName, isotopeID, Q = '0.0', HL = 0.0, Ex = '0.0'):
        if HL == 0: return
        self.isotope = self.root.createElement('isotope')
        self.isotope.setAttribute('name', str(istpName))
        self.isotope.setAttribute('isotope', str(isotopeID))
        self.isotope.setAttribute('Q', str(Q))
        self.isotope.setAttribute('HL', str(HL))
        self.isotope.setAttribute('Ex', str(Ex))
        self.DB.appendChild(self.isotope)

    # decay branches of the isotope
    def editBranch(self, fraction, end_point_E, dJpi, sigma_frac='0.0', sigma_E0='0.0'):
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

# class to save parent isotope information
class ParentIstp:
    def __init__(self, line):
        self.A = int(line[:3].strip())
        element = line[3:5].capitalize().strip()
        self.name = element+line[:3].strip()
        self.Z = int(elementdict[element])
        leveltxt = line[9:19].strip()

        for a in leveltxt:
            if a in ["X", "Y", "Z"]:
                print("+"+a)
                leveltxt=line[9:19].strip().replace("+"+a, "")
                print(leveltxt)
                break

        print(leveltxt, line[9:19].strip().isdigit())
        if not any(char.isdigit() for char in leveltxt):
            leveltxt = '0'
        if not float(leveltxt):
            leveltxt = '0'
        self.level = float(leveltxt)/1e3#convert to MeV
        self.HL = 0.
        hl_str = line[39:49]
        hl_array = hl_str.split()
        if hl_array:
            life = float(hl_array[0])
            unit = tu[hl_array[1]]
            self.HL = life*unit
        self.Emax = float(line[64:74].strip())/1e3+float(leveltxt)/1e3
        self.d_Emax = transUncert(line[64:74].strip(), line[74:76].strip())/1e3
        self.pi, self.J = convert_J(line[21:39])

# class to save information of the level before decay
class DecayLevel:
    def __init__(self, line):
        if (is_number(line[9:19].strip())):
            self.level = float(line[9:19].strip())/1e3 #convert to MeV
        elif (is_number(line[9:19].strip().replace('+X',''))):
            self.level = float(line[9:19].strip().replace('+X',''))/1e3
        else:
            self.level = 0
        self.d_level = transUncert(line[9:19].strip(), line[19:21].strip())/1e3 if is_number(line[19:21]) else 0
        # hl_str = line[39:49]
        # hl_array = hl_str.split()
        # life = hl_array[0]
        # unit = tu[hl_array[1]]
        # self.HL = life*unit
        self.pi, self.J = convert_J(line[21:39])

# class to save information of decay branches
class DecayBranch:
    def __init__(self, line):
        self.E0 = float(line[9:19].strip())/1e3 if is_number(line[9:19]) else 0
        self.sigma_E0 = transUncert(line[9:19].strip(), line[19:21].strip())/1e3 if is_number(line[19:21]) else 0
        self.fraction = float(line[21:29].strip())/100 if is_number(line[21:29]) else 1
        self.sigma_frac = transUncert(line[21:29].strip(), line[29:31].strip())/100 if is_number(line[29:31]) else 0
        self.forbidden = line[77:79].strip() if not line[77:79].isspace() else "0"

def ENSDFbeta(fileList):
    xmloutput = XMLedit()
    for filename in tqdm(fileList):
        inputfile = open(dirName+filename, "r", errors='replace')

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
            # betatag = ["B- DECAY", 
            #            "B-N DECAY", 
            #            "B-2N DECAY", 
            #            "B-P DECAY",
            #            "B-_ DECAY"]
            betatag = ["B- DECAY"]
            excludetag = ["B- DECAY:NEUTRON", "2B- DECAY"] # providing a list of exceptions
            if MT == "   " and contains_any(line, betatag) and not contains_any(line, excludetag):
                lastline = line
                betabool = True

            if betabool == True:
                # save information of the parent isotope
                if MT == "  P":
                    print(lastline)
                    print(line)
                    decayparent = ParentIstp(line)
                    if 1e3 in decayparent.J:
                        decayparent.J = np.array([decayparent.A/2 - int(decayparent.A/2)])
                        print(line, decayparent.J)
                    if decayparent.level > beginlevel:
                        beginlevel = decayparent.level
                        isomer+=1
                    else:
                        beginlevel = 0
                        isomer = 0
                    ZAI = int(decayparent.Z*1e4 + decayparent.A*10 + isomer)
                    # save the parent isotope information
                    istpname = decayparent.name
                    if isomer>0: istpname = decayparent.name+f"m{isomer}"
                    xmloutput.createIsotope(decayparent.name, ZAI, decayparent.Emax, decayparent.HL, decayparent.level)
                    lastline = line
                elif MT == "  L":
                    lastline = line
                elif MT == "  B" and "  B" not in lastline:
                    decaybranch = DecayBranch(line)
                    decaydaughter = DecayLevel(lastline)
                    if decaybranch.E0 == 0:
                        decaybranch.E0 = decayparent.Emax - decaydaughter.level
                        decaybranch.sigma_E0 = decayparent.d_Emax + decaydaughter.d_level
                    # calculate angular momentum difference
                    #for j in decaydaughter.J:
                    for j in range(len(decaydaughter.J)):
                        if decaydaughter.J[j] == 1e3:
                            decaydaughter.J[j] = decayparent.J.min()
                    decaybranch.forbidden = [decayparent.pi.min()*decaydaughter.pi, abs(decayparent.J.min() - decaydaughter.J)]
                    forbid_str = ""
                    for i in range(len(decaybranch.forbidden[1])):
                        dspin_str = str(int(decaybranch.forbidden[1][i]))
                        if decaybranch.forbidden[0][i] < 0:
                            dspin_str = dspin_str+'-'
                        if i == 0:
                            forbid_str = forbid_str+dspin_str
                        else:
                            forbid_str = forbid_str+","+dspin_str
                    # save the decay branch information
                    xmloutput.editBranch(str("{:.3f}".format(decaybranch.fraction)), 
                                         str("{:.3f}".format(decaybranch.E0)), 
                                         forbid_str, 
                                         str("{:.3f}".format(decaybranch.sigma_frac)), 
                                         str("{:.3f}".format(decaybranch.sigma_E0)))
                    lastline = line
                    #print(lastline)

    xmloutput.saveXML("ENSDFbetaDB3.xml")

    return 0


if __name__ == "__main__":
    dirName = sys.argv[1]
    fileList = listdir(dirName)
    fileList.sort()
    ENSDFbeta(fileList)
