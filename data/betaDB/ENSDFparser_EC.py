#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Xianyi Zhang, LLNL <zhang39@llnl.gov>
# License: MIT License
#
# Description:
#   Generate CONFLUX xml format nuclear structure dataset for EC and B+ DECAY ONLY
#
# Usage:
#   python ENSDFparser.py <ENSDF folder path>

# # universal modules

# universal modules
import sys
from xml.dom import minidom
from os import listdir, environ
import csv
import numpy as np
import re
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
def convert_J2(chars):

    chars = chars.replace('(', '')
    chars = chars.replace(')', '')
    chars = chars.replace('[', '')
    chars = chars.replace(']', '')
    chars = chars.replace('J', '')

    tags = ['LT', 'GT', 'LE', 'GE', 'AP', 'CA', 'SY', '|>', '<|', '< ', ' >']
    links = ['TO', 'AND', ':', "&", "to", ","]

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
        # print(whichtag)

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

            parts = re.split(r'(?i)\bTO\b|AND|:|&|,', chars)
            Jlist = []
            for p in parts:
                p = p.strip()
                if not p:
                    continue
                if '/' in p:
                    num, denom = p.split('/')
                    Jlist.append(float(num)/float(denom))
                else:
                    Jlist.append(float(p))

            J = np.arange(float(Jlist[0]), float(Jlist[1]), 1)

            pi = np.zeros(len(J))+pi

            return [pi, J]
        else:
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

        parts = re.split(r'(?i)\bTO\b|AND|:|&|,', chars)
        Jlist = []
        for p in parts:
            p = p.strip()
            if not p:
                continue
            if '/' in p:
                num, denom = p.split('/')
                Jlist.append(float(num)/float(denom))
            else:
                Jlist.append(float(p))

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

# 1) Define your tags and linking‐word patterns once
TAG_LIST    = ['LT','GT','LE','GE','AP','CA','SY','|>','<|','<','>']
# sort by length so that "<|" is found before "<"
TAG_PATTERN = '|'.join(map(re.escape, sorted(TAG_LIST, key=len, reverse=True)))

# splits on TO, AND (case‐insensitive) or on comma, colon or ampersand
LINK_PATTERN = r'(?i)\b(?:TO|AND)\b|[:,&]'

from fractions import Fraction
def parse_tok(tok: str) -> tuple[float,int]:
    """
    Split off any leading or trailing +/- as parity,
    then parse the remainder as an integer or fraction.
    
    Returns (magnitude, parity)
      e.g. "3/2-" -> (1.5, -1)
            "-5/3+"-> (5/3, -1)
            "+4"   -> (4.0, +1)
    """
    tok = tok.strip()
    parity = 1
    
    # leading sign?
    if tok and tok[0] in '+-':
        if tok[0] == '-': parity = -parity
        tok = tok[1:]
    
    # trailing sign?
    if tok and tok[-1] in '+-':
        if tok[-1] == '-': parity = -parity
        tok = tok[:-1]
    
    # now tok is purely digits or a/b
    mag = float(Fraction(tok))
    return mag, parity

def convert_J(chars: str) -> tuple[np.ndarray,np.ndarray]:
    """
    Parse a spin‐parity Jπ spec.  
    Special cases:
      • empty → (pi=[+1], J=[1000.0])
      • single '+' or '-' → (pi=[±1], J=[1.0])
    """

    # 1) strip parentheses, brackets, literal 'J'
    s = re.sub(r'[()\[\]J]', '', chars).strip()
    # print(f"s = {s}")

    # 2) empty‐field guard
    if not s:
        return np.array([+1], dtype=int), np.array([1e3], dtype=float)

    # 3) single‐sign guard → force J=1.0
    if s in ('+', '-'):
        sign = +1 if s == '+' else -1
        return np.array([sign], dtype=int), np.array([1.0], dtype=float)

    # 4) pull off any comparison tag
    m   = re.search(TAG_PATTERN, s)
    tag = m.group(0) if m else ''
    if tag:
        s = s.replace(tag, '').strip()

    # 5) tokenise
    parts = [p.strip() for p in re.split(LINK_PATTERN, s) if p.strip()]
    if not parts:
        # e.g. someone wrote only a tag like "<"
        return np.array([+1], dtype=int), np.array([1e3], dtype=float)

    # 6) parse each into (magnitude, parity)
    mags, parities = zip(*(parse_tok(p) for p in parts))

    # 7) assemble as before
    # if len(mags) == 2:
    #     start, end = mags
    #     parity      = parities[0]
    #     # print(mags) 
    #     pi          = np.full(int((end - start)/1.0), parity, dtype=int)
    #     J           = np.arange(start, end, 1.0)

    if len(mags) == 1:
        Jval = mags[0]
        if tag in ('LT','<'): Jval -= 1
        if tag in ('GT','>'): Jval += 1

        pi = np.array([parities[0]], dtype=int)
        J  = np.array([Jval], dtype=float)

    else:
        # literal list of levels
        pi = np.array(parities, dtype=int)
        J  = np.array(mags, dtype=float)

    return pi, J

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
    def editBranch(self, IB, IE, end_point_E, dJpi, sigma_IB='0.0', sigma_IE='0.0', sigma_E0='0.0'):
        branch = self.root.createElement('branch')
        branch.setAttribute('IB', str(IB))
        branch.setAttribute('sigma_IB', str(sigma_IB))
        branch.setAttribute('IE', str(IE))
        branch.setAttribute('sigma_IE', str(sigma_IE))
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

        leveltxt = re.sub(r"[^0-9.]", "", leveltxt)

        # print(leveltxt, line[9:19].strip().isdigit())
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
        self.E0 = float(line[9:19].strip())/1e3 if is_number(line[9:19]) else 0.
        self.sigma_E0 = transUncert(line[9:19].strip(), line[19:21].strip())/1e3 if is_number(line[19:21]) else 0.
        self.IB = float(line[21:29].strip())/100 if is_number(line[21:29]) else 0.
        self.IE = float(line[31:39].strip())/100 if is_number(line[31:39]) else 0.
        # if not is_number(line[21:29]) and not is_number(line[31:39]):
        #     self.IB = 0.5
        #     self.IE = 0.5
        self.sigma_IB = transUncert(line[21:29].strip(), line[29:31].strip())/100 if is_number(line[29:31]) else 0
        self.sigma_IE = transUncert(line[31:39].strip(), line[39:41].strip())/100 if is_number(line[39:41]) else 0

        self.forbidden = line[77:79].strip() if not line[77:79].isspace() else "0"

def ENSDFbeta(fileList):
    xmloutput = XMLedit()
    for filename in tqdm(fileList):
        inputfile = open(f"{dirName}/{filename}", "r", errors='replace')

        lastline = ""
        betabool = False
        beginlevel = 0.
        isomer = 0
        for line in inputfile.readlines():

            if line.isspace():
                # if (betabool == True): print("EMPTY LINE")
                lastline = ""
                betabool = False

            # find parent isotope of beta decay (ignore other types of decay)
            MT = line[5:8]      # Check datatype
            # betatag = ["B- DECAY", 
            #            "B-N DECAY", 
            #            "B-2N DECAY", 
            #            "B-P DECAY",
            #            "B-_ DECAY"]
            betatag = ["EC DECAY",
                       "B+ DECAY"]
            excludetag = [] # providing a list of exceptions
            if MT == "   " and contains_any(line, betatag) and not contains_any(line, excludetag):
                lastline = line
                betabool = True

            if betabool == True:
                # save information of the parent isotope
                if MT == "  P":
                    decayparent = ParentIstp(line)
                    if decayparent.name == "Sr73": 
                        print(lastline)
                        print(line)
                    if 1e3 in decayparent.J:
                        decayparent.J = np.array([decayparent.A/2 - int(decayparent.A/2)])
                        # print(line, decayparent.J)
                    if decayparent.level > beginlevel:
                        beginlevel = decayparent.level
                        isomer+=1
                    else:
                        beginlevel = 0
                        isomer = 0
                    ZAI = int(decayparent.Z*1e4 + decayparent.A*10 + isomer)
                    # save the parent isotope information
                    istpname = decayparent.name
                    if isomer>0: istpname = istpname+f"m{isomer}"
                    xmloutput.createIsotope(istpname, ZAI, decayparent.Emax, decayparent.HL, decayparent.level)
                    lastline = line
                elif MT == "  L":
                    if decayparent.name == "Sr73": 
                        print(line)
                    lastline = line
                elif MT == "  E" and "  E" not in lastline:
                    if decayparent.name == "Sr73": 
                        print(line)
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
                    if any([decaybranch.IB, decaybranch.IE, decaybranch.sigma_IB, decaybranch.sigma_IE]) != 0:
                        xmloutput.editBranch(str("{:.4f}".format(decaybranch.IB)),
                                            str("{:.4f}".format(decaybranch.IE)), 
                                            str("{:.4f}".format(decaybranch.E0)), 
                                            forbid_str, 
                                            str("{:.4f}".format(decaybranch.sigma_IB)), 
                                            str("{:.4f}".format(decaybranch.sigma_IE)), 
                                            str("{:.4f}".format(decaybranch.sigma_E0)))
                    lastline = line
                    #print(lastline)

    xmloutput.saveXML("ENSDFbetaDB_EC.xml")

    return 0


if __name__ == "__main__":
    dirName = sys.argv[1]
    fileList = listdir(dirName)
    fileList.sort()
    ENSDFbeta(fileList)
