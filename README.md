CONFLUX: A reactor neutrino flux calculation framework
======================================================
CONFLUX, Calculation Of Neutrino FLUX, is a framework that allow users to
calculate reactor neutrino flux with flexible and time dependent inputs of
reactor models. 

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Features
The framework provides three different modes of neutrino flux
calculation:
- Summation mode,
- Beta-conversion mode,
- Experiment mode.
  
## Installation:
Follow these steps to set up the project locally:
1. Clone the repository
2. Execute: 
`pip3 install ./conflux`
3. In the system environment setup, add:
`export CONFLUX_DB=</path/to/conflux>/data` on Linux or MacOS
`set CONFLUX_DB=</path/to/conflux>/data` on Windows 
this will setup the nuclear databases necessary for the reactor neutrino calculation.

Program structure:
==================

Databases:
-----------
cumulative fission yield, beta decay branch

Summation calculation modules:
------------------------------
fission yield tally, beta/neutrino spectral shape
calculation, spectrum summation

Usage:
======

Summation calculation:
----------------------
Define fission isotope and beta decaying isotope

Define reactor model with reactor power, isotope fractions and uncertainties of
fractions

Load beta decay database and calculate beta spectra.

Sum beta spectra with respect to total isotope fractions.
