CONFLUX: A reactor neutrino flux calculation framework
======================================================

CONFLUX, Calculation Of Neutrino FLUX, is a framework that allow users to
calculate reactor neutrino flux with flexible and time dependent inputs of
reactor models. The framework provides three different modes of neutrino flux
calculation:
- Summation mode,
- Beta-conversion mode,
- Experiment mode.

Installation:
=============
Execute: 
`pip3 install ./conflux`
In the system environment setup, add:
`export CONFLUX_DB=</path/to/conflux>/data`, OR<br>
`set CONFLUX_DB=</path/to/conflux>/data`
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
