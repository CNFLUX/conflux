intro
*****


How Does CONFLUX Work?
======================

The CONFLUX (Calculation Of Neutrino FLUX) software framework is built with the goal of simplifying and standardizing neutrino flux calculations. CONFLUX is packaged with three prediction modes:

- Summation mode
- :math:`\beta` Conversion mode
- Direct Experimental Measurement mode

All of which calculate the Neutrino spectrum in different ways. A block diagram of how the calculation is run is provided below.

.. image:: block_diagram.jpg

Modes
=====

Summation
---------

Summation mode is an ab-initio calculation that takes in Fission product information, either from ``ENDF``, ``JEFF``, or a ``user-defined DB``, and combines it with the spectral shape of each individual 
:math:`\beta`-branch. Thus, we sum the product of the individual branch spectra and their contributions to form the total neutrino spectrum for a given isotope. A block diagram of how the mode works,
as well as a graphical representation of the calculation is provided below.

.. image:: Summation_block.jpg

.. image:: Summation_figure.png







Databases
---------


Beta Spectrum Generator
-----------------------
