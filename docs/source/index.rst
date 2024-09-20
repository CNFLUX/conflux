.. CONFLUX documentation master file, created by Anosh on Friday Aug 30 2024.


CONFLUX Software Package
========================

CONFLUX, Calculation Of Neutrino FLUX, is a framework that allow users to calculate reactor neutrino flux with flexible and time dependent inputs of reactor models. The framework provides three different modes of neutrino flux calculation:

|  Summation mode
|  Beta-conversion mode
|  Direct measurement mode
|
The formalism is described in this publication_, while the code is published here_.


.. _publication: www.google.com 
.. _here: https://github.com/CNFLUX/conflux

Installation & Basic execution
------------

.. toctree::
   :maxdepth: 2
   
   installation

Overview
--------

.. toctree::
   :maxdepth: 2
    
   Intro
   Modes
   Input
    

API Reference
-------------
.. toctree::
   :maxdepth: 2

   Basic
   BetaEngine
   FPYEngine
   SumEngine
   ConversionEngine


Status
------

.. include:: status.rst

Indices and tables
===================

* :ref:`genindex`
