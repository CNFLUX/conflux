.. CONFLUX documentation master file, created by Anosh on Friday Aug 30 2024.


CONFLUX Software Package
========================

CONFLUX, Calculation Of Neutrino FLUX, is a framework that allows users to calculate reactor neutrino flux and uncertainties with flexible and time dependent inputs of reactor models. It uses the most up-to-date Nuclear libraries, as well as a robust Beta spectrum shape generator to create a neutrino spectrum.

The formalism is described in this publication_, while the code is published here_. Finally, this manual_ should be used as a general reference for how to use the framework.

.. _publication:
.. _here: https://github.com/CNFLUX/conflux/tree/master
.. _manual: https://conflux.readthedocs.io/en/latest/


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
