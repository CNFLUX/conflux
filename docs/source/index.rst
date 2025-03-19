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


Status
------

.. include:: status.rst

Acknowledgement
===============
This work was supported by the Lawrence Livermore National Laboratory LDRD Program under Project No. 20-SI-005, the U.S. Department of Energy Office of Science, Office of High Energy Physics under Award No. DE-SC0020262 to Virginia Polytechnic Institute and State University and under Work Proposal Number SCW1504 to Lawrence Livermore National Laboratory, and by the U.S. Department of Energy Office of Defense Nuclear Nonproliferation Research and Development.  This work was supported by the Consortium for Monitoring, Technology, and Verification under DOE-NNSA award number DE-NA0003920. The authors thank Daniel Nestares from the University of California, Merced, for his work to test the software. The authors thank Mitchel Crockett from the University of Tennessee, Knoxville, for his reactor simulation output to aid CONFLUX on understanding the reactor simulation input. At last, The authors thank Eric F. Matthew from University of California, Berkeley for providing a reference fission product covariance dataset. This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. LLNL-SM-872132

Indices and tables
===================

* :ref:`genindex`
