Installation & Basic execution
******************************

Installation & Dependencies
===========================

CONFLUX is ready to use from download and does not need to be built, however it does require python_ version 3.6 or greater

Databases are automatically included in every CONFLUX installation. however, you can download the Fission Covariance matrices here_, 
the Fission Yield Databases here_, and the Beta Decay Databases here_. 

.. _python: https://www.python.org/
.. _pip: https://pypi.org/project/pip/
.. _here: https://nucleardata.berkeley.edu/FYCoM/index.html
.. _here: https://www.nndc.bnl.gov/ensdf/
.. _here: https://www.nndc.bnl.gov/endf/

Dependencies - Python Libraries
-------------------------------

CONFLUX requires the following packages to be installed

- Numpy_
- Scipy_ (version 1.8.1 or greater)
- TQDM_
- matplotlib_
- iminuit_
- fortranformat_

All of these can be installed using ''pip'' with the name of the library to be installed

.. code-block:: bash

    pip install "Library"

Note, all of these libraries will be automatically installed on installation of CONFLUX. 

Installation
------------

CONFLUX can be installed using pip

.. code-block:: bash
    pip install ./conflux

Make sure that after installation, the ``CONFLUX_DB`` environmental variable points to the data directory in the CONFLUX installation.
Alternatively, you can define it in your ``.bashrc`` file

Execution
=========

Execution of CONFLUX is performed as with any other python file

.. code-block:: bash
    python3 
