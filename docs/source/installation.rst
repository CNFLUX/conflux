Installation & Basic execution
******************************

Installation & Dependencies
===========================

CONFLUX is ready to use from download and does not need to be built, however it does require python_ version 3.6 or greater.

CONFLUX is installed using pip_. 

.. _python: https://www.python.org/
.. _pip: https://pypi.org/project/pip/

Dependencies - Python Libraries
-------------------------------

CONFLUX requires the following packages to be installed

- Numpy_
- Scipy_ (version 1.8.1 or greater)
- TQDM_
- matplotlib_
- iminuit_
- fortranformat_


.. _Numpy: https://numpy.org/
.. _Scipy: https://scipy.org/
.. _TQDM: https://github.com/tqdm/tqdm
.. _matplotlib: https://matplotlib.org/
.. _iminuit: https://pypi.org/project/iminuit/
.. _fortranformat: https://pypi.org/project/fortranformat/


All of these can be installed using ``pip`` with the included requirements.txt file. 

.. code-block:: bash

    pip install -r requirements.txt


Note that all of these libraries will be automatically installed on installation of CONFLUX. 

Databases
---------

CONFLUX comes with up to date Nuclear data databases for Fission Yields (ENDF/JEFF), Beta Decay (ENSDF), and for Covariance Matrices (FYCOM), all of which can be found here:

- ENDF_
- JEFF_
- ENSDF_
- FYCOM_

.. _ENDF: https://www.nndc.bnl.gov/endf/
.. _JEFF:
.. _ENSDF: https://www.nndc.bnl.gov/ensdf/
.. _FYCOM: https://nucleardata.berkeley.edu/FYCoM/index.html

More on what each database does in the context of CONFLUX can be found in the introduction_ of this manual.
Additionally, ``ENSDFparser.py`` , ``betdaDBxmlGen.py`` , and ``CovMatDownloader.py`` are all files that either download these databases to CONFLUX, or 
help to make the databases more readable. 

.. _introduction: 


Installation
------------

CONFLUX can be installed using pip

.. code-block:: bash
    pip install /path/to/conflux

Make sure that after installation, the ``CONFLUX_DB`` environmental variable points to the ``data`` directory in the CONFLUX installation.
Alternatively, you can define it in your ``.bashrc`` file

Execution
=========

Execution of CONFLUX is performed as with any other python file. an example running a Summation Engine is given below.

.. code-block:: bash
    python3 SumEngineExample.py

Running this specific example will output a ``.png`` file with a U235 neutrino spectrum. How to extract other spectral information can be found within the various examples provided with CONFLUX.
