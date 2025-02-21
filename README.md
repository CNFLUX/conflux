CONFLUX: A reactor neutrino flux calculation framework
======================================================
CONFLUX, Calculation Of Neutrino FLUX, is a framework that allow users to
calculate reactor neutrino flux with flexible and time dependent inputs of
reactor models. 

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Databases](#database)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Features
The framework provides three different modes of neutrino flux
calculation:
- Summation mode,
- Beta-conversion mode,
- Neutrino data mode.
  
## Installation:
Follow these steps to set up the project locally:
1. Clone the repository
2. Execute:
`cd <cloned_repository>`\
`pip3 install ./conflux`
4. In the system environment setup, such as `$HOME/.bashrc`, add:\
`export CONFLUX_DB=</path/to/conflux>/data` on Linux or MacOS\
`set CONFLUX_DB=</path/to/conflux>/data` on Windows\
this will setup the nuclear databases necessary for the reactor neutrino calculation.

## Databases:
CONFLUX provides python executables to download and parse nucleaer databases including 
[ENDF](https://www.nndc.bnl.gov/endf-releases/?version=B-VIII.1), 
[JEFF](https://www.oecd-nea.org/dbdata/jeff/jeff33/index.html), 
[JENDL](https://wwwndc.jaea.go.jp/jendl/jendl.html),
for fission for fission product yield calculation, and
[ENSDF](https://www.nndc.bnl.gov/ensdfarchivals/),
for beta decay and neutrino spectrum measurement.
All databases are saved in the `$CONFLUX_DB` folder, in `xml` format. For databases different from the CONFLUX prepackaged version, download the database and run:\
`python3 $CONFLUX_DB/fissionDB/ENDF/FPYParserENDF.py <ENDF-6 format database folder>`\
to parse fission product yield data into the CONFLUX xml format, and\
`python3 $CONFLUX_DB/betaDB/ENSDFparser.py <ENSDF database folder>`\
to parse beta decay data into the CONFLUX xml format.

CONFLUX also uses correlation and covariance matrix from [FYCOM](https://nucleardata.berkeley.edu/FYCoM/). Run\
`python3 CovMatDownloader.py` to download the database. This data is in csv file due to its size.

##  Usage:
### Executable
CONFLUX contains a python executable at `<conflux>/exec/quickflux.py`, which takes `json` file as macros to calculate reactor or beta decay neutrino productions with simple source term configurations. An example `<conflux>/exec/example.json` contains all basic json entries and sequences needed to execute the calculation by running:
`python3 <conflux>/exec/quickflux.py <conflux>/exec/example.json`

### Libraries
Users can import CONFLUX libraries in their own python scripts for neutrino flux calculations. The major libraries include:
- `conflux.bsg`: the beta spectrum generation functions through beta decay calculation with theoretical corrections \
- `conflux.BetaEngine` tallies beta decay branches of beta-unstable isotopes to calculate each individual beta/neutrino spectrum\
- `conflux.FPYEngine` tallies fission products to calculate the spectrum and uncertainty of each individual fissile isotopes\
- `conflux.SumEngine` sums neutrino/beta spectra with respect to fission fraction and non-fissile contributions in a reactor model\
- `conflux.Conversion` converts beta spectra of fissile isotopes to the corresponding neutrino spectra using best fit virtual beta branches
 
### Examples
A large list of example python scripts are saved in `<conflux>/examples/`. Users can run the examples or write calculation programs based off the examples for most common reactor neutrino production modeling. 

## Contributing:

## License:

## Acknowledgement:
This work was supported by the Lawrence Livermore National Laboratory LDRD Program under Project No. 20-SI-005, the U.S. Department of Energy Office of Science, Office of High Energy Physics under Award No. DE-SC0020262 to Virginia Polytechnic Institute and State University and under Work Proposal Number SCW1504 to Lawrence Livermore National Laboratory, and by the U.S. Department of Energy Office of Defense Nuclear Nonproliferation Research and Development.  This work was supported by the Consortium for Monitoring, Technology, and Verification under DOE-NNSA award number DE-NA0003920. The authors thank Daniel Nestares from the University of California, Merced, for his work to test the software. The authors thank Mitchel Crockett from the University of Tennessee, Knoxville, for his reactor simulation output to aid CONFLUX on understanding the reactor simulation input. At last, The authors thank Eric F. Matthew from University of California, Berkeley for providing a reference fission product covariance dataset. This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. LLNL-SW-872132
