![GitHub](https://img.shields.io/github/license/CNFLUX/conflux)

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
  
## Installation

### Quick Start

For most users, this is all you need:

```bash
# 1. Install CONFLUX
pip install .

# 2. Run setup wizard
conflux-setup

# 3. Verify installation
python -c "import conflux; print('CONFLUX installed successfully!')"
```

The `conflux-setup` wizard will guide you through database configuration automatically.

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

## Databases

CONFLUX uses nuclear databases for reactor neutrino flux calculations:

### Fission Product Yields
- **[ENDF](https://www.nndc.bnl.gov/endf-releases/?version=B-VIII.1)** - Primary database (auto-downloaded by setup)
- **[JEFF](https://www.oecd-nea.org/dbdata/jeff/jeff33/index.html)** - Alternative European database
- **[JENDL](https://wwwndc.jaea.go.jp/jendl/jendl.html)** - Alternative Japanese database

### Beta Decay Data
- **[ENSDF](https://www.nndc.bnl.gov/ensdfarchivals/)** - Evaluated Nuclear Structure Data File

### Covariance Data
- **[FYCoM](https://nucleardata.berkeley.edu/FYCoM/)** - Fission Yield Covariance Matrices

The `conflux-setup` wizard handles database downloads and parsing automatically. All databases are stored in `$CONFLUX_DB` in XML format.

For advanced usage and manual parsing, see [INSTALL.md](INSTALL.md).

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
- `conflux.ConversionEngine` converts beta spectra of fissile isotopes to the corresponding neutrino spectra using best fit virtual beta branches
 
### Examples

A large list of example python scripts are saved in `<conflux>/examples/`. Users can run the examples or write calculation programs based off the examples for most common reactor neutrino production modeling. 

### Documentation

Documentation from code comments is built using `Sphinx` (https://www.sphinx-doc.org)

Install if needed: `pip3 install -U sphinx sphinx-autoapi`

```
cd <conflux>/docs

# print list of options
make

# HTML pages --- open _build/html/index.html in browser
make html

# PDF --- requires `latexmk`
make latexpdf
```

## Contributing:

CONFLUX is distributed under the terms of the MIT license. All new contributions must be made under this license.
To contribute to the CONFLUX project, fork the repository and submit the pull request. 

For more information and collaboration, please contact Xianyi Zhang (zhang39@llnl.gov).

## Acknowledgement:

This work was supported by the Lawrence Livermore National Laboratory LDRD Program under Project No. 20-SI-005, the U.S. Department of Energy Office of Science, Office of High Energy Physics under Award No. DE-SC0020262 to Virginia Polytechnic Institute and State University and under Work Proposal Number SCW1504 to Lawrence Livermore National Laboratory, and by the U.S. Department of Energy Office of Defense Nuclear Nonproliferation Research and Development.  This work was supported by the Consortium for Monitoring, Technology, and Verification under DOE-NNSA award number DE-NA0003920. The authors thank Daniel Nestares from the University of California, Merced, for his work to test the software. The authors thank Mitchel Crockett from the University of Tennessee, Knoxville, for his reactor simulation output to aid CONFLUX on understanding the reactor simulation input. At last, The authors thank Eric F. Matthew from University of California, Berkeley for providing a reference fission product covariance dataset. This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. 

## License:

CONFLUX is distributed under the terms of the **MIT License** . All new contributions must be made under this license.

See [LICENSE](https://github.com/CNFLUX/conflux/LICENSE) and [NOTICE](https://github.com/CNFLUX/conflux/conflux/NOTICE) for details.

`SPDX-License-Identifier: MIT`

``LLNL-CODE-2003431``

