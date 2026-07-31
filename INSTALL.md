# CONFLUX Installation Guide

Complete installation instructions for the CONFLUX neutrino flux calculation package.

## Table of Contents

- [Quick Start](#quick-start)
- [Detailed Installation](#detailed-installation)
- [Manual Installation](#manual-installation)
- [Troubleshooting](#troubleshooting)
- [Uninstallation](#uninstallation)

---

## Quick Start

For most users, this is all you need:

```bash
# Install from source
cd /path/to/conflux
pip install .

# Run setup wizard
conflux-setup

# Verify installation
python -c "import conflux; print('CONFLUX version:', conflux.__version__)"
```

That's it! The setup wizard will guide you through database configuration.

---

## Detailed Installation

### Prerequisites

- Python 3.6 or higher
- pip package manager
- Internet connection (for database downloads)
- ~50 MB free disk space

### Step 1: Install CONFLUX Package

#### Option A: From Source (Recommended for Development)

```bash
cd /path/to/conflux
pip install -e .  # Editable install for development
```

#### Option B: From Source (Standard)

```bash
cd /path/to/conflux
pip install .
```

#### Option C: From GitHub (Future)

```bash
pip install git+https://github.com/CNFLUX/conflux.git
```

#### Option D: From PyPI (Future)

```bash
pip install conflux
```

### Step 2: Run Setup Wizard

The setup wizard will configure your installation:

```bash
conflux-setup
```

The wizard will:

1. **Configure Database Path**
   - Default: `~/.conflux/data` or package installation directory
   - Custom: You can specify any directory

2. **Set Environment Variable**
   - Adds `CONFLUX_DB` to your shell configuration (`.zshrc`, `.bashrc`, etc.)
   - Makes databases accessible to CONFLUX

3. **Download ENDF Database**
   - Downloads ENDF-B-VIII.0 fission product yields (~3 MB)
   - Parses to XML format automatically
   - Required for fission calculations

4. **Optional: ENSDF Database**
   - Guide for installing ENSDF decay database
   - Useful for decay calculations
   - Can be downloaded automatically

5. **Optional: JEFF/JENDL Databases**
   - Setup JEFF-3.3 and JENDL-5 fission product yield databases
   - Requires manual download from NEA/JAEA
   - See [JEFF_JENDL_SETUP.md](JEFF_JENDL_SETUP.md) for details

6. **Optional: Covariance Matrices**
   - Downloads covariance matrices from FYCoM project (~300-500 MB)
   - Required for uncertainty calculations
   - See [COVARIANCE_MATRICES.md](COVARIANCE_MATRICES.md) for details

7. **Verify Installation**
   - Checks Python version
   - Verifies all dependencies
   - Confirms database setup

### Step 3: Activate Environment

Reload your shell configuration:

```bash
# For zsh
source ~/.zshrc

# For bash
source ~/.bashrc
```

Or simply restart your terminal.

### Step 4: Verify Installation

```bash
# Check CONFLUX can be imported
python -c "import conflux; print('Success!')"

# Check database path
echo $CONFLUX_DB

# Verify setup
conflux-setup --verify
```

---

## Manual Installation

If you prefer to set up manually without the wizard:

### 1. Install Package

```bash
pip install .
```

### 2. Create Database Directory

```bash
mkdir -p ~/.conflux/data
```

### 3. Set Environment Variable

Add to your shell config (`~/.zshrc` or `~/.bashrc`):

```bash
export CONFLUX_DB="$HOME/.conflux/data"
```

### 4. Download and Parse ENDF Database

```bash
# Using the provided script
cd /path/to/conflux
./update_endf_database.py --output-dir ~/.conflux/data/fissionDB/ENDF

# Or use the CLI tool
conflux-update-endf --output-dir ~/.conflux/data/fissionDB/ENDF
```

### 5. (Optional) Install ENSDF Database

Download ENSDF from https://www.nndc.bnl.gov/ensdf/ and place in:

```bash
~/.conflux/data/ENSDF/
```

---

## Installation Options

### Database Location Options

#### Option 1: User Home Directory (Recommended)

```bash
~/.conflux/data/
```

**Pros:**
- Separate from package (survives package updates)
- User-specific (no permission issues)
- Easy to back up

**Cons:**
- Not shared between users

#### Option 2: Package Directory

```bash
/path/to/site-packages/conflux/data/
```

**Pros:**
- Bundled with package
- Works immediately after install

**Cons:**
- Lost on package reinstall
- May require admin permissions

#### Option 3: System-Wide

```bash
/opt/conflux/data/
# or
/usr/local/share/conflux/data/
```

**Pros:**
- Shared between all users
- Central management

**Cons:**
- Requires admin permissions
- More complex setup

### Advanced Setup Options

#### Specify Custom Database Path

```bash
conflux-setup --set-db-path /custom/path/to/databases
```

#### Download ENDF Only

```bash
conflux-setup --download-endf --db-path ~/.conflux/data
```

#### Verify Installation Only

```bash
conflux-setup --verify
```

---

## Dependencies

CONFLUX will automatically install these Python packages:

- `numpy` - Numerical computations
- `scipy>=1.8.1` - Scientific computing
- `tqdm` - Progress bars
- `matplotlib` - Plotting
- `iminuit` - Minimization
- `fortranformat` - ENDF parsing
- `pandas` - Data manipulation
- `xraydb` - X-ray data

---

## Database Information

### ENDF-B-VIII.0 Database

- **Size:** ~3 MB (compressed), ~10 MB (extracted)
- **Contents:** Fission product yields for 31 isotopes
- **Format:** XML (parsed from ENDF-6 format)
- **Required for:** Fission reactor calculations
- **Auto-downloaded:** Yes (via setup wizard)

### ENSDF Database

- **Size:** ~100 MB
- **Contents:** Nuclear decay data (beta, gamma, alpha)
- **Format:** ENSDF text files
- **Required for:** Decay calculations, summation methods
- **Auto-downloaded:** No (must obtain from NNDC)

---

## Post-Installation

### Test Your Installation

Create a test script `test_conflux.py`:

```python
#!/usr/bin/env python3
import conflux
import numpy as np

# Test basic import
print(f"CONFLUX version: {conflux.__version__}")

# Test database path
import os
db_path = os.environ.get('CONFLUX_DB')
print(f"Database path: {db_path}")

# Test creating an engine
from conflux import BetaEngine

engine = BetaEngine()
print("BetaEngine created successfully!")

# Test basic calculation
e_spec = np.arange(0, 10, 0.1)
print(f"Energy spectrum points: {len(e_spec)}")

print("\n✓ All tests passed!")
```

Run it:

```bash
python test_conflux.py
```

### Example Usage

See the `examples/` directory for complete examples:

```bash
cd examples/
python BetaExample.py
python ReactorExample.py
```

---

## Troubleshooting

### Issue: Module Not Found

**Error:**
```
ModuleNotFoundError: No module named 'conflux'
```

**Solution:**
```bash
# Verify installation
pip list | grep conflux

# If not found, reinstall
pip install --force-reinstall .
```

### Issue: CONFLUX_DB Not Set

**Error:**
```
Warning: CONFLUX_DB environment variable not set
```

**Solution:**
```bash
# Check if set
echo $CONFLUX_DB

# If empty, add to shell config
echo 'export CONFLUX_DB="$HOME/.conflux/data"' >> ~/.zshrc
source ~/.zshrc
```

### Issue: Database Files Not Found

**Error:**
```
FileNotFoundError: Could not find ENDF database
```

**Solution:**
```bash
# Re-download database
conflux-setup --download-endf --db-path ~/.conflux/data

# Or manually
./update_endf_database.py --output-dir ~/.conflux/data/fissionDB/ENDF
```

### Issue: Permission Denied

**Error:**
```
PermissionError: [Errno 13] Permission denied
```

**Solution:**
```bash
# Use user install
pip install --user .

# Or use a virtual environment
python -m venv venv
source venv/bin/activate
pip install .
```

### Issue: Dependency Installation Failed

**Error:**
```
ERROR: Failed building wheel for scipy
```

**Solution:**
```bash
# Install build dependencies
pip install --upgrade pip setuptools wheel

# On macOS, install Xcode Command Line Tools
xcode-select --install

# On Linux, install build essentials
sudo apt-get install build-essential python3-dev
```

### Issue: Covariance Matrix Download Slow or Fails

**Problem:**
Covariance matrix download takes a long time or fails partway through.

**Solution:**
```bash
# Skip during initial setup, download later
conflux-setup  # Answer 'n' to covariance matrices

# Download covariance matrices separately when you need them
conflux-setup --download-covariance --db-path $CONFLUX_DB

# Or use the manual downloader script
python $CONFLUX_DB/CovMatDownloader.py
```

**Note:** Covariance matrices are large (~300-500 MB) but optional. Skip them if:
- You don't need uncertainty calculations
- You have slow internet connection
- You have limited disk space

See [COVARIANCE_MATRICES.md](COVARIANCE_MATRICES.md) for more details.

---

## Virtual Environment (Recommended)

Using a virtual environment keeps CONFLUX isolated:

```bash
# Create virtual environment
python -m venv conflux-env

# Activate it
source conflux-env/bin/activate  # Linux/Mac
conflux-env\Scripts\activate     # Windows

# Install CONFLUX
cd /path/to/conflux
pip install .

# Run setup
conflux-setup

# When done, deactivate
deactivate
```

---

## Uninstallation

To completely remove CONFLUX:

```bash
# 1. Uninstall package
pip uninstall conflux

# 2. Remove databases (optional)
rm -rf ~/.conflux

# 3. Remove environment variable from shell config
# Edit ~/.zshrc or ~/.bashrc and remove:
#   export CONFLUX_DB="..."

# 4. Reload shell
source ~/.zshrc
```

---

## Docker Installation (Future)

Coming soon: Docker container with pre-configured environment.

---

## Support

- **Issues:** https://github.com/CNFLUX/conflux/issues
- **Documentation:** See `README.md` and `examples/`
- **Email:** zhang39@llnl.gov

---

## Next Steps

After installation:

1. Read the `README.md` for package overview
2. Explore `examples/` directory for usage examples
3. Check `ENDF_UPDATE_README.md` for database management
4. Join the CONFLUX community on GitHub

---

**Installation successful? Start calculating neutrino fluxes!**
