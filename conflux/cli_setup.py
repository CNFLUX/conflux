#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONFLUX Setup and Configuration Tool

This module provides interactive setup for CONFLUX installation,
including database downloads, environment configuration, and verification.

Author: Xianyi Zhang, LLNL
Date: 2026-07-16
"""

import os
import sys
import argparse
import json
from pathlib import Path
import urllib.request
import zipfile
import shutil
import subprocess

try:
    import fortranformat
    HAS_PARSER_DEPS = True
except ImportError:
    HAS_PARSER_DEPS = False

# ANSI color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_header(text):
    """Print colored header"""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text.center(70)}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}\n")

def print_success(text):
    """Print success message"""
    print(f"{Colors.OKGREEN}✓ {text}{Colors.ENDC}")

def print_error(text):
    """Print error message"""
    print(f"{Colors.FAIL}✗ {text}{Colors.ENDC}")

def print_warning(text):
    """Print warning message"""
    print(f"{Colors.WARNING}⚠ {text}{Colors.ENDC}")

def print_info(text):
    """Print info message"""
    print(f"{Colors.OKBLUE}ℹ {text}{Colors.ENDC}")

def get_user_input(prompt, default=None):
    """Get user input with optional default"""
    if default:
        prompt = f"{prompt} [{default}]: "
    else:
        prompt = f"{prompt}: "

    response = input(prompt).strip()
    return response if response else default

def get_yes_no(prompt, default=True):
    """Get yes/no response from user"""
    default_str = "Y/n" if default else "y/N"
    response = input(f"{prompt} [{default_str}]: ").strip().lower()

    if not response:
        return default
    return response in ['y', 'yes']

def parse_endf_inline(endf_dir, output_dir):
    """
    Inline ENDF parser as fallback when script is not found
    Uses the unified FPYParser module
    """
    if not HAS_PARSER_DEPS:
        print_error("fortranformat not installed - cannot parse ENDF files")
        print_info("Install with: pip install fortranformat")
        return False

    try:
        # Try to import the unified parser from the data directory
        sys.path.insert(0, str(Path(__file__).parent / "data" / "fissionDB"))

        try:
            from FPYParser import ENDF6FPYreader

            endf_files = sorted(Path(endf_dir).glob("*.endf"))
            parsed_count = 0
            failed_count = 0

            # Change to output directory so XML files are saved there
            original_cwd = os.getcwd()
            os.chdir(output_dir)

            for endf_file in endf_files:
                try:
                    ENDF6FPYreader(str(endf_file))
                    parsed_count += 1
                except Exception as e:
                    print_warning(f"Failed to parse {endf_file.name}")
                    failed_count += 1

            os.chdir(original_cwd)

            if parsed_count > 0:
                print_success(f"Parsed {parsed_count}/{len(endf_files)} files")
                return True
            else:
                print_error("No files were successfully parsed")
                return False

        except ImportError:
            # Fall back to basic XML creation if parser not available
            print_warning("FPYParser not found, creating placeholder XMLs")

            endf_files = sorted(Path(endf_dir).glob("*.endf"))
            for endf_file in endf_files:
                # Create minimal XML as placeholder
                from xml.dom import minidom
                root = minidom.Document()
                db = root.createElement(endf_file.stem)
                db.setAttribute("status", "placeholder")
                root.appendChild(db)

                output_path = output_dir / f"{endf_file.stem}.xml"
                with open(output_path, "w") as f:
                    f.write(root.toprettyxml(indent="\t"))

            print_warning("Created placeholder XMLs - run conflux-update-endf to properly parse")
            return True

    except Exception as e:
        print_error(f"Inline parsing error: {e}")
        return False


def copy_package_data_to_db(db_path):
    """
    Copy ALL git-tracked data files from package installation to CONFLUX_DB

    This ensures basic functionality even if downloads/parsing fails.

    Copies everything from the installation's data/ directory:
    - XML files (fission yield data, beta decay data) - ALL git-tracked files
    - CSV files (reference data, example models) - small files only
    - Parser scripts (.py files - needed for parsing ENSDF)

    Skips only:
    - Covariance matrices (*cov*.csv, *corr*.csv) - these are huge (~3.5 GB)
    - __pycache__ and temporary files

    For editable installs: Copies from git repo
    For PyPI installs: No data to copy (will download instead)

    Returns:
        bool: True if successful, False if no data available
    """
    print_info("Copying all git-tracked data from installation to CONFLUX_DB...")

    try:
        # Find package data directory (works for both editable and regular installs)
        import conflux

        if conflux.__file__ is None:
            print_error("Cannot locate conflux installation")
            return False

        package_dir = Path(conflux.__file__).parent / "data"

        if not package_dir.exists():
            print_warning(f"No data directory found at installation location: {package_dir}")
            print_info("For PyPI installs, data will be downloaded via setup wizard")
            return False

        print_success(f"Found data directory: {package_dir}")

        # Copy entire data directory structure, excluding only huge covariance files
        total_copied = 0
        total_skipped = 0

        for src_file in package_dir.rglob("*"):
            if not src_file.is_file():
                continue

            # Skip covariance/correlation matrices (huge files)
            # These are identified by having "cov" or "corr" in the filename
            if src_file.suffix == ".csv":
                if "cov" in src_file.name.lower() or "corr" in src_file.name.lower():
                    total_skipped += 1
                    continue

            # Skip __pycache__ and bytecode
            if "__pycache__" in src_file.parts or src_file.suffix in [".pyc", ".pyo"]:
                continue

            # Compute destination path
            rel_path = src_file.relative_to(package_dir)
            dest_file = db_path / rel_path

            # Create destination directory
            dest_file.parent.mkdir(parents=True, exist_ok=True)

            # Skip if already exists and same size
            if dest_file.exists() and dest_file.stat().st_size == src_file.stat().st_size:
                continue

            # Copy file
            shutil.copy2(src_file, dest_file)
            total_copied += 1

        if total_copied > 0:
            print_success(f"Copied {total_copied} data files from installation")
            if total_skipped > 0:
                print_info(f"Skipped {total_skipped} large covariance files (download separately if needed)")
            return True
        else:
            print_info("All data files already present in CONFLUX_DB")
            return True

    except Exception as e:
        print_error(f"Failed to copy data: {e}")
        import traceback
        print_error(traceback.format_exc())
        return False


def get_default_db_path():
    """
    Get default database path

    Always returns the installation's data directory if it exists.
    This way CONFLUX_DB defaults to where the data is bundled.

    For editable installs: Returns git repo's data directory
    For regular installs: Returns site-packages/conflux/data
    Fallback: ~/.conflux/data (only if no installation found)
    """
    try:
        import conflux

        if conflux.__file__ is None:
            return Path.home() / ".conflux" / "data"

        package_dir = Path(conflux.__file__).parent
        data_dir = package_dir / "data"

        # Always use installation's data directory if it exists
        if data_dir.exists():
            return data_dir
        else:
            # No data bundled - fallback to home directory
            # (This shouldn't happen with dumb-proof installation)
            return Path.home() / ".conflux" / "data"
    except:
        # Fallback to user's home directory
        return Path.home() / ".conflux" / "data"

def setup_environment_variable(db_path):
    """Setup CONFLUX_DB environment variable"""
    print_info("Setting up CONFLUX_DB environment variable...")

    shell = os.environ.get('SHELL', '/bin/bash')
    shell_name = Path(shell).name

    # Determine shell config file
    home = Path.home()
    if shell_name == 'zsh':
        config_file = home / ".zshrc"
    elif shell_name == 'bash':
        config_file = home / ".bashrc"
        if not config_file.exists():
            config_file = home / ".bash_profile"
    else:
        config_file = home / ".profile"

    export_line = f'export CONFLUX_DB="{db_path}"'

    # Check if already set
    if config_file.exists():
        with open(config_file, 'r') as f:
            content = f.read()
            if 'CONFLUX_DB' in content:
                print_warning(f"CONFLUX_DB already set in {config_file}")
                if not get_yes_no("Overwrite?", default=False):
                    return False

    # Add to config file
    try:
        with open(config_file, 'a') as f:
            f.write(f"\n# CONFLUX database path\n")
            f.write(f"{export_line}\n")

        # Also set for current session
        os.environ['CONFLUX_DB'] = str(db_path)

        print_success(f"Added CONFLUX_DB to {config_file}")
        print_info(f"Run: source {config_file}")
        return True
    except Exception as e:
        print_error(f"Failed to update {config_file}: {e}")
        return False

def download_endf_database(db_path, version="VIII.0"):
    """Download ENDF database"""
    print_info(f"Downloading ENDF-B-{version} database...")

    db_url = f"https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-{version}_nfy.zip"
    temp_dir = Path("/tmp/conflux_setup_endf")
    temp_dir.mkdir(exist_ok=True)

    zip_path = temp_dir / f"ENDF-B-{version}_nfy.zip"

    try:
        # Download with progress
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            sys.stdout.write(f"\r  Progress: {percent:.1f}% ({downloaded / 1024 / 1024:.1f} MB)")
            sys.stdout.flush()

        urllib.request.urlretrieve(db_url, zip_path, reporthook=report_progress)
        print()  # New line after progress
        print_success("Download complete")

        # Extract
        print_info("Extracting files...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)

        # Find .endf files
        endf_files = list(temp_dir.rglob("*.endf"))
        if not endf_files:
            print_error("No .endf files found in archive")
            return False

        endf_dir = endf_files[0].parent
        print_success(f"Found {len(endf_files)} ENDF files")

        # Parse to XML
        print_info("Parsing ENDF files to XML...")
        output_dir = db_path / "fissionDB" / "ENDF"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Find the update_endf_database script
        # Try multiple possible locations
        possible_paths = [
            Path(__file__).parent / "update_endf_database.py",  # In package directory
            Path(__file__).parent.parent / "update_endf_database.py",  # In site-packages
            Path(__file__).parent.parent.parent / "update_endf_database.py",  # In lib
            Path.cwd() / "update_endf_database.py",  # Current directory
        ]

        parser_script = None
        for path in possible_paths:
            if path.exists():
                parser_script = path
                break

        if parser_script:
            result = subprocess.run([
                sys.executable,
                str(parser_script),
                "--no-download",
                "--input-dir", str(endf_dir),
                "--output-dir", str(output_dir)
            ], capture_output=True, text=True)

            if result.returncode == 0:
                print_success(f"Parsed files to {output_dir}")
            else:
                print_error(f"Parsing failed: {result.stderr}")
                return False
        else:
            print_warning(f"Parser script not found in standard locations")
            print_info("Trying inline parser...")
            # Fall back to inline parsing
            success = parse_endf_inline(endf_dir, output_dir)
            if success:
                print_success(f"Parsed files to {output_dir}")
            else:
                print_error("Inline parsing failed")
                return False

        # Cleanup
        shutil.rmtree(temp_dir)

        return True

    except Exception as e:
        print_error(f"Failed to download ENDF database: {e}")
        return False


def download_jeff_database(db_path, version="3.3"):
    """
    Download JEFF database automatically

    Downloads JEFF-3.3 fission product yields from NEA Data Bank
    """
    print_info(f"Downloading JEFF-{version} fission product yields...")

    # NEA Data Bank URLs for JEFF-3.3 fission yields
    base_url = "https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/"
    files_to_download = [
        ("JEFF33-nfy.asc", "Neutron-induced fission yields"),
        ("JEFF33-sfy.asc", "Spontaneous fission yields"),
    ]

    temp_dir = Path("/tmp/conflux_setup_jeff")
    temp_dir.mkdir(exist_ok=True)

    try:
        downloaded_files = []

        # Download both files
        for filename, description in files_to_download:
            url = base_url + filename
            local_path = temp_dir / filename

            print_info(f"Downloading {description}...")

            def report_progress(block_num, block_size, total_size):
                downloaded = block_num * block_size
                if total_size > 0:
                    percent = min(100, downloaded * 100 / total_size)
                    sys.stdout.write(f"\r  Progress: {percent:.1f}% ({downloaded / 1024 / 1024:.2f} MB)")
                else:
                    sys.stdout.write(f"\r  Downloaded: {downloaded / 1024 / 1024:.2f} MB")
                sys.stdout.flush()

            urllib.request.urlretrieve(url, local_path, reporthook=report_progress)
            print()  # New line after progress
            downloaded_files.append(local_path)

        print_success(f"Downloaded {len(downloaded_files)} JEFF files")

        # Parse to XML
        jeff_output = db_path / "fissionDB" / "JEFF"
        jeff_output.mkdir(parents=True, exist_ok=True)

        print_info(f"Parsing JEFF files to {jeff_output}...")

        # Use JEFF parser
        parser_path = Path(__file__).parent / "data" / "fissionDB" / "FPYParser.py"

        if not parser_path.exists():
            print_error("JEFF parser not found")
            print_info(f"Expected at: {parser_path}")
            print_warning("Files downloaded to {temp_dir} but not parsed")
            return False

        # JEFF parser expects .dat files in a directory
        # Rename .asc files to .dat
        for jeff_file in downloaded_files:
            new_name = jeff_file.parent / (jeff_file.stem + '.dat')
            jeff_file.rename(new_name)

        parsed_count = 0
        original_cwd = os.getcwd()
        os.chdir(jeff_output)

        # JEFF parser expects a directory path as argument
        result = subprocess.run([
            sys.executable,
            str(parser_path),
            str(temp_dir) + "/"  # Parser expects directory with trailing slash
        ], capture_output=True, text=True, timeout=120)

        os.chdir(original_cwd)

        # Count XML files generated
        xml_files = list(jeff_output.glob("*.xml"))
        parsed_count = len(xml_files)

        if result.returncode != 0 and parsed_count == 0:
            print_warning("Parser execution had issues")
            if result.stderr:
                print(f"    Error: {result.stderr[:300]}")
        elif parsed_count > 0:
            print_success(f"Parsed {parsed_count} JEFF files successfully")

        # Cleanup
        shutil.rmtree(temp_dir)

        if parsed_count > 0:
            print_success(f"Successfully set up JEFF-{version} database")
            return True
        else:
            print_error("No files were successfully parsed")
            return False

    except urllib.error.HTTPError as e:
        print_error(f"Download failed: HTTP {e.code}")
        print_info("URL may have changed. Check: https://www.oecd-nea.org/dbdata/jeff/jeff33/")
        return False
    except Exception as e:
        print_error(f"Failed to download JEFF database: {e}")
        return False


def download_jendl_database(db_path, version="5"):
    """
    Download JENDL database automatically

    Downloads JENDL-5 fission product yields from JAEA
    Note: Uses SSL context that disables certificate verification due to JAEA's self-signed cert
    """
    print_info(f"Downloading JENDL-{version} fission product yields...")

    # JAEA URL for JENDL-5 fission yields
    jendl_url = f"https://wwwndc.jaea.go.jp/ftpnd/ftp/JENDL/jendl{version}-fpy_upd8.tar.gz"

    temp_dir = Path("/tmp/conflux_setup_jendl")
    temp_dir.mkdir(exist_ok=True)
    tar_path = temp_dir / f"jendl{version}-fpy.tar.gz"

    try:
        print_warning("Note: JAEA uses self-signed SSL certificates")
        print_info("Creating SSL context with certificate verification disabled...")

        # Create SSL context that accepts self-signed certificates
        # This is necessary for JAEA's website
        import ssl
        ssl_context = ssl.create_default_context()
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.CERT_NONE

        print_info(f"Downloading from JAEA...")

        # Download with progress and custom SSL context
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                percent = min(100, downloaded * 100 / total_size)
                sys.stdout.write(f"\r  Progress: {percent:.1f}% ({downloaded / 1024 / 1024:.2f} MB)")
            else:
                sys.stdout.write(f"\r  Downloaded: {downloaded / 1024 / 1024:.2f} MB")
            sys.stdout.flush()

        # Use custom opener with SSL context
        opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl_context))
        urllib.request.install_opener(opener)

        urllib.request.urlretrieve(jendl_url, tar_path, reporthook=report_progress)
        print()  # New line after progress
        print_success("Download complete")

        # Extract tar.gz
        print_info("Extracting files...")
        import tarfile
        with tarfile.open(tar_path, 'r:gz') as tar_ref:
            tar_ref.extractall(temp_dir)

        # Find JENDL files (they use ENDF-6 format)
        jendl_files = []
        for ext in ['*.txt', '*.dat', '*.endf', '']:
            jendl_files.extend(temp_dir.rglob(ext) if ext else [])

        # Filter to actual data files (not README, etc.)
        jendl_files = [f for f in jendl_files if f.is_file() and f.stat().st_size > 1000]

        if not jendl_files:
            print_error("No JENDL data files found in archive")
            return False

        print_success(f"Found {len(jendl_files)} JENDL files")

        # Parse using ENDF parser (JENDL uses ENDF-6 format)
        jendl_output = db_path / "fissionDB" / "JENDL"
        jendl_output.mkdir(parents=True, exist_ok=True)

        print_info(f"Parsing JENDL files to {jendl_output}...")
        print_info("JENDL uses ENDF-6 format - using ENDF parser")

        if not HAS_PARSER_DEPS:
            print_error("fortranformat not installed - cannot parse JENDL files")
            print_info("Install with: pip install fortranformat")
            return False

        parsed_count = 0
        original_cwd = os.getcwd()
        os.chdir(jendl_output)

        # Import unified FPY parser
        fpy_parser_path = Path(__file__).parent / "data" / "fissionDB"
        sys.path.insert(0, str(fpy_parser_path))

        try:
            from FPYParser import ENDF6FPYreader

            for jendl_file in jendl_files:
                try:
                    ENDF6FPYreader(str(jendl_file))
                    parsed_count += 1
                    print_success(f"Parsed {jendl_file.name}")
                except Exception as e:
                    print_warning(f"Failed to parse {jendl_file.name}: {str(e)[:100]}")
        except ImportError as e:
            print_error(f"Could not import unified FPY parser: {e}")
            os.chdir(original_cwd)
            return False

        os.chdir(original_cwd)

        # Cleanup
        shutil.rmtree(temp_dir)

        if parsed_count > 0:
            print_success(f"Successfully set up JENDL-{version} database ({parsed_count} files)")
            return True
        else:
            print_error("No files were successfully parsed")
            return False

    except urllib.error.URLError as e:
        if "CERTIFICATE_VERIFY_FAILED" in str(e):
            print_error("SSL certificate verification failed (expected with self-signed cert)")
            print_info("The SSL workaround may not be working correctly")
        else:
            print_error(f"Download failed: {e}")
        print_info("Check: https://wwwndc.jaea.go.jp/jendl/j5/j5.html")
        return False
    except Exception as e:
        print_error(f"Failed to download JENDL database: {e}")
        import traceback
        print_warning(f"Details: {traceback.format_exc()[:300]}")
        return False


def find_latest_ensdf_url():
    """
    Find the latest ENSDF distribution URL from NNDC

    Returns:
        str: URL of latest ENSDF distribution
    """
    import datetime

    # Known working URLs to try (most recent first)
    fallback_urls = [
        "https://www.nndc.bnl.gov/ensdfarchivals/distributions/dist26/ensdf_260707.zip",
        "https://www.nndc.bnl.gov/ensdfarchivals/distributions/dist26/ensdf_260601.zip",
        "https://www.nndc.bnl.gov/ensdfarchivals/distributions/dist25/ensdf_251201.zip",
    ]

    # Start from current date and work backwards
    today = datetime.date.today()

    # Try current year and previous year (limited search to save time)
    for year_offset in range(0, 2):  # Try current year and 1 year back
        year = today.year - year_offset
        year_short = year % 100  # Get last 2 digits

        # Try months from most recent backwards
        start_month = 12 if year_offset > 0 else today.month
        for month in range(start_month, max(start_month - 6, 0), -1):  # Last 6 months
            # Try first and 15th of each month
            for day in [15, 1, 7]:  # Common release dates
                # Skip future dates
                try:
                    check_date = datetime.date(year, month, day)
                except ValueError:
                    continue

                if check_date > today:
                    continue

                # Construct URL
                date_str = f"{year_short:02d}{month:02d}{day:02d}"
                url = f"https://www.nndc.bnl.gov/ensdfarchivals/distributions/dist{year_short:02d}/ensdf_{date_str}.zip"

                try:
                    # Quick HEAD request to check if URL exists
                    req = urllib.request.Request(url, method='HEAD')
                    with urllib.request.urlopen(req, timeout=3) as response:
                        if response.status == 200:
                            return url
                except:
                    continue

    # Try fallback URLs
    for url in fallback_urls:
        try:
            req = urllib.request.Request(url, method='HEAD')
            with urllib.request.urlopen(req, timeout=3) as response:
                if response.status == 200:
                    return url
        except:
            continue

    # Last resort
    return fallback_urls[0]


def download_ensdf_database(db_path):
    """Download ENSDF database"""
    print_info("ENSDF database setup...")
    print_info("ENSDF contains nuclear decay data (beta, gamma, alpha)")
    print_warning("ENSDF database is ~30-50 MB compressed")
    print()

    # Ask if user wants to download
    if get_yes_no("Download ENSDF database from NNDC?", default=True):
        print_info("Finding latest ENSDF distribution...")

        # Find the latest ENSDF URL
        ensdf_url = find_latest_ensdf_url()

        if ensdf_url:
            # Extract date from URL for display
            import re
            date_match = re.search(r'ensdf_(\d{6})\.zip', ensdf_url)
            if date_match:
                date_str = date_match.group(1)
                year = "20" + date_str[0:2]
                month = date_str[2:4]
                day = date_str[4:6]
                print_success(f"Found ENSDF distribution: {year}-{month}-{day}")

            print_info("Downloading ENSDF database from NNDC...")
        temp_dir = Path("/tmp/conflux_setup_ensdf")
        temp_dir.mkdir(exist_ok=True)
        zip_path = temp_dir / "ensdf.zip"

        try:
            # Download with progress
            def report_progress(block_num, block_size, total_size):
                downloaded = block_num * block_size
                if total_size > 0:
                    percent = min(100, downloaded * 100 / total_size)
                    sys.stdout.write(f"\r  Progress: {percent:.1f}% ({downloaded / 1024 / 1024:.1f} MB)")
                else:
                    sys.stdout.write(f"\r  Downloaded: {downloaded / 1024 / 1024:.1f} MB")
                sys.stdout.flush()

            print_info("This may take a few minutes...")
            urllib.request.urlretrieve(ensdf_url, zip_path, reporthook=report_progress)
            print()  # New line after progress
            print_success("Download complete")

            # Extract
            print_info("Extracting files...")
            extract_dir = temp_dir / "ensdf_extracted"
            extract_dir.mkdir(exist_ok=True)

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

            # Find ENSDF files (ensdf.001, ensdf.002, etc.)
            ensdf_files = list(extract_dir.rglob("ensdf.*"))
            # Filter to only numbered files
            ensdf_files = [f for f in ensdf_files if f.stem.split('.')[-1].isdigit() or f.stem == 'ensdf']

            if ensdf_files:
                ensdf_dir = ensdf_files[0].parent
                print_success(f"Found {len(ensdf_files)} ENSDF files")

                # Copy to database location
                dest = db_path / "ENSDF"
                if dest.exists():
                    print_warning(f"ENSDF directory already exists at {dest}")
                    if not get_yes_no("Replace?", default=False):
                        shutil.rmtree(temp_dir)
                        return True
                    shutil.rmtree(dest)

                # Copy the directory
                shutil.copytree(ensdf_dir, dest)
                print_success(f"Installed ENSDF to {dest}")

                # Cleanup
                shutil.rmtree(temp_dir)

                # Offer to parse ENSDF to betaDB
                print()
                print_info("ENSDF has been downloaded and installed")
                if get_yes_no("Parse ENSDF to beta decay databases (betaDB)?", default=True):
                    parse_ensdf_to_betadb(dest, db_path)

                return True
            else:
                print_error("No ENSDF files found in archive")
                shutil.rmtree(temp_dir)
                # Fall back to manual setup
                print_info("Trying manual setup...")
                if get_yes_no("Do you have a local copy of ENSDF?", default=False):
                    return setup_ensdf_from_local(db_path)
                return False

        except urllib.error.URLError as e:
            print_error(f"Download failed: {e}")
            print_warning("Could not download from NNDC")
            print_info("You can download manually from: https://www.nndc.bnl.gov/ensdfarchivals/")

            # Offer manual setup option
            if get_yes_no("Do you have a local copy of ENSDF?", default=False):
                return setup_ensdf_from_local(db_path)
            return False

        except zipfile.BadZipFile as e:
            print_error(f"Invalid ZIP file: {e}")
            # Fall back to manual
            if get_yes_no("Do you have a local copy of ENSDF?", default=False):
                return setup_ensdf_from_local(db_path)
            return False

        except Exception as e:
            print_error(f"Failed to setup ENSDF: {e}")
            # Offer manual option
            if get_yes_no("Do you have a local copy of ENSDF?", default=False):
                return setup_ensdf_from_local(db_path)
            return False

    # Manual setup option
    elif get_yes_no("Do you have a local copy of ENSDF?", default=False):
        return setup_ensdf_from_local(db_path)

    else:
        print_info("Skipping ENSDF database setup")
        print_info("ENSDF is optional - CONFLUX will work without it")
        print_info("You can set it up later with: conflux-setup --download-ensdf")
        return True


def setup_ensdf_from_local(db_path):
    """Setup ENSDF from local directory"""
    ensdf_path = get_user_input("Enter path to ENSDF directory")
    if ensdf_path:
        ensdf_path = Path(ensdf_path).expanduser()
        if ensdf_path.exists():
            dest = db_path / "ENSDF"
            if dest.exists():
                print_warning(f"ENSDF directory already exists at {dest}")
                if not get_yes_no("Replace?", default=False):
                    return True
                shutil.rmtree(dest)

            shutil.copytree(ensdf_path, dest)
            print_success(f"Copied ENSDF to {dest}")

            # Offer to parse ENSDF to betaDB
            print()
            print_info("ENSDF has been copied to the database directory")
            if get_yes_no("Parse ENSDF to beta decay databases (betaDB)?", default=True):
                return parse_ensdf_to_betadb(dest, db_path)

            return True
        else:
            print_error(f"Directory not found: {ensdf_path}")
            return False

    print_info("Skipping ENSDF database setup")
    return True


def _copy_default_betadb_as_fallback(db_path, version_id=None):
    """
    Copy existing beta decay databases from installation as fallback

    This is called when ENSDF parsing fails. It attempts to copy any
    existing beta databases from the installation's data directory.

    For editable installs: Copies from git repo's data/betaDB/
    For PyPI installs: No defaults available (returns False gracefully)

    Args:
        db_path: CONFLUX database path
        version_id: Version identifier (if None, uses any available database)

    Returns:
        bool: True if successful
    """
    try:
        import conflux

        # Find installation's data directory
        install_data_dir = Path(conflux.__file__).parent / "data" / "betaDB"

        if not install_data_dir.exists():
            print_warning("No data directory at installation location")
            print_info("For PyPI installs, you need to download ENSDF and parse manually")
            return False

        # Find any existing beta database XMLs
        xml_files = list(install_data_dir.glob("ENSDF_betaDB_*.xml"))

        if not xml_files:
            print_warning("No beta database XML files found in installation")
            print_info("Run conflux-setup again and download ENSDF to parse")
            return False

        # Create output directory
        output_dir = db_path / "betaDB"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Copy all available beta databases
        copied = 0
        for src_file in xml_files:
            dest_file = output_dir / src_file.name

            # Don't overwrite if already exists
            if dest_file.exists():
                continue

            shutil.copy2(src_file, dest_file)
            print_success(f"Copied: {src_file.name}")
            copied += 1

        if copied > 0:
            print_warning("Using databases from installation location as fallback")
            print_info("Tip: To use latest ENSDF, download from NNDC and parse again")
            return True
        else:
            print_info("Beta databases already exist in CONFLUX_DB")
            return True

    except Exception as e:
        print_error(f"Failed to copy fallback databases: {e}")
        import traceback
        print_error(traceback.format_exc())
        return False


def parse_ensdf_to_betadb(ensdf_dir, db_path, version_id=None):
    """
    Parse ENSDF files to beta decay databases and update config.py

    Args:
        ensdf_dir: Directory containing ENSDF files
        db_path: CONFLUX database path
        version_id: Version identifier (e.g., "260707"), auto-detected if None

    Returns:
        bool: True if successful
    """
    print_info("Parsing ENSDF to beta decay databases...")

    # Auto-detect version from directory name if not provided
    if version_id is None:
        import re
        import datetime
        # Try to extract date from directory name (e.g., "ensdf_260707" -> "260707")
        match = re.search(r'(\d{6})', str(ensdf_dir))
        if match:
            version_id = match.group(1)
        else:
            # Fall back to current date
            today = datetime.date.today()
            version_id = f"{today.year % 100:02d}{today.month:02d}{today.day:02d}"

        print_info(f"Detected version: {version_id}")

    # Output paths
    beta_db_output = db_path / "betaDB" / f"ENSDF_betaDB_{version_id}.xml"
    beta_ec_output = db_path / "betaDB" / f"ENSDF_betaDB_EC_{version_id}.xml"

    # Create output directory
    beta_db_output.parent.mkdir(parents=True, exist_ok=True)

    # Find the parser scripts
    parser_dir = Path(__file__).parent / "data" / "betaDB"
    b_parser = parser_dir / "ENSDFparser.py"
    ec_parser = parser_dir / "ENSDFparser_EC.py"

    if not b_parser.exists():
        print_error(f"B- parser not found: {b_parser}")
        return False

    if not ec_parser.exists():
        print_error(f"EC/B+ parser not found: {ec_parser}")
        return False

    # Parse B- decay
    print_info("Parsing B- decay database...")
    try:
        result = subprocess.run([
            sys.executable,
            str(b_parser),
            str(ensdf_dir),
            str(beta_db_output)
        ], capture_output=True, text=True, timeout=300)

        if result.returncode == 0:
            print_success(f"Created: {beta_db_output.name}")
        else:
            print_error(f"B- parsing failed: {result.stderr[:200]}")
            print_info("Will use default database from package as fallback")
            return _copy_default_betadb_as_fallback(db_path, version_id)
    except Exception as e:
        print_error(f"B- parsing error: {e}")
        print_info("Will use default database from package as fallback")
        return _copy_default_betadb_as_fallback(db_path, version_id)

    # Parse EC/B+ decay
    print_info("Parsing EC/B+ decay database...")
    try:
        result = subprocess.run([
            sys.executable,
            str(ec_parser),
            str(ensdf_dir),
            str(beta_ec_output)
        ], capture_output=True, text=True, timeout=300)

        if result.returncode == 0:
            print_success(f"Created: {beta_ec_output.name}")
        else:
            print_error(f"EC/B+ parsing failed: {result.stderr[:200]}")
            print_info("Will use default database from package as fallback")
            return _copy_default_betadb_as_fallback(db_path, version_id)
    except Exception as e:
        print_error(f"EC/B+ parsing error: {e}")
        print_info("Will use default database from package as fallback")
        return _copy_default_betadb_as_fallback(db_path, version_id)

    # Update config.py with new version
    print()
    if get_yes_no(f"Update config.py to use version {version_id}?", default=True):
        success = update_config_beta_version(version_id)
        if success:
            print_success(f"Updated config.py: BETA_DB_VERSION = \"{version_id}\"")
            print_info("All CONFLUX engines will now use the new database automatically!")
        else:
            print_warning("Could not update config.py automatically")
            print_info(f"Manually update conflux/config.py:")
            print_info(f"  BETA_DB_VERSION = \"{version_id}\"")

    return True


def update_config_beta_version(version_id):
    """
    Update BETA_DB_VERSION in config.py

    Args:
        version_id: New version string (e.g., "270101")

    Returns:
        bool: True if successful
    """
    try:
        # Find config.py
        import conflux
        config_path = Path(conflux.__file__).parent / "config.py"

        if not config_path.exists():
            print_error(f"config.py not found: {config_path}")
            return False

        # Read current config
        with open(config_path, 'r') as f:
            content = f.read()

        # Replace version line
        import re
        pattern = r'BETA_DB_VERSION\s*=\s*["\'](\d{6})["\']'

        if not re.search(pattern, content):
            print_error("Could not find BETA_DB_VERSION in config.py")
            return False

        # Replace with new version
        new_content = re.sub(
            pattern,
            f'BETA_DB_VERSION = "{version_id}"',
            content
        )

        # Write back
        with open(config_path, 'w') as f:
            f.write(new_content)

        return True

    except Exception as e:
        print_error(f"Failed to update config.py: {e}")
        return False


def setup_jeff_database(db_path):
    """
    Setup JEFF-3.3 fission product yield database

    JEFF data must be downloaded manually from NEA website
    """
    print_info("JEFF-3.3 database setup...")
    print_info("JEFF contains fission product yields from NEA Data Bank")
    print_warning("JEFF database must be downloaded manually from NEA")
    print()

    print_info("To download JEFF-3.3:")
    print("  1. Visit: https://www.oecd-nea.org/dbdata/jeff/jeff33/")
    print("  2. Look for 'Fission Yields' section")
    print("  3. Download the JEFF-3.3 fission yield files")
    print("  4. Extract to a local directory")
    print()

    if get_yes_no("Do you have JEFF-3.3 files downloaded locally?", default=False):
        jeff_path = get_user_input("Enter path to JEFF-3.3 directory (containing .jeff files)")
        if not jeff_path:
            print_info("Skipping JEFF database setup")
            return True

        jeff_path = Path(jeff_path).expanduser()

        if not jeff_path.exists():
            print_error(f"Directory not found: {jeff_path}")
            return False

        # Check for JEFF files
        jeff_files = list(jeff_path.glob("*.jeff"))
        if not jeff_files:
            # Try other possible extensions
            jeff_files = list(jeff_path.glob("*.asc")) + list(jeff_path.glob("*.dat"))

        if not jeff_files:
            print_error("No JEFF format files found in directory")
            print_info("Expected files with .jeff, .asc, or .dat extensions")
            return False

        print_success(f"Found {len(jeff_files)} JEFF files")

        # Setup parser
        jeff_output = db_path / "fissionDB" / "JEFF"
        jeff_output.mkdir(parents=True, exist_ok=True)

        parser_path = Path(__file__).parent / "data" / "fissionDB" / "FPYParser.py"

        if parser_path.exists():
            print_info(f"Parsing JEFF files to {jeff_output}...")

            try:
                # Try to parse using the JEFF parser
                sys.path.insert(0, str(parser_path.parent))

                parsed_count = 0
                original_cwd = os.getcwd()
                os.chdir(jeff_output)

                for jeff_file in jeff_files:
                    try:
                        # Run parser on each file
                        result = subprocess.run([
                            sys.executable,
                            str(parser_path),
                            str(jeff_file)
                        ], capture_output=True, text=True, timeout=30)

                        if result.returncode == 0:
                            parsed_count += 1
                        else:
                            print_warning(f"Failed to parse {jeff_file.name}")

                    except Exception as e:
                        print_warning(f"Error parsing {jeff_file.name}: {e}")

                os.chdir(original_cwd)

                if parsed_count > 0:
                    print_success(f"Parsed {parsed_count}/{len(jeff_files)} JEFF files")
                    return True
                else:
                    print_error("No files were successfully parsed")
                    return False

            except Exception as e:
                os.chdir(original_cwd)
                print_error(f"JEFF parsing failed: {e}")
                print_info("You may need to parse files manually using FPYParser.py")
                return False
        else:
            print_warning("JEFF parser not found")
            print_info(f"Copy .jeff files manually to {jeff_output}")
            return True
    else:
        print_info("Skipping JEFF database setup")
        print_info("You can set it up later with: conflux-setup --setup-jeff")
        return True


def setup_jendl_database(db_path):
    """
    Setup JENDL-5 fission product yield database

    JENDL data must be downloaded manually from JAEA website
    """
    print_info("JENDL-5 database setup...")
    print_info("JENDL contains fission product yields from JAEA")
    print_warning("JENDL database must be downloaded manually from JAEA")
    print()

    print_info("To download JENDL-5:")
    print("  1. Visit: https://wwwndc.jaea.go.jp/jendl/jendl.html")
    print("  2. Navigate to JENDL-5 fission product yields")
    print("  3. Download the JENDL-5 fpy files")
    print("  4. Extract to a local directory")
    print()

    if get_yes_no("Do you have JENDL-5 files downloaded locally?", default=False):
        jendl_path = get_user_input("Enter path to JENDL-5 directory (containing fpy files)")
        if not jendl_path:
            print_info("Skipping JENDL database setup")
            return True

        jendl_path = Path(jendl_path).expanduser()

        if not jendl_path.exists():
            print_error(f"Directory not found: {jendl_path}")
            return False

        # Check for JENDL files (they use ENDF-6 format)
        jendl_files = list(jendl_path.glob("*.jendl")) + list(jendl_path.glob("*.txt"))

        if not jendl_files:
            print_error("No JENDL format files found in directory")
            print_info("Expected files with .jendl or .txt extensions")
            return False

        print_success(f"Found {len(jendl_files)} JENDL files")

        # JENDL uses ENDF-6 format, so we can use ENDF parser
        jendl_output = db_path / "fissionDB" / "JENDL"
        jendl_output.mkdir(parents=True, exist_ok=True)

        print_info(f"Parsing JENDL files to {jendl_output}...")
        print_info("JENDL uses ENDF-6 format - using ENDF parser")

        try:
            # Use inline ENDF parser since JENDL uses ENDF-6 format
            if not HAS_PARSER_DEPS:
                print_error("fortranformat not installed - cannot parse JENDL files")
                print_info("Install with: pip install fortranformat")
                return False

            parsed_count = 0
            original_cwd = os.getcwd()
            os.chdir(jendl_output)

            # Try to import unified FPY parser
            fpy_parser_path = Path(__file__).parent / "data" / "fissionDB" / "FPYParser.py"
            if fpy_parser_path.exists():
                sys.path.insert(0, str(fpy_parser_path.parent))

                try:
                    from FPYParser import ENDF6FPYreader

                    for jendl_file in jendl_files:
                        try:
                            ENDF6FPYreader(str(jendl_file))
                            parsed_count += 1
                        except Exception as e:
                            print_warning(f"Failed to parse {jendl_file.name}")

                except ImportError:
                    print_error("Could not import ENDF parser")
                    os.chdir(original_cwd)
                    return False

            os.chdir(original_cwd)

            if parsed_count > 0:
                print_success(f"Parsed {parsed_count}/{len(jendl_files)} JENDL files")
                return True
            else:
                print_error("No files were successfully parsed")
                return False

        except Exception as e:
            os.chdir(original_cwd)
            print_error(f"JENDL parsing failed: {e}")
            print_info("You may need to parse files manually")
            return False
    else:
        print_info("Skipping JENDL database setup")
        print_info("You can set it up later with: conflux-setup --setup-jendl")
        return True


def download_covariance_matrices(db_path):
    """
    Download covariance matrices from FYCoM project

    These matrices are too large to distribute with the package (~20-30 MB each).
    They are downloaded from the FYCoM GitHub repository.
    """
    print_info("Covariance matrix setup...")
    print_info("Covariance matrices from FYCoM (https://nucleardata.berkeley.edu/FYCoM/)")
    print_warning("These files are large (~300-500 MB total) and take time to download")
    print_info("Covariance matrices are optional but recommended for uncertainty calculations")
    print()

    if not get_yes_no("Download covariance matrices?", default=True):
        print_info("Skipping covariance matrix download")
        print_info("You can download them later with: conflux-setup --download-covariance")
        return True

    print_info("Downloading covariance and correlation matrices from FYCoM...")

    # GitHub raw URL base
    url_base = 'https://raw.githubusercontent.com/efmatthews/FYCoM/master/matrices/'

    # Database paths
    db_categories = ['/ENDF/', '/JEFF/']
    data_type = 'cumulative/'
    matrix_types = ['cov', 'corr', 'normed_cov']
    energies = ['T', 'F', 'H', 'SF']

    # Element mapping
    elements = {
        90: 'Th',
        92: 'U',
        94: 'Pu'
    }

    total_downloaded = 0
    total_skipped = 0
    total_failed = 0

    try:
        for z in elements:
            for mass in range(z*2+47, z*2+56):
                for energy in energies:
                    for mat_type in matrix_types:
                        for category in db_categories:
                            # Construct URL and local path
                            url_name = f"{url_base}{category[1:]}{data_type}{elements[z]}{mass}{energy}_cml_{mat_type}.csv"

                            # Local filename format
                            local_name = db_path / f"fissionDB{category}{mat_type}_nfy_{z}_{elements[z]}_{mass}_{energy}.csv"

                            # Skip if already exists
                            if local_name.exists():
                                total_skipped += 1
                                continue

                            # Ensure directory exists
                            local_name.parent.mkdir(parents=True, exist_ok=True)

                            try:
                                # Try to download
                                with urllib.request.urlopen(url_name, timeout=30) as response:
                                    if response.status == 200:
                                        data = response.read()
                                        with open(local_name, 'wb') as f:
                                            f.write(data)
                                        total_downloaded += 1
                                        size_mb = len(data) / 1024 / 1024
                                        print(f"  Downloaded: {local_name.name} ({size_mb:.1f} MB)")

                            except urllib.error.HTTPError as err:
                                if err.code == 404:
                                    # File doesn't exist on server - this is expected for some combinations
                                    pass
                                else:
                                    total_failed += 1
                                    print_warning(f"Failed to download {local_name.name}: HTTP {err.code}")

                            except Exception as e:
                                total_failed += 1
                                print_warning(f"Failed to download {local_name.name}: {e}")

        # Summary
        print()
        if total_downloaded > 0:
            print_success(f"Downloaded {total_downloaded} covariance matrix files")
        if total_skipped > 0:
            print_info(f"Skipped {total_skipped} existing files")
        if total_failed > 0:
            print_warning(f"Failed to download {total_failed} files")

        if total_downloaded > 0 or total_skipped > 0:
            print_success("Covariance matrices setup complete")
            return True
        else:
            print_error("No covariance matrices were downloaded")
            return False

    except Exception as e:
        print_error(f"Failed to download covariance matrices: {e}")
        return False

def verify_installation():
    """Verify CONFLUX installation"""
    print_info("Verifying installation...")

    checks = []

    # Check Python version
    py_version = sys.version_info
    if py_version >= (3, 6):
        print_success(f"Python version: {py_version.major}.{py_version.minor}")
        checks.append(True)
    else:
        print_error(f"Python version {py_version.major}.{py_version.minor} < 3.6")
        checks.append(False)

    # Check if conflux is importable
    try:
        import conflux
        print_success(f"CONFLUX package found: {conflux.__file__}")
        checks.append(True)
    except ImportError as e:
        print_error(f"Cannot import conflux: {e}")
        checks.append(False)
        return False

    # Check required dependencies
    # Note: Some packages may have numpy compatibility issues during import
    # but are still installed. Check via pip/conda first.
    required_packages = ['numpy', 'scipy', 'tqdm', 'matplotlib', 'iminuit', 'fortranformat', 'pandas', 'xraydb']

    # Try to check via pip list (more reliable than importing with numpy conflicts)
    try:
        result = subprocess.run(
            [sys.executable, '-m', 'pip', 'list', '--format=json'],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode == 0:
            import json
            installed = {pkg['name'].lower().replace('-', '_'): pkg['version']
                        for pkg in json.loads(result.stdout)}

            for package in required_packages:
                pkg_key = package.lower().replace('-', '_')
                if pkg_key in installed:
                    print_success(f"Found dependency: {package} ({installed[pkg_key]})")
                    checks.append(True)
                else:
                    print_error(f"Missing dependency: {package}")
                    checks.append(False)
        else:
            # Fallback to import check
            for package in required_packages:
                try:
                    __import__(package)
                    print_success(f"Found dependency: {package}")
                    checks.append(True)
                except ImportError:
                    print_error(f"Missing dependency: {package}")
                    checks.append(False)

    except Exception as e:
        # Fallback to import check if pip check fails
        print_warning(f"Could not check via pip: {e}")
        for package in required_packages:
            try:
                __import__(package)
                print_success(f"Found dependency: {package}")
                checks.append(True)
            except ImportError as err:
                print_error(f"Missing dependency: {package}")
                checks.append(False)

    # Check CONFLUX_DB
    db_path = os.environ.get('CONFLUX_DB')
    if db_path:
        print_success(f"CONFLUX_DB set: {db_path}")
        if Path(db_path).exists():
            print_success(f"Database directory exists")
            checks.append(True)
        else:
            print_warning(f"Database directory does not exist yet")
            checks.append(True)
    else:
        print_warning("CONFLUX_DB not set (optional)")
        checks.append(True)

    # Check database files
    if db_path:
        db_path = Path(db_path)
        endf_dir = db_path / "fissionDB" / "ENDF"
        if endf_dir.exists():
            xml_files = list(endf_dir.glob("*.xml"))
            if xml_files:
                print_success(f"Found {len(xml_files)} ENDF database files")
                checks.append(True)
            else:
                print_warning("No ENDF database files found")
                checks.append(False)
        else:
            print_warning("ENDF database not installed")
            checks.append(False)

    print()
    if all(checks):
        print_success("All checks passed!")
        return True
    else:
        print_warning("Some checks failed")
        return False

def interactive_setup():
    """Run interactive setup wizard"""
    print_header("CONFLUX Installation Setup")

    print("This wizard will help you set up CONFLUX:")
    print("  1. Configure database path (CONFLUX_DB)")
    print("  2. Download ENDF fission product yield database")
    print("  3. Set up ENSDF decay database (optional)")
    print("  4. Set up JEFF/JENDL fission databases (optional)")
    print("  5. Download covariance matrices from FYCoM (optional)")
    print("  6. Verify installation")
    print()

    # Step 1: Database path
    print_header("Step 1: Database Path Configuration")

    default_path = get_default_db_path()

    # Check if default is the installation's data directory
    import conflux
    install_data_dir = None
    if conflux.__file__:
        install_data_dir = Path(conflux.__file__).parent / "data"
        if install_data_dir.exists() and default_path == install_data_dir:
            print_success(f"Using installation's data directory: {default_path}")
            print_info("(Editable install detected - no copying needed)")
        else:
            print_info(f"Default database path: {default_path}")
    else:
        print_info(f"Default database path: {default_path}")

    if get_yes_no("Use default path?", default=True):
        db_path = default_path
    else:
        custom_path = get_user_input("Enter custom database path", str(default_path))
        db_path = Path(custom_path).expanduser()

    # Create directory if it doesn't exist
    db_path.mkdir(parents=True, exist_ok=True)
    print_success(f"Database path: {db_path}")

    # Setup environment variable
    if get_yes_no("Set CONFLUX_DB environment variable?", default=True):
        setup_environment_variable(db_path)

    # Copy package data only if not using installation's data directory
    copy_needed = False
    if install_data_dir and db_path == install_data_dir:
        print_info("Using installation's data directly - no copying needed")
        copy_needed = False
    else:
        print_header("Step 1.5: Copy Package Data")
        print_info("Copying all git-tracked data as fallback...")
        print_info("(This ensures basic functionality even if downloads fail)")
        copy_success = copy_package_data_to_db(db_path)
        copy_needed = True

        if not copy_success:
            print_warning("Data copying failed - some files may be missing")
            print_info("You may need to download databases manually")

    # Step 2: ENDF database
    print_header("Step 2: ENDF Database Download")

    if get_yes_no("Download ENDF-B-VIII.0 database? (~3 MB)", default=True):
        endf_success = download_endf_database(db_path)
        if not endf_success and copy_needed:
            print_info("Download failed - using data from installation as fallback")
    elif copy_needed and install_data_dir and (install_data_dir / "fissionDB" / "ENDF").exists():
        print_info("Skipped download - using data from installation")

    # Step 3: ENSDF database
    print_header("Step 3: ENSDF Database (Optional)")

    if get_yes_no("Set up ENSDF database?", default=False):
        ensdf_success = download_ensdf_database(db_path)
        if not ensdf_success and copy_needed:
            print_info("ENSDF setup failed - using beta databases from installation as fallback")
    elif copy_needed and install_data_dir and (install_data_dir / "betaDB").exists():
        print_info("Skipped ENSDF - using beta databases from installation")

    # Step 4: JEFF and JENDL databases
    print_header("Step 4: JEFF/JENDL Databases (Optional)")

    print_info("JEFF and JENDL provide alternative fission product yield data")
    print_info("These databases can now be downloaded automatically!")
    print()

    if get_yes_no("Download JEFF-3.3 database?", default=False):
        print_info("Attempting automatic download from NEA Data Bank...")
        success = download_jeff_database(db_path)

        if not success:
            print_warning("Automatic download failed")
            if get_yes_no("Try manual setup instead?", default=True):
                setup_jeff_database(db_path)

    if get_yes_no("Download JENDL-5 database?", default=False):
        print_info("Attempting automatic download from JAEA...")
        print_warning("Note: JAEA uses self-signed SSL certificates (security workaround applied)")
        success = download_jendl_database(db_path)

        if not success:
            print_warning("Automatic download failed")
            if get_yes_no("Try manual setup instead?", default=True):
                setup_jendl_database(db_path)

    # Step 5: Covariance matrices
    print_header("Step 5: Covariance Matrices (Optional)")

    if get_yes_no("Download covariance matrices from FYCoM?", default=True):
        download_covariance_matrices(db_path)

    # Step 6: Verification
    print_header("Step 6: Verification")

    # Check if data was copied/available
    if copy_needed and db_path.exists():
        print_info("Checking CONFLUX_DB contents...")
        beta_db = db_path / "betaDB"
        fission_db = db_path / "fissionDB" / "ENDF"

        beta_xmls = len(list(beta_db.glob("*.xml"))) if beta_db.exists() else 0
        fission_xmls = len(list(fission_db.glob("*.xml"))) if fission_db.exists() else 0

        if beta_xmls > 0 or fission_xmls > 0:
            print_success(f"✓ Found {beta_xmls} beta databases and {fission_xmls} fission yield files")
            print_success("✓ Data successfully copied to CONFLUX_DB")
        else:
            print_warning("No database files found in CONFLUX_DB")
            print_info("You may need to download databases manually")

    verify_installation()

    # Final message
    print_header("Setup Complete!")
    print_success("CONFLUX is ready to use!")
    print()
    print("Next steps:")
    print(f"  1. Restart your terminal or run: source ~/.zshrc")
    print(f"  2. Test installation: python -c 'import conflux; print(conflux.__version__)'")
    print(f"  3. Read documentation: {Path(__file__).parent.parent / 'README.md'}")
    print()

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="CONFLUX Setup and Configuration Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run interactive setup wizard
  conflux-setup

  # Verify installation only
  conflux-setup --verify

  # Set database path only
  conflux-setup --set-db-path /path/to/database

  # Download ENDF database only
  conflux-setup --download-endf --db-path /path/to/database

  # Download ENSDF database from NNDC
  conflux-setup --download-ensdf --db-path /path/to/database

  # Parse your own ENSDF to betaDB and update config.py
  conflux-setup --parse-ensdf /path/to/ensdf_270101 --db-path /path/to/database
  conflux-setup --parse-ensdf /path/to/ensdf --version-id 270101 --db-path /path/to/database

  # Download covariance matrices from FYCoM
  conflux-setup --download-covariance --db-path /path/to/database

  # Download JEFF or JENDL databases automatically
  conflux-setup --download-jeff --db-path /path/to/database
  conflux-setup --download-jendl --db-path /path/to/database

  # Setup JEFF/JENDL from manually downloaded files (fallback)
  conflux-setup --setup-jeff --db-path /path/to/database
  conflux-setup --setup-jendl --db-path /path/to/database
        """
    )

    parser.add_argument(
        "--interactive",
        action="store_true",
        default=True,
        help="Run interactive setup wizard (default)"
    )

    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify installation only"
    )

    parser.add_argument(
        "--set-db-path",
        metavar="PATH",
        help="Set CONFLUX_DB path and exit"
    )

    parser.add_argument(
        "--download-endf",
        action="store_true",
        help="Download ENDF database"
    )

    parser.add_argument(
        "--download-ensdf",
        action="store_true",
        help="Download ENSDF database from NNDC"
    )

    parser.add_argument(
        "--download-covariance",
        action="store_true",
        help="Download covariance matrices from FYCoM"
    )

    parser.add_argument(
        "--download-jeff",
        action="store_true",
        help="Download JEFF-3.3 database automatically from NEA"
    )

    parser.add_argument(
        "--download-jendl",
        action="store_true",
        help="Download JENDL-5 database automatically from JAEA (uses SSL workaround)"
    )

    parser.add_argument(
        "--setup-jeff",
        action="store_true",
        help="Setup JEFF-3.3 database from manually downloaded files"
    )

    parser.add_argument(
        "--setup-jendl",
        action="store_true",
        help="Setup JENDL-5 database from manually downloaded files"
    )

    parser.add_argument(
        "--parse-ensdf",
        metavar="ENSDF_DIR",
        help="Parse ENSDF directory to betaDB and update config.py"
    )

    parser.add_argument(
        "--version-id",
        metavar="YYMMDD",
        help="Version identifier for parsed database (default: auto-detect from directory name)"
    )

    parser.add_argument(
        "--db-path",
        metavar="PATH",
        help="Database path (for all database setup options)"
    )

    args = parser.parse_args()

    # Handle non-interactive modes
    if args.verify:
        verify_installation()
        return

    if args.set_db_path:
        db_path = Path(args.set_db_path).expanduser()
        db_path.mkdir(parents=True, exist_ok=True)
        setup_environment_variable(db_path)
        return

    if args.download_endf:
        if not args.db_path:
            print_error("--db-path required with --download-endf")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        download_endf_database(db_path)
        return

    if args.download_ensdf:
        if not args.db_path:
            print_error("--db-path required with --download-ensdf")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        download_ensdf_database(db_path)
        return

    if args.download_covariance:
        if not args.db_path:
            print_error("--db-path required with --download-covariance")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        download_covariance_matrices(db_path)
        return

    if args.download_jeff:
        if not args.db_path:
            print_error("--db-path required with --download-jeff")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        success = download_jeff_database(db_path)
        if not success:
            print_warning("Automatic download failed. Use --setup-jeff for manual setup.")
            sys.exit(1)
        return

    if args.download_jendl:
        if not args.db_path:
            print_error("--db-path required with --download-jendl")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        success = download_jendl_database(db_path)
        if not success:
            print_warning("Automatic download failed. Use --setup-jendl for manual setup.")
            sys.exit(1)
        return

    if args.setup_jeff:
        if not args.db_path:
            print_error("--db-path required with --setup-jeff")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        setup_jeff_database(db_path)
        return

    if args.setup_jendl:
        if not args.db_path:
            print_error("--db-path required with --setup-jendl")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        setup_jendl_database(db_path)
        return

    if args.parse_ensdf:
        if not args.db_path:
            print_error("--db-path required with --parse-ensdf")
            sys.exit(1)
        db_path = Path(args.db_path).expanduser()
        ensdf_dir = Path(args.parse_ensdf).expanduser()

        if not ensdf_dir.exists():
            print_error(f"ENSDF directory not found: {ensdf_dir}")
            sys.exit(1)

        version_id = args.version_id if args.version_id else None
        success = parse_ensdf_to_betadb(ensdf_dir, db_path, version_id)
        sys.exit(0 if success else 1)

    # Run interactive setup by default
    interactive_setup()

if __name__ == "__main__":
    main()
