#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ENDF Database Downloader and Parser

This script downloads the latest ENDF-B-VIII database from NNDC
and parses the fission product yield (FPY) files into XML format.

Author: Xianyi Zhang, LLNL
Date: 2026-07-15
"""

import os
import sys
import argparse
import urllib.request
import zipfile
import shutil
from pathlib import Path
from xml.dom import minidom
import fortranformat as ff

# ENDF database URL
# Note: VIII.1 is distributed in GNDS format, not traditional ENDF-6
# and must be obtained separately from NNDC
ENDF_BASE_URL = "https://www.nndc.bnl.gov/endf-b8.0/zips/"
ENDF_VERSIONS = {
    "VIII.0": "ENDF-B-VIII.0_nfy.zip"
}

class XMLedit:
    """Class to create and write XML files for fission product yield data"""

    def __init__(self, DBname):
        self.root = minidom.Document()
        self.DB = self.root.createElement(DBname)
        self.root.appendChild(self.DB)
        self.outputName = DBname

    def createHEAD(self, ZA, AWR, LE, MT):
        """
        Create HEAD element with fission isotope metadata

        Args:
            ZA: Z*1000+A number of fission isotope
            AWR: Atomic weight ratio
            LE: Number of fission triggering energies
            MT: ENDF data branch (IFP or CFP)
        """
        self.HEAD = self.root.createElement('HEAD')
        self.HEAD.setAttribute('FissionZA', str(ZA))
        self.HEAD.setAttribute('AWR', str(AWR))
        self.HEAD.setAttribute('LE', str(LE))
        self.HEAD.setAttribute('MT', str(MT))
        self.DB.appendChild(self.HEAD)

    def createLIST(self, Ei, I, NFP):
        """
        Create LIST element for a specific fission energy

        Args:
            Ei: Fission triggering energy
            I: Interpolation scheme
            NFP: Total number of fission products
        """
        self.LIST = self.root.createElement('LIST')
        self.LIST.setAttribute('Ei', str(Ei))
        self.LIST.setAttribute('Ii', str(I))
        self.LIST.setAttribute('NFPi', str(NFP))
        self.HEAD.appendChild(self.LIST)

    def editCONT(self, data):
        """
        Add CONT element with fission product data

        Args:
            data: [ZA, FPS, Y, DY] - isotope data, isomeric state, yield, uncertainty
        """
        CONT = self.root.createElement('CONT')
        CONT.setAttribute('ZA', str(data[0]))
        CONT.setAttribute('FPS', str(data[1]))
        CONT.setAttribute('Y', str(data[2]))
        CONT.setAttribute('DY', str(data[3]))
        self.LIST.appendChild(CONT)

    def saveXML(self, output_dir=None):
        """Save XML to file"""
        xml_str = self.root.toprettyxml(indent="\t")
        output_path = self.outputName + ".xml"
        if output_dir:
            output_path = os.path.join(output_dir, os.path.basename(output_path))
        with open(output_path, "w") as f:
            f.write(xml_str)
        return output_path


def download_endf_database(version="VIII.0", download_dir="endf_downloads"):
    """
    Download ENDF database from NNDC

    Args:
        version: ENDF version to download (VIII.0 or VIII.1)
        download_dir: Directory to save downloaded files

    Returns:
        Path to extracted nfy directory
    """
    if version not in ENDF_VERSIONS:
        raise ValueError(f"Version {version} not supported. Choose from: {list(ENDF_VERSIONS.keys())}")

    # Create download directory
    download_path = Path(download_dir)
    download_path.mkdir(parents=True, exist_ok=True)

    zip_filename = ENDF_VERSIONS[version]
    zip_path = download_path / zip_filename
    url = ENDF_BASE_URL + zip_filename

    print(f"Downloading ENDF-B-{version} from {url}...")

    try:
        # Download with progress
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            sys.stdout.write(f"\rProgress: {percent:.1f}% ({downloaded / 1024 / 1024:.1f} MB)")
            sys.stdout.flush()

        urllib.request.urlretrieve(url, zip_path, reporthook=report_progress)
        print("\nDownload complete!")

    except urllib.error.URLError as e:
        print(f"\nError downloading database: {e}")
        print("Please check your internet connection or try a different version.")
        sys.exit(1)

    # Extract zip file
    print(f"Extracting {zip_filename}...")
    extract_dir = download_path / f"ENDF-B-{version}"

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    print(f"Extraction complete!")

    # Find nfy directory - the zip extracts to ENDF-B-VIII.0_nfy/
    # Search for directory containing .endf files
    nfy_dir = None
    for root, dirs, files in os.walk(extract_dir):
        # Check if this directory contains .endf files
        endf_files = [f for f in files if f.endswith('.endf')]
        if endf_files:
            nfy_dir = Path(root)
            break

    if not nfy_dir:
        raise FileNotFoundError(f"Could not find directory with .endf files in extracted archive")

    return nfy_dir


def parse_endf_file(filename, output_dir=None):
    """
    Parse ENDF fission product yield file and create XML

    Args:
        filename: Path to .endf file
        output_dir: Directory to save XML output

    Returns:
        Path to created XML file
    """
    record_reader = ff.FortranRecordReader('(A66,I4,I2,I3,I5)')
    DBname = os.path.basename(filename).rsplit('.', 1)[0]
    xmloutput = XMLedit(DBname)

    with open(filename) as inputfile:
        linenum = 0
        ZA = -1
        LEplus1 = -1
        Ei = []
        Ii = []
        NNi = []
        NFPi = []
        i = 0
        CONT = []

        for line in inputfile:
            TEXT, MAT, MF, MT, NS = record_reader.read(line)

            # Independent fission products
            if MT == 454:
                data_reader = ff.FortranRecordReader('(6E11.0)')
                data = data_reader.read(TEXT)

                if linenum == 0:
                    ZA = int(data[0])
                    AWR = float(data[1])
                    LEplus1 = int(data[2])
                    Ei = [0] * LEplus1
                    Ii = [0] * LEplus1
                    NNi = [0] * LEplus1
                    NFPi = [0] * LEplus1
                    xmloutput.createHEAD(ZA, AWR, LEplus1, 'IFP')
                else:
                    # Reading CONT data
                    if len(CONT) <= NFPi[i] and NFPi[i] > 0:
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 4:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 2:
                            CONT[-1].append(data[0])
                            CONT[-1].append(data[1])
                            xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                        if len(CONT[-1]) != 4:
                            continue

                    if len(CONT) == NFPi[i] > 0:
                        CONT = []
                        i += 1
                        continue

                    # Reading LIST data
                    if i < LEplus1 and len(CONT) == NFPi[i] == 0:
                        Ei[i] = data[0]
                        Ii[i] = int(data[2])
                        NNi[i] = int(data[4])
                        NFPi[i] = int(data[5])
                        xmloutput.createLIST(Ei[i], Ii[i], NFPi[i])

                linenum += 1

            # Data type separator
            if linenum != 0 and MT == 0:
                linenum = 0

            # Cumulative fission products
            if MT == 459:
                data = data_reader.read(TEXT)
                if linenum == 0:
                    i = 0
                    CONT = []
                    ZA = int(data[0])
                    AWR = float(data[1])
                    LEplus1 = int(data[2])
                    Ei = [0] * LEplus1
                    Ii = [0] * LEplus1
                    NNi = [0] * LEplus1
                    NFPi = [0] * LEplus1
                    xmloutput.createHEAD(ZA, AWR, LEplus1, 'CFP')
                else:
                    # Reading CONT
                    if len(CONT) <= NFPi[i] and NFPi[i] > 0:
                        if len(CONT) == 0:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 4:
                            datacache = [int(data[0]), float(data[1]), float(data[2]), float(data[3])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                        elif len(CONT[-1]) == 2:
                            CONT[-1].append(data[0])
                            CONT[-1].append(data[1])
                            xmloutput.editCONT(CONT[-1])
                            datacache = [int(data[2]), float(data[3]), float(data[4]), float(data[5])]
                            if datacache[0] != 0:
                                CONT.append(datacache)
                                xmloutput.editCONT(CONT[-1])
                        if len(CONT[-1]) != 4:
                            continue

                    if len(CONT) == NFPi[i] > 0:
                        CONT = []
                        i += 1
                        continue

                    # Reading LIST
                    if i < LEplus1 and len(CONT) == NFPi[i] == 0:
                        Ei[i] = data[0]
                        Ii[i] = int(data[2])
                        NNi[i] = int(data[4])
                        NFPi[i] = int(data[5])
                        xmloutput.createLIST(Ei[i], Ii[i], NFPi[i])

                linenum += 1

    return xmloutput.saveXML(output_dir)


def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(
        description="Download and parse ENDF fission product yield database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download ENDF-B-VIII.0 and parse to current directory
  %(prog)s

  # Parse existing ENDF-B-VIII.0 files without downloading
  %(prog)s --no-download --input-dir ./ENDF-B-VIII.0/nfy --output-dir ./parsed_xml

  # Parse ENDF-B-VIII.1 GNDS files (must have local copy)
  # Note: VIII.1 is distributed in GNDS XML format, different from VIII.0
  %(prog)s --no-download --input-dir ./ENDF-B-VIII.1-GNDS/nfy --output-dir ./parsed_xml
        """
    )

    parser.add_argument(
        "--version",
        default="VIII.0",
        choices=["VIII.0"],
        help="ENDF database version to download (default: VIII.0). Note: VIII.1 uses GNDS format and must be parsed from local files using --no-download"
    )

    parser.add_argument(
        "--download-dir",
        default="endf_downloads",
        help="Directory to save downloaded files (default: endf_downloads)"
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory to save parsed XML files (default: same as script location)"
    )

    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip download and parse existing files from --input-dir"
    )

    parser.add_argument(
        "--input-dir",
        default=None,
        help="Directory containing .endf files to parse (used with --no-download)"
    )

    parser.add_argument(
        "--keep-downloads",
        action="store_true",
        help="Keep downloaded files after parsing (default: delete)"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.no_download and not args.input_dir:
        parser.error("--input-dir is required when using --no-download")

    # Set output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path.cwd()

    print("=" * 70)
    print("ENDF Fission Product Yield Database Parser")
    print("=" * 70)

    # Download or use existing files
    if args.no_download:
        nfy_dir = Path(args.input_dir)
        if not nfy_dir.exists():
            print(f"Error: Input directory '{nfy_dir}' does not exist")
            sys.exit(1)
        print(f"\nUsing existing ENDF files from: {nfy_dir}")
    else:
        nfy_dir = download_endf_database(args.version, args.download_dir)
        print(f"ENDF files location: {nfy_dir}")

    # Parse all .endf files
    print(f"\nParsing ENDF files to XML...")
    print(f"Output directory: {output_dir}")
    print("-" * 70)

    endf_files = sorted(nfy_dir.glob("*.endf"))

    if not endf_files:
        print(f"No .endf files found in {nfy_dir}")
        sys.exit(1)

    parsed_count = 0
    failed_files = []

    for endf_file in endf_files:
        try:
            print(f"Parsing: {endf_file.name}...", end=" ")
            xml_path = parse_endf_file(str(endf_file), str(output_dir))
            print(f"✓ Created: {os.path.basename(xml_path)}")
            parsed_count += 1
        except Exception as e:
            print(f"✗ Failed: {e}")
            failed_files.append((endf_file.name, str(e)))

    # Summary
    print("-" * 70)
    print(f"\nParsing complete!")
    print(f"  Successfully parsed: {parsed_count} files")

    if failed_files:
        print(f"  Failed: {len(failed_files)} files")
        for filename, error in failed_files:
            print(f"    - {filename}: {error}")

    # Cleanup
    if not args.no_download and not args.keep_downloads:
        print(f"\nCleaning up downloads from {args.download_dir}...")
        try:
            shutil.rmtree(args.download_dir)
            print("Download directory removed.")
        except Exception as e:
            print(f"Warning: Could not remove download directory: {e}")

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
