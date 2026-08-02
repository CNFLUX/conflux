# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

import os

# Database location configuration
# CONFLUX_DB environment variable is OPTIONAL:
#   - If SET: Use the specified custom database location
#   - If NOT SET: Auto-detect from package installation
#
# This allows users to:
#   1. Use the package immediately after install (no setup needed)
#   2. Optionally point to custom/updated databases
#   3. Test different database versions easily
#
# To use a custom location:
#   export CONFLUX_DB="/path/to/your/databases"
#
# To check current location:
#   python -c "from conflux.config import CONFLUX_DB; print(CONFLUX_DB)"
try:
    CONFLUX_DB = os.environ["CONFLUX_DB"]
except KeyError:
    # CONFLUX_DB not set - use installation's data directory as default
    import sys
    from pathlib import Path
    if hasattr(sys.modules[__name__], '__file__') and __file__:
        CONFLUX_DB = str(Path(__file__).parent / "data")
    else:
        CONFLUX_DB = str(Path.home() / ".conflux" / "data")

# Database version configuration
# Update BETA_DB_VERSION when switching to a new database version
BETA_DB_VERSION = "260707"

# Auto-generated filenames based on version
BETA_DB_FILENAME = f"ENSDF_betaDB_{BETA_DB_VERSION}.xml"
BETA_DB_EC_FILENAME = f"ENSDF_betaDB_EC_{BETA_DB_VERSION}.xml"

# Full paths to default database files
BETA_DB_PATH = os.path.join(CONFLUX_DB, "betaDB", BETA_DB_FILENAME)
BETA_DB_EC_PATH = os.path.join(CONFLUX_DB, "betaDB", BETA_DB_EC_FILENAME)

# Optional: Named database registry for managing multiple databases
# Uncomment and customize if you need to work with multiple database versions simultaneously
"""
NAMED_DATABASES = {
    "ensdf_260707": {
        "beta": "ENSDF_betaDB_260707.xml",
        "ec": "ENSDF_betaDB_EC_260707.xml",
        "description": "Official ENSDF July 2026"
    },
    "xundl_260101": {
        "beta": "XUNDL_betaDB_260101.xml",
        "ec": "XUNDL_betaDB_EC_260101.xml",
        "description": "XUNDL January 2026"
    },
}

def get_database_path(db_name="ensdf_260707", db_type="beta"):
    '''Get path to a named database.

    Args:
        db_name: Name from NAMED_DATABASES dictionary
        db_type: Either "beta" (B-) or "ec" (EC/B+)

    Returns:
        Full path to database file

    Example:
        from conflux.config import get_database_path
        from conflux.BetaEngine import BetaEngine

        # Use XUNDL database instead of default
        engine = BetaEngine(targetDB=get_database_path("xundl_260101", "beta"))
    '''
    if db_name not in NAMED_DATABASES:
        raise ValueError(f"Unknown database: {db_name}. Available: {list(NAMED_DATABASES.keys())}")
    if db_type not in ["beta", "ec"]:
        raise ValueError(f"Invalid db_type: {db_type}. Must be 'beta' or 'ec'")

    filename = NAMED_DATABASES[db_name][db_type]
    return os.path.join(CONFLUX_DB, "betaDB", filename)
"""
