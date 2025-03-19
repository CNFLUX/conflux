# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

import os

# Load environment variables and provide default values if needed
CONFLUX_DB = os.getenv("CONFLUX_DB", os.environ["CONFLUX_DB"])
