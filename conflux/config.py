import os

# Load environment variables and provide default values if needed
CONFLUX_DB = os.getenv("CONFLUX_DB", os.environ["CONFLUX_DB"])
