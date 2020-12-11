"""Configuration variables including stnadard data directory"""
from pathlib import Path

STANDARD_DATA_DIRECTORY = Path.home() / "ixdat"

if not STANDARD_DATA_DIRECTORY.exists():
    STANDARD_DATA_DIRECTORY.mkdir()
