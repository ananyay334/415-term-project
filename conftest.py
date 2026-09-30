# piv_pressure/conftest.py
"""pytest configuration: make the project root importable (fields, core, ...)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
