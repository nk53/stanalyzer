__version__ = "0.1.0"

from .stanalyzer import get_settings
from .stanalyzer import main as analyze

__all__ = ["analyze", "config", "get_settings"]
