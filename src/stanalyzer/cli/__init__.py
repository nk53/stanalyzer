import argparse
from .stanalyzer import get_settings
from .stanalyzer import main as analyze

__all__ = ["analyze", "config", "get_settings", "FakeParser"]


class FakeParser(argparse.ArgumentParser):
    """An ArgumentParser that ignores all arguments."""

    def __init__(self, *args, **kwargs):
        pass

    def parse_args(self, args, **kwargs):  # type: ignore[override]
        """Ignores arguments and returns an empty namespace.

        Note that get_settings() will still read project.json (if present)
        and append its settings to the namespace.
        """
        return argparse.Namespace()
