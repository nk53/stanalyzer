import argparse

__all__ = ["analyze", "config", "get_settings", "FakeParser"]


def __getattr__(name):
    """Lazily re-export analyze/get_settings; importing cli must not execute stanalyzer.py."""
    if name in ("analyze", "get_settings"):
        from importlib import import_module
        return getattr(import_module(f"{__name__}.stanalyzer"),
                       "main" if name == "analyze" else name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


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
