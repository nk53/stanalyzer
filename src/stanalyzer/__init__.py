try:
    from importlib.metadata import version
    __version__ = version(__package__)
    _release = bool(__version__ != '0.0.0' and 'dev' not in __version__)
except Exception:
    __version__ = '0.0.0'
    _release = False
