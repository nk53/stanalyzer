import shutil
import typing as t
from pathlib import Path

__all__ = ['p_int', 'nn_int', 'p_float', 'nn_float']
T = t.TypeVar('T')


#
# Numeric validators
#

def p_int(value: t.Any) -> int:
    """Checks if value is a positive integer.

    Returns valid int or raises TypeError.
    """
    value = int(value)
    if value <= 0:
        raise TypeError(f"Not a positive integer: '{value}'")
    return value


def nn_int(value: t.Any) -> int:
    """Checks if value is a non-negative integer.

    Returns valid int or raises TypeError.
    """
    value = int(value)
    if value < 0:
        raise TypeError(f"Not a non-negative integer: '{value}'")
    return value


def p_float(value: t.Any) -> float:
    """Checks if value is a positive float.

    Returns valid float or raises TypeError.
    """
    value = float(value)
    if value <= 0:
        raise TypeError(f"Not a positive float: '{value}'")
    return value


def nn_float(value: t.Any) -> float:
    """Checks if value is a non-negative float.

    Returns valid float or raises TypeError.
    """
    value = float(value)
    if value < 0:
        raise TypeError(f"Not a non-negative float: '{value}'")
    return value


#
# Misc
#

def exec_name(value: t.Any, error=True) -> str | None:
    """Searches PATH for exec. The path is returned if found; else None.

    If error is True and the path was not found, raises a ValueError.
    """
    exec_path = Path(value)
    result: str | None

    match len(exec_path.parts):
        case 0:
            result = None
        case 1:
            result = shutil.which(value)
        case _:
            result = value

    if error and result is None:
        if value and exec_path.exists():
            raise ValueError(f"Not executable: {value}")
        if value:
            raise ValueError(f"No such file: {value}")
        raise ValueError(f"Not a path: {value}")

    return result
