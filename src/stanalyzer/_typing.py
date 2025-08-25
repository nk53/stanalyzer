from collections.abc import Callable, Sized  # noqa: F401
from typing import (  # noqa: F401
    TYPE_CHECKING,
    Annotated,
    Any,
    Generic,
    Literal,
    Optional,
    TypeAlias,
    TypeVar,
)

# the 'type' keyword was added in 3.11; stanalyzer aims to support 3.10
StrDict: TypeAlias = dict[str, Any]
StrDictList: TypeAlias = list[dict[str, Any]]

# error types; no idea why pydantic doesn't have these out of the box
CustomError: TypeAlias = Literal[
    'dir_not_creatable',
    'dir_not_writable',
    'file_not_found',
    'glob_failed',
    'is_a_directory',
    'not_a_regular_file',
    'not_a_directory',
    'not_abspath',
    'field_missing',
]

error_templates: dict[CustomError, str] = {
    'dir_not_creatable': "No permission to create '{}'",
    'dir_not_writable': "No permission to write to '{}'",
    'file_not_found': "No such file: '{}'",
    'glob_failed': "No files matched by pattern: '{}'",
    'is_a_directory': "'{}' is not a regular file",
    'not_a_regular_file': "'{}' is not a regular file",
    'not_a_directory': "'{}' is not a directory",
    'not_abspath': "'{}' is not an absolute path",
    'field_missing': "This field is required",
}
