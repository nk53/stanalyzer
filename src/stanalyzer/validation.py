import os
import sys
from pathlib import Path
from ._typing import (
    Any,
    Annotated,
    CustomError,
    Literal,
    Optional,
    TypeVar,
    TYPE_CHECKING,
    error_templates,
)
from typing_extensions import Self

from pydantic import (
    BaseModel,
    EmailStr,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)
from pydantic.functional_validators import (
    AfterValidator as AV,
    BeforeValidator as BV,
)
from pydantic_core import PydanticCustomError
import stanalyzer

T = TypeVar('T', bound=type)
BM = TypeVar('BM', bound=BaseModel)


def maybe_raise(condition: object, error_type: CustomError,
                *template_vars: Any, context: dict[str, Any] | None = None) -> None:
    if condition:
        return
    if TYPE_CHECKING:
        assert callable(PydanticCustomError)

    ctx = context if context is not None else {}
    tpl = error_templates[error_type]
    raise PydanticCustomError(error_type, tpl.format(*template_vars, **ctx))


def exists(p: Path | str) -> Path:
    p = Path(p)
    maybe_raise(p.exists(), 'file_not_found', p)
    return p


def is_file(p: Path | str) -> Path:
    p = Path(p)
    maybe_raise(p.is_file(), 'not_a_regular_file', p)
    return p


def is_dir(p: Path | str) -> Path:
    p = Path(p)
    maybe_raise(p.is_dir(), 'not_a_directory', p)
    return p


def path_is_absolute(p: Path | str) -> Path:
    p = Path(p)
    maybe_raise(p.is_absolute(), 'not_abspath', p)
    return p


def not_empty(v: T, info: ValidationInfo) -> T:
    maybe_raise(v, 'field_missing')
    return v


def dir_is_writable(p: Path | str) -> Path:
    """Check that we can write to OR create a dir"""
    p = Path(p)
    if p.is_dir():
        # , f"No permission to write to {p}"
        maybe_raise(os.access(p, os.W_OK), 'dir_not_writable', p)
    else:
        # this approach doesn't require abspath args
        j = p.resolve()
        while (j != j.parent and not j.exists()):
            j = j.parent

        maybe_raise(os.access(j, os.W_OK), 'dir_not_creatable', p)
    return p


def str_to_args(cls: type[BM]) -> Any:
    fields = cls.model_fields
    keys = tuple(fields.keys())

    def _str_to_args(value: Any) -> Any:
        if isinstance(value, str) and len(keys) > 1:
            model_dict = dict(zip(keys, value.split()))
            return cls(**model_dict)

        return value

    return _str_to_args


ExistingFile = Annotated[Path, AV(not_empty), AV(exists), AV(is_file)]
ExistingDir = Annotated[Path, AV(not_empty), AV(exists), AV(is_dir)]
WritableDir = Annotated[Path, AV(not_empty), AV(dir_is_writable)]
ExistingFileAbs = Annotated[Path, AV(not_empty), AV(path_is_absolute), AV(exists), AV(is_file)]
ExistingDirAbs = Annotated[Path, AV(not_empty), AV(path_is_absolute), AV(exists), AV(is_dir)]
WritableDirAbs = Annotated[Path, AV(not_empty), AV(path_is_absolute), AV(dir_is_writable)]

# Field(min_length=1) error message confuses end-user; this just changes the error message
StrNotEmpty = Annotated[str, AV(not_empty)]
PathNotEmpty = Annotated[Path, AV(not_empty)]


class Timestep(BaseModel):
    num: float
    scale: Literal['ms', 'µs', 'us', 'ns', 'ps', 'fs']


class Project(BaseModel):
    id: int | str | None = Field(default=None)
    title: StrNotEmpty
    input_path: ExistingDirAbs
    output_path: WritableDirAbs
    python_path: ExistingFileAbs = Path(sys.executable)
    application_path: ExistingDirAbs = Path(stanalyzer.__path__[0])
    shell_path: ExistingFileAbs = Path('/bin/bash')
    traj: StrNotEmpty
    time_step: Annotated[Timestep, BV(str_to_args(Timestep))] = \
        Field(examples=['1 ns', '4.5 us', '.01 ms'])
    psf: PathNotEmpty
    scheduler: Literal['interactive', 'SLURM', 'PBS'] = 'interactive'
    SLURM: str = Field(default='')
    PBS: str = Field(default='')

    @field_validator('traj', mode='after')
    @classmethod
    def traj_is_valid_glob(cls, traj: str, info: ValidationInfo) -> str:
        if 'input_path' not in info.data:
            return ''  # input path raises its own error

        # check traj pattern matches at least one file
        traj_path = Path(traj)
        if traj_path.is_absolute():
            traj_relpathstr = '/'.join(traj_path.parts[1:])
            traj_glob = sorted(Path('/').glob(traj_relpathstr))
            as_relpath = Path(traj)
        else:
            as_relpath = traj_path
            traj_glob = sorted(info.data['input_path'].glob(traj))
        maybe_raise(traj_glob, 'glob_failed', traj)

        # enforce schema: traj is relative to input_path
        return str(as_relpath)

    @field_validator('psf', mode='after')
    @classmethod
    def psf_is_valid_path(cls, psf: Path, info: ValidationInfo) -> Path:
        if 'input_path' not in info.data:
            return Path()  # input path raises its own error

        psf_abspath: Path
        if psf.is_absolute():
            psf_abspath = psf
            as_relpath = psf.relative_to(info.data['input_path'])
        else:
            as_relpath = psf
            psf_abspath = info.data['input_path'] / psf
        maybe_raise(psf_abspath.is_file(), 'file_not_found', psf_abspath)

        # enforce schema: psf is relative to input_path
        return as_relpath


class User(BaseModel):
    name: str
    email: Optional[EmailStr]
