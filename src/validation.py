import os
import sys
from pathlib import Path
from typing import (
    Any,
    Annotated,
    Literal,
    Optional,
    TypeAlias,
    TypeVar,
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
import stanalyzer

T = TypeVar('T', bound=type)
BM = TypeVar('BM', bound=BaseModel)


def exists(p: Path) -> Path:
    assert p.exists(), f"No such file or directory: {p}"
    return p


def is_file(p: Path) -> Path:
    assert p.is_file(), f"{p} is not a regular file"
    return p


def is_dir(p: Path) -> Path:
    assert p.is_dir(), f"{p} is not a directory"
    return p


def path_is_absolute(p: Path) -> Path:
    assert p.is_absolute(), f"Not an absolute path: {p}"
    return p


def dir_is_writable(p: Path) -> Path:
    """Check that we can write to OR create a dir"""
    if p.is_dir():
        assert os.access(p, os.W_OK), f"No permission to write to {p}"
    else:
        # this approach doesn't require abspath args
        parent = (p / '..').resolve()
        assert os.access(parent, os.W_OK), f"No permission to create {p}"
    return p


def str_to_args(cls: type[BM]):
    fields = cls.model_fields
    keys = tuple(fields.keys())

    def _str_to_args(value: Any, info: ValidationInfo) -> Any:
        if isinstance(value, str) and len(keys) > 1:
            model_dict = dict(zip(keys, value.split()))
            return cls(**model_dict)

        return value

    return _str_to_args


ExistingFile = Annotated[Path, AV(exists), AV(is_file)]
ExistingDir = Annotated[Path, AV(exists), AV(is_dir)]
WritableDir = Annotated[Path, AV(dir_is_writable)]
ExistingFileAbs = Annotated[Path, AV(path_is_absolute), AV(exists), AV(is_file)]
ExistingDirAbs = Annotated[Path, AV(path_is_absolute), AV(exists), AV(is_dir)]
WritableDirAbs = Annotated[Path, AV(path_is_absolute), AV(dir_is_writable)]


class Depends:
    """Annotation metadata indicating that a field depends on another field"""
    fields: tuple[str, ...]

    def __init__(self, *fields: str):
        self.fields = fields

    def __class_getitem__(cls, *key: str):
        """Not a real type parameterization; just convenience"""
        return cls(*key)


class Timestep(BaseModel):
    num: float
    scale: Literal['ms', 'µs', 'us', 'ns', 'ps', 'fs']


class Project(BaseModel):
    id: int | str | None = Field(default=None)
    title: str = Field(min_length=1)
    input_path: ExistingDirAbs
    output_path: WritableDirAbs
    python_path: ExistingFileAbs = Path(sys.executable)
    application_path: ExistingDirAbs = Path(stanalyzer.__path__[0])
    shell_path: ExistingFileAbs = Path('/bin/bash')
    traj: Annotated[str, Depends['input_path']] = Field(min_length=1)
    time_step: Annotated[Timestep, BV(str_to_args(Timestep))] = \
        Field(examples=['1 ns', '4.5 us', '.01 ms'])
    psf: Annotated[Path, Depends['input_path']]
    scheduler: Literal['interactive', 'SLURM', 'PBS'] = 'interactive'
    SLURM: str = Field(default='')
    PBS: str = Field(default='')

    @model_validator(mode='after')
    def traj_is_valid_glob(self) -> Self:
        # check traj pattern matches at least one file
        traj = Path(self.traj)
        if traj.is_absolute():
            traj_relpathstr = '/'.join(traj.parts[1:])
            traj_glob = sorted(Path('/').glob(traj_relpathstr))
            as_relpath = Path(traj)
        else:
            as_relpath = traj
            traj_glob = sorted(self.input_path.glob(self.traj))
        assert traj_glob, "No files matched by pattern: {self.traj}"

        # enforce schema: traj is relative to input_path
        self.traj = str(as_relpath)
        return self

    @model_validator(mode='after')
    def psf_is_valid_path(self) -> Self:
        psf = Path(self.psf)
        if psf.is_absolute():
            psf_abspath = psf
            as_relpath = psf.relative_to(self.input_path)
        else:
            as_relpath = psf
            psf_abspath = self.input_path / psf
        assert psf_abspath.is_file(), "No such file or directory: {psf_abspath}"

        # enforce schema: psf is relative to input_path
        self.psf = as_relpath
        return self

    @model_validator(mode='after')
    def valid_paths_to_str(self) -> Self:
        """Convert Path -> str for SQLModel"""
        for field, value in self:
            if isinstance(value, Path):
                setattr(self, field, str(value))

        return self


class User(BaseModel):
    name: str
    email: Optional[EmailStr]
