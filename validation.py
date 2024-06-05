import os
from pathlib import Path
from typing import Annotated, Optional
from typing_extensions import Self

from pydantic import BaseModel, EmailStr, model_validator
from pydantic.functional_validators import AfterValidator as AV


def exists(p: Path) -> Path:
    assert p.exists(), f"No such file or directory: {p}"
    return p


def is_file(p: Path) -> Path:
    assert p.is_file(), f"{p} is not a regular file"
    return p


def is_dir(p: Path) -> Path:
    assert p.is_dir(), f"{p} is not a directory"
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


ExistingFile = Annotated[Path, AV(exists), AV(is_file)]
ExistingDir = Annotated[Path, AV(exists), AV(is_dir)]
WritableDir = Annotated[Path, AV(dir_is_writable)]


class Timestep(BaseModel):
    num: float
    scale: str


class Project(BaseModel):
    id: Optional[int | str]
    title: str
    input_path: ExistingDir
    output_path: WritableDir
    python_path: ExistingFile
    application_path: ExistingDir
    shell_path: ExistingFile
    traj: str
    time_step: Timestep
    psf: Path
    scheduler: str
    SLURM: Optional[str]
    PBS: Optional[str]

    @model_validator(mode='after')
    def valid_relpaths(self) -> Self:
        # check PSF relpath is a valid file
        psf_abspath = self.input_path / self.psf
        assert psf_abspath.is_file(), "No such file or directory: {psf_abspath}"

        # check traj pattern matches at least one file
        traj_glob = sorted(self.input_path.glob(self.traj))
        assert traj_glob, "No files matched by pattern: {self.traj}"

        # convert Path -> str b/c SQLModel apparently doesn't do it for us
        for field, value in self:
            if isinstance(value, Path):
                setattr(self, field, str(value))

        return self


class User(BaseModel):
    name: str
    email: Optional[EmailStr]
