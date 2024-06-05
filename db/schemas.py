from typing import Optional
from sqlmodel import Field, SQLModel


class User(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    email: Optional[str]


class Project(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    uid: int = Field(foreign_key="user.id")
    title: str
    input_path: str
    output_path: str
    python_path: str
    application_path: str
    shell_path: str
    psf: str
    traj: str
    time_step: str
    scheduler: str
    SLURM: Optional[str]
    PBS: Optional[str]


class Analysis(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    uid: int = Field(foreign_key="user.id")
    project_id: int = Field(foreign_key="project.id")
    process_id: Optional[int] = None
