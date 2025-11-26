import os
import psutil
from pathlib import Path
from typing import Literal, Optional, overload
from sqlmodel import Field, Session, SQLModel, select


class User(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    email: Optional[str] = None


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
    args: str = ''

    __basepath: str = ''
    __table_args__ = {"sqlite_autoincrement": True}

    @overload
    def read(self, which: Literal['both'] = 'both') -> tuple[str, str]: ...

    @overload
    def read(self, which: Literal['out', 'err']) -> str: ...

    def read(self, which: Literal['out', 'err', 'both'] = 'both') \
            -> str | tuple[str, str]:
        def _read(p: Path) -> str:
            # don't read all of a large file (>2 pages)
            if not p.exists():
                return ''
            if os.stat(p).st_size > 0x2000:
                # first page, '[...]', & last page
                result: list[str] = []
                with p.open() as file:
                    result.extend([file.read(0x1000), '[...]'])
                    file.seek(-0x1000, os.SEEK_END)
                    result.append(file.read())
                return '\n'.join(result)
            # return all of a small file (<2 pages)
            return p.read_text()

        analysis_name = self.args.split()[1]
        project_dir = Path(self.get_basepath()) / analysis_name
        outfile = project_dir / f"analysis_{self.id}.out"
        errfile = outfile.with_suffix(".err")

        result_paths = outfile, errfile
        to_read = which in ('out', 'both'), which in ('err', 'both')

        results = tuple(_read(rp) for rp, tr in zip(result_paths, to_read) if tr)

        if len(results) == 1:
            return results[0]

        assert len(results) == 2
        return results

    def set_basepath(self, session: Session | None = None, path: str = '') -> None:
        """Either session or path must be specified. If both are given, path is ignored"""
        if session:
            statement = select(Project).where(Project.uid == self.uid)
            project = session.exec(statement).one()
            self.__basepath = project.output_path
        else:
            self.__basepath = path

    def get_basepath(self, session: Session | None = None) -> str:
        if session is not None:
            self.set_basepath(session)
        return self.__basepath

    def done(self, session: Session | None = None) -> bool:
        """Returns False if the associated process is still running.

        If session is not None, also updates DB if process has ended.
        """
        # There is a slight chance of false positives if a lot of time has
        # passed and another submitted analysis has the same arguments.
        # This *should* be exceedingly rare.
        if not self.process_id:
            return True

        username = psutil.Process().username()
        match = any(proc for proc in psutil.process_iter(['pid', 'username'])
                    if (info := proc.info)['pid'] == self.process_id and
                    info['username'] == username)

        if session is not None and not match:
            self.process_id = 0
            session.add(self)

        return not match
