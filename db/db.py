import os
from sqlmodel import Session, SQLModel, create_engine, select
from .schemas import (
    User,
    Project,
)

DB_PATH = 'database.db'
engine = create_engine(f"sqlite:///{DB_PATH}")


def setup_db():
    if not os.path.exists(DB_PATH):
        SQLModel.metadata.create_all(engine)

        with Session(engine) as session:
            session.add(User(name="admin"))
            session.commit()


def get_user_projects(uid: int):
    with Session(engine) as session:
        statement = select(Project).where(Project.uid == uid)
        projects = [p.model_dump() for p in session.exec(statement).all()]

        project_dict = {p['id']: p for p in projects}
        for pid, model in project_dict.items():
            # user does not need to know their id
            del project_dict[pid]['uid']

        return project_dict


def update_project(settings: Project):
    with Session(engine) as session:
        # get matching project
        statement = select(Project).where(Project.uid == settings.uid, Project.id == settings.id)
        project = session.exec(statement).one()

        # update new values
        for field in project.model_fields.keys():
            value = getattr(settings, field)
            setattr(project, field, value)

        # commit change
        session.add(project)
        session.commit()


def insert_project(settings: Project):
    with Session(engine) as session:
        session.add(settings)
        session.commit()


def delete_project(pid: int, uid: int):
    with Session(engine) as session:
        statement = select(Project).where(Project.uid == uid, Project.id == pid)
        project = session.exec(statement).one()
        session.delete(project)
        session.commit()
