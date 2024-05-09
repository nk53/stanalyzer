# built-ins
import json
import os
import sys
from typing import Optional

# webserver modules
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates
from pydantic import BaseModel, EmailStr

# CLI modules
import invoke

# intra-package
from . import utils
from .bin import analyze
from .db import db

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

templates = Jinja2Templates(directory="templates",
                            extensions=['jinja2.ext.do', 'jinja2.ext.debug'])
MENU = utils.read_yaml('static/menu.yml')
PAGES = utils.read_yaml('static/pages.yml')
ANALYSIS = utils.read_yaml('static/analysis.yml')


db.setup_db()


class Timestep(BaseModel):
    num: float
    scale: str


class Project(BaseModel):
    id: Optional[int | str]
    title: str
    input_path: str
    output_path: str
    python_path: str
    application_path: str
    shell_path: str
    traj_pattern: str
    time_step: Timestep
    psf: str
    scheduler: str
    SLURM: Optional[str]
    PBS: Optional[str]


class User(BaseModel):
    name: str
    email: Optional[EmailStr]


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return await page(request, "home")


@app.get("/{page}", response_class=HTMLResponse)
async def page(request: Request, page: str):
    if page not in PAGES:
        page = "404"

    context = {
        "request": request, "page_id": page, "menu": MENU,
        "page": PAGES[page], "auto_tooltip": utils.auto_tooltip,
        "call_macro_by_name": utils.call_macro_by_name,
    }

    if page == 'analysis':
        context['analysis_form'] = ANALYSIS
        context['projects_menu'] = [
            (p['id'], p['title']) for p in db.get_user_projects(uid=1).values()
        ]

    elif page == 'project':
        projects = db.get_user_projects(uid=1)
        projects_json = projects.copy()
        projects_json['add_new'] = {'PBS': '#!/bin/bash\nhi'}
        context['projects'] = json.dumps(projects_json)
        context['projects_menu'] = [
            *[(p['id'], p['title']) for p in projects.values()],
            ('add_new', 'Add New Project'),
        ]
        context['application_path'] = os.getcwd()
        context['python_path'] = sys.executable

    return templates.TemplateResponse("layouts/single_page.html", context)


@app.post("/analysis")
async def insert_analysis(settings: dict):
    job_settings = utils.get_active_settings(settings, ANALYSIS)

    # this should ONLY be used when localhost is allow
    invoke.run(f'python bin/stanalyzer.py')

    return job_settings


@app.put("/project")
async def update_project(settings: Project):
    time_step = settings.time_step

    model_dict = settings.model_dump()
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    db.update_project(db.Project(**model_dict))

    return {
        'request_type': 'put',
        'data': db.get_user_projects(uid=1)
    }


@app.post("/project")
async def insert_project(request: Request, settings: Project):
    time_step = settings.time_step

    model_dict = settings.model_dump()
    model_dict.pop('id', None)  # user does not get to choose id
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    db.insert_project(db.Project(**model_dict))

    return {
        'request_type': 'post',
        'data': db.get_user_projects(uid=1),
    }


@app.delete("/project")
async def delete_project(settings: Project):
    model_dict = settings.model_dump()
    model_dict['uid'] = 1   # TODO: get user id

    db.delete_project(pid=model_dict['id'], uid=1)

    return {
        'request_type': 'delete',
        'data': db.get_user_projects(uid=1),
    }
