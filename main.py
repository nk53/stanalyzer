# built-ins
import json
import os
import signal
import sys
from pathlib import Path

# webserver modules
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates

# CLI modules
import invoke

# intra-package
from stanalyzer import utils
from stanalyzer.db import db
from stanalyzer import validation

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

templates = Jinja2Templates(directory="templates",
                            extensions=['jinja2.ext.do', 'jinja2.ext.debug'])
MENU = utils.read_yaml('static/menu.yml')
PAGES = utils.read_yaml('static/pages.yml')
ANALYSIS = utils.read_yaml('static/analysis.yml')


db.setup_db()


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
    project_id = int(settings['project'])
    project_settings = db.get_user_projects(uid=1).get(project_id)

    in_dir = Path(project_settings['input_path'])
    out_dir = Path(project_settings['output_path'])
    settings_path = Path(in_dir / "project.json")
    utils.write_settings(settings_path, project_settings)

    ctx = invoke.Context()
    processes = []
    for analysis, analysis_settings in job_settings.items():
        # setup streams
        analysis_id = db.insert_analysis(db.Analysis(uid=1, project_id=project_id))
        out_file = Path(out_dir / f"analysis_{analysis_id}.out")
        err_file = Path(out_dir / f"analysis_{analysis_id}.err")

        # this should ONLY be used when localhost is allowed
        with ctx.cd(in_dir):
            # setup invocation args
            program = f'python -m stanalyzer {analysis}'
            args = ' '.join(f'--{k} {v}' for k,v in analysis_settings.items())
            args = f"{program} {args}"

            # log invocation for user
            out_stream = out_file.open('w')
            print('args:', args, file=out_stream)
            out_stream.flush()  # avoids race cond. w/ ctx.run

            # ignore end of child process; prevents "zombie" proc status
            signal.signal(signal.SIGCHLD, signal.SIG_IGN)

            # init program
            promise = ctx.run(args, asynchronous=True,
                              out_stream=out_stream, err_stream=err_file.open('w'))

            # get PID and update analysis info
            if isinstance(promise, invoke.Promise) and isinstance(promise.runner, invoke.Local):
                pid = promise.runner.process.pid
                db.set_analysis_pid(analysis_id=analysis_id, process_id=pid)
                processes.append({'args': args, 'pid': pid})
            else:
                print('something has gone wrong with ', analysis)

    print(processes)
    return processes


@app.put("/project")
async def update_project(settings: validation.Project):
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
async def insert_project(request: Request, settings: validation.Project):
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
async def delete_project(settings: validation.Project):
    model_dict = settings.model_dump()
    model_dict['uid'] = 1   # TODO: get user id

    db.delete_project(pid=model_dict['id'], uid=1)

    return {
        'request_type': 'delete',
        'data': db.get_user_projects(uid=1),
    }
