# built-ins
import json
import os
import signal
import sys
from pathlib import Path

# webserver modules
from fastapi import FastAPI, Request, Response, status
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates

# CLI modules
import invoke

# intra-package
from stanalyzer import utils, _release
from stanalyzer.db import db
from stanalyzer import validation
from stanalyzer._typing import Any, StrDict, StrDictList

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

templates = Jinja2Templates(directory="templates",
                            extensions=['jinja2.ext.do', 'jinja2.ext.debug'])
MENU = utils.read_yaml('static/menu.yml')
PAGES = utils.read_yaml('static/pages.yml')
ANALYSIS = utils.read_yaml('static/analysis.yml')
ANALYSIS_CATEGORIES = utils.read_yaml('static/analysis.category.yml')

if _release:
    MENU = utils.filter_unreleased(MENU)
    ANALYSIS = utils.filter_unreleased(ANALYSIS)

db.setup_db()


@app.get("/", response_class=HTMLResponse)
async def index(request: Request) -> HTMLResponse:
    return await page(request, "home")


@app.get("/{page}", response_class=HTMLResponse)
async def page(request: Request, page: str) -> HTMLResponse:
    context = {
        "request": request, "page_id": page, "menu": MENU,
        "page": PAGES.get(page, None), "auto_tooltip": utils.auto_tooltip,
        "call_macro_by_name": utils.call_macro_by_name,
    }
    kwargs: dict[str, Any] = {}

    if page not in PAGES:
        context['page'] = PAGES["404"]
        context['page_id'] = "404"
        kwargs['status_code'] = 404
    elif page == 'analysis':
        context['analysis_form'] = ANALYSIS
        context['projects_menu'] = [
            (p['id'], p['title']) for p in db.get_user_projects(uid=1).values()
        ]
        context['categories'] = ANALYSIS_CATEGORIES
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
    elif page == 'results':
        projects = db.get_user_projects(uid=1)
        results = {project_id: db.get_analysis_results(project_id)
                   for project_id in projects}
        context['analysis'] = ANALYSIS
        context['projects'] = projects
        context['results'] = results

    return templates.TemplateResponse("layouts/single_page.html", context, **kwargs)


@app.delete("/analysis/{analysis_id}")
async def delete_analysis(analysis_id: int) -> int:
    # TODO: auth
    return db.delete_analysis(analysis_id)


@app.post("/analysis")
async def insert_analysis(settings: dict, response: Response) -> StrDict | StrDictList:
    job_settings = utils.get_active_settings(settings, ANALYSIS)

    if not job_settings:
        response.status_code = status.HTTP_422_UNPROCESSABLE_ENTITY
        return []

    project_id = int(settings['project'])
    project_settings = db.get_user_projects(uid=1).get(project_id)

    assert project_settings

    in_dir = Path(project_settings['input_path'])
    out_dir = Path(project_settings['output_path'])

    ctx = invoke.Context()
    processes = []
    for analysis, analysis_settings in job_settings.items():
        # setup invocation args
        program = f'stanalyzer {analysis}'
        bool_args = tuple(f'--{k} "{v}"' for k, v in analysis_settings.items()
                          if not isinstance(v, bool))
        other_args = tuple(f'--{k}' for k, v in analysis_settings.items() if v is True)
        args = ' '.join(bool_args + other_args)

        args = f"{program} {args}"

        analysis_id = db.insert_analysis(db.Analysis(uid=1, project_id=project_id, args=args))

        # setup streams
        out_file = Path(out_dir / f"analysis_{analysis_id}.out")
        err_file = Path(out_dir / f"analysis_{analysis_id}.err")

        if analysis_id is None:
            return {'error': 'failed to add analysis to database'}

        # TODO: this should ONLY be used when localhost is allowed
        with ctx.cd(in_dir):
            # CLI opts expect '-' word sep, not '_'
            analysis_settings = {k.replace('_', '-'): v
                                 for k, v in analysis_settings.items()}

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
                processes.append({'args': args, 'pid': -1})
                print('Something unexpected went wrong with', analysis)

    print(processes)
    return processes


@app.put("/project")
async def update_project(settings: validation.Project) -> StrDict:
    time_step = settings.time_step

    model_dict = settings.model_dump()
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    db.update_project(db.Project(**model_dict))

    project_id = int(model_dict['id'])
    project_settings = db.get_user_projects(uid=1).get(project_id)
    assert project_settings

    in_dir = Path(project_settings['input_path'])
    settings_path = Path(in_dir / "project.json")
    utils.write_settings(settings_path, project_settings)

    return {
        'request_type': 'put',
        'data': db.get_user_projects(uid=1)
    }


@app.post("/project")
async def insert_project(request: Request,
                         settings: validation.Project) -> StrDict:
    time_step = settings.time_step

    model_dict = settings.model_dump()
    model_dict.pop('id', None)  # user does not get to choose id
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    project_id = db.insert_project(db.Project(**model_dict))

    project_settings = db.get_user_projects(uid=1).get(project_id)
    assert project_settings

    # create settings file
    in_dir = Path(project_settings['input_path'])
    settings_path = Path(in_dir / "project.json")
    utils.write_settings(settings_path, project_settings)

    return {
        'request_type': 'post',
        'data': db.get_user_projects(uid=1),
    }


@app.delete("/project")
async def delete_project(settings: validation.Project) -> StrDict:
    model_dict = settings.model_dump()
    model_dict['uid'] = 1   # TODO: get user id

    db.delete_project(pid=model_dict['id'], uid=1)

    return {
        'request_type': 'delete',
        'data': db.get_user_projects(uid=1),
    }
