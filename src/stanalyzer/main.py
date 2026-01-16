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
from stanalyzer import utils, _release, __version__
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
ANALYSIS_CATEGORIES: dict[str, utils.AnalysisCategory] = \
    utils.read_yaml('static/analysis.category.yml')

if _release:
    MENU = utils.filter_unreleased(MENU)
    ANALYSIS = utils.filter_unreleased(ANALYSIS)
    ANALYSIS_CATEGORIES = utils.filter_empty_categories(ANALYSIS_CATEGORIES, ANALYSIS)

db.setup_db()


@app.get("/", response_class=HTMLResponse)
async def index(request: Request) -> HTMLResponse:
    return await page(request, "home")


@app.get("/{page}", response_class=HTMLResponse)
async def page(request: Request, page: str) -> HTMLResponse:
    context = {
        "request": request, "page_id": page, "menu": MENU,
        "page": PAGES.get(page, None), "auto_tooltip": utils.auto_tooltip,
        "call_macro_by_name": utils.call_macro_by_name, "__version__": __version__,
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

    ctx = invoke.Context(
        invoke.Config(overrides={'shell': project_settings['shell_path']})
    )
    processes = []
    for analysis, analysis_settings in job_settings.items():
        # setup invocation args
        program = f'stanalyzer {analysis}'

        # CLI opts expect '-' word sep, not '_'
        analysis_settings = {k.replace('_', '-'): v for k, v in analysis_settings.items()}

        bool_args = tuple(f'--{k} "{v}"' for k, v in analysis_settings.items()
                          if not isinstance(v, bool) and v != '')
        other_args = tuple(f'--{k}' for k, v in analysis_settings.items() if v is True)
        args = ' '.join(bool_args + other_args)

        args = f"{program} {args}"

        analysis_id = db.insert_analysis(db.Analysis(uid=1, project_id=project_id, args=args))

        # setup streams
        out_file = out_dir / analysis / f"analysis_{analysis_id}.out"
        err_file = out_dir / analysis / f"analysis_{analysis_id}.err"

        # ensure out_dir exists (default mode=0o777)
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True)

        if analysis_id is None:
            return {'error': 'failed to add analysis to database'}

        # auto-migration: settings from in_dir -> out_dir
        in_settings = in_dir / 'project.json'
        out_settings = out_dir / 'project.json'
        if not out_settings.exists():
            if in_settings.exists():
                in_settings.rename(out_settings)
            else:
                # user deleted project.json?
                response.status_code = status.HTTP_422_UNPROCESSABLE_ENTITY
                return [{'pid': -1, 'args': 'missing project.json'}]

        # TODO: this should ONLY be used when localhost is allowed
        with ctx.cd(out_dir):
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

    model_dict = settings.model_dump(mode="json")
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    db.update_project(db.Project(**model_dict))

    project_id = int(model_dict['id'])
    project_settings = db.get_user_projects(uid=1).get(project_id)
    assert project_settings

    out_dir = Path(project_settings['output_path'])
    settings_path = Path(out_dir / "project.json")
    utils.write_settings(settings_path, project_settings)

    return {
        'request_type': 'put',
        'data': db.get_user_projects(uid=1)
    }


@app.post("/project")
async def insert_project(request: Request,
                         settings: validation.Project) -> StrDict:
    time_step = settings.time_step

    model_dict = settings.model_dump(mode="json")
    model_dict.pop('id', None)  # user does not get to choose id
    model_dict['uid'] = 1   # TODO: get user id
    model_dict['time_step'] = f"{time_step.num} {time_step.scale}"

    project_id = db.insert_project(db.Project(**model_dict))

    project_settings = db.get_user_projects(uid=1).get(project_id)
    assert project_settings

    # create settings file
    out_dir = Path(project_settings['output_path'])
    settings_path = Path(out_dir / "project.json")
    if not settings_path.parent.exists():
        settings_path.parent.mkdir(parents=True)
    utils.write_settings(settings_path, project_settings)

    return {
        'request_type': 'post',
        'data': db.get_user_projects(uid=1),
    }


@app.delete("/project")
async def delete_project(settings: validation.Project) -> StrDict:
    model_dict = settings.model_dump(mode="json")
    model_dict['uid'] = 1   # TODO: get user id

    db.delete_project(pid=model_dict['id'], uid=1)

    return {
        'request_type': 'delete',
        'data': db.get_user_projects(uid=1),
    }
