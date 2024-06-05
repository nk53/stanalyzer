import json
import pathlib
import re
from typing import Any, List, Optional

import jinja2
import yaml

# groups strings[like][this] into their individual tokens
P_RE = re.compile(r'^([^[]*)|\[([^]]*)\]')
TOOLTIPS = None


def read_yaml(filename: str) -> Any:
    """Load a YAML file by name"""
    with open(filename) as file_obj:
        return yaml.load(file_obj, Loader=yaml.FullLoader)


def read_json(filename: str) -> Any:
    """Load a JSON file by name"""
    with open(filename) as file_obj:
        return json.load(file_obj)


def write_settings(path: pathlib.Path, data: Any, mode: str = 'w'):
    if not path.parent.exists():
        path.parent.mkdir()
    with path.open(mode=mode) as file_obj:
        json.dump(data, file_obj)


@jinja2.pass_context
def call_macro_by_name(ctx, name, *args, **kwargs):
    for varname in name.split('.'):
        getter = getattr(ctx, 'get', None)
        if getter is not None:
            ctx = ctx.get(varname)
        elif isinstance(ctx, jinja2.environment.TemplateModule):
            ctx = getattr(ctx, varname)
    return ctx(*args, **kwargs)


def auto_tooltip(name: str, debug: bool = False) -> Optional[str]:
    """Return tooltip string if it exists, else None

    Tooltip names derive from the field name, which is subscripted with
    either ``underscores_like_this`` or ``brackets[like][this]``. Brackets are
    checked first, then subscripts.

    If both of those checks fail, then one last check is performed by
    inserting ``_`` between keys from right to left, e.g.:

    .. code-block:: text

        brackets[like][this] --> KeyError
        brackets[like_this] --> KeyError
        brackets_like_this --> match

    Thus, only leaf nodes may contain underscores in tooltips.yml.

    Finally, if no match is found, then :python:`None` is returned.
    """
    global TOOLTIPS

    if TOOLTIPS is None:
        TOOLTIPS = read_yaml('static/tooltips.yml')

    def lookup_keys(keys):
        if debug:
            print(f'lookup_keys({keys}) = ', end='')

        # descend dict tree until reaching leaf, null, or keys are exhausted
        tooltip = TOOLTIPS
        while isinstance(tooltip, dict) and keys:
            key, keys = keys[0], keys[1:]
            tooltip = tooltip.get(key, None)

        if debug:
            print(f'{tooltip}')

        if not isinstance(tooltip, str):
            return None

        return tooltip

    def get_tooltip(keys):
        tooltip = lookup_keys(keys)

        if tooltip is not None:
            return tooltip

        # reconstruct with underscores from right to left
        key = ''
        while keys and tooltip is None:
            if key:
                key = f"{keys.pop()}_{key}"
            elif len(keys) > 1:
                key, keys = '_'.join(keys[-2:]), keys[:-2]
            else:
                break

            tooltip = lookup_keys(keys + [key])

        if debug:
            print(f'get_tooltip() = {tooltip}')

        return tooltip

    if '[' in name:
        keys = [''.join(n) for n in re.findall(P_RE, name)]

        return get_tooltip(keys)

    if '_' in name:
        # transform name_like_this to [name][like][this]
        keys = name.split('_')
        return get_tooltip(keys)

    return get_tooltip([name])


def get_active_settings(settings: dict, analysis: dict, path: Optional[List[str]] = None) -> dict:
    if path is None:
        path = []

    def get_setting(settings, option):
        return settings.get('_'.join(path + [option]), False)

    results = {}
    for analysis_type in analysis.keys():
        analysis_key = '_'.join(path + [analysis_type])
        if settings.get(analysis_key, None):
            result = {}
            if sub_opts := {k: v for (k,v) in analysis[analysis_type]['options'].items()
                            if isinstance(v, dict) and 'options' in v}:
                result.update(
                    get_active_settings(settings, sub_opts, path + [analysis_type])
                )
            for opt, attrs in analysis[analysis_type]['options'].items():
                if 'options' not in attrs:
                    result[opt] = get_setting(settings, analysis_type+'_'+opt)

            results[analysis_type] = result

    return results
