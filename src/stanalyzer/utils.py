import json
import pathlib
import re
from glob import glob
from pathlib import Path
from typing import Any, List, TypeVar

import bracex
import yaml

T = TypeVar('T')

try:
    import jinja2

    @jinja2.pass_context
    def call_macro_by_name(ctx, name, *args, **kwargs):
        for varname in name.split('.'):
            getter = getattr(ctx, 'get', None)
            if getter is not None:
                ctx = ctx.get(varname)
            elif isinstance(ctx, jinja2.environment.TemplateModule):
                ctx = getattr(ctx, varname)
        return ctx(*args, **kwargs)
except ModuleNotFoundError:
    def call_macro_by_name(*args, **kwargs):
        raise NotImplementedError(
            "'client' build does not support jinja2 templating")


# groups strings[like][this] into their individual tokens
P_RE = re.compile(r'^([^[]*)|\[([^]]*)\]')
TOOLTIPS: dict[str, Any] | None = None


def load_tooltips() -> dict[str, Any]:
    tooltips = read_yaml('static/tooltips.yml')
    assert isinstance(tooltips, dict)
    return tooltips


def read_yaml(filename: str) -> Any:
    """Load a YAML file by name"""
    with open(filename) as file_obj:
        return yaml.load(file_obj, Loader=yaml.FullLoader)


def read_json(filename: str) -> Any:
    """Load a JSON file by name"""
    with open(filename) as file_obj:
        return json.load(file_obj)


def write_settings(path: pathlib.Path, data: Any, mode: str = 'w') -> None:
    if not path.parent.exists():
        path.parent.mkdir()
    with path.open(mode=mode) as file_obj:
        json.dump(data, file_obj, indent=4)


def auto_tooltip(name: str, debug: bool = False) -> str | None:
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
        TOOLTIPS = load_tooltips()

    def lookup_keys(keys: list[str]) -> str | None:
        if debug:
            print(f'lookup_keys({keys}) = ', end='')

        # descend dict tree until reaching leaf, null, or keys are exhausted
        tooltip: dict[str, Any] | str | None = TOOLTIPS
        while isinstance(tooltip, dict) and keys:
            key, keys = keys[0], keys[1:]
            tooltip = tooltip.get(key, None)

        if debug:
            print(f'{tooltip}')

        if not isinstance(tooltip, str):
            return None

        return tooltip

    def get_tooltip(keys: list[str]) -> str | None:
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


def get_active_settings(settings: dict[str, Any], analysis: dict[str, Any],
                        path: List[str] | None = None) -> dict:
    if path is None:
        path = []

    def get_setting(settings: dict[str, Any], option: str) -> Any:
        nonlocal path
        return settings.get('_'.join(path + [option]), False)

    def has_sub_opts(values: Any) -> bool:
        if not isinstance(values, dict):
            return False
        if 'options' not in values:
            return False
        if 'select' in values.get('type', ''):
            return False
        return True

    def allowed(setting_name: str, required_value: object) -> bool:
        nonlocal settings
        value = settings.get(setting_name)
        if isinstance(required_value, bool):
            value = bool(value)
        return value == required_value

    results: dict = {}
    for analysis_type in analysis.keys():
        analysis_key = '_'.join(path + [analysis_type])
        if settings.get(analysis_key, None):
            if path:
                result = results
            else:
                result = {}

            if sub_opts := {k: v for (k, v) in analysis[analysis_type]['options'].items()
                            if has_sub_opts(v)}:
                result.update(
                    get_active_settings(settings, sub_opts, path + [analysis_type])
                )
            for opt, attrs in analysis[analysis_type]['options'].items():
                opt_long = f"{analysis_type}_{opt}"
                if path:
                    opt = opt_long
                opt_type = attrs.get('type', '')
                if not has_sub_opts(attrs):
                    setting = get_setting(settings, opt_long)
                    if 'checkbox' in opt_type:
                        setting = bool(setting)
                elif 'select' in opt_type:
                    setting = get_setting(settings, opt_long)
                else:
                    continue

                _requires = attrs.get('_requires', {})
                if all(allowed(name, req) for name, req in _requires.items()):
                    result[opt] = setting

            if not path:
                results[analysis_type] = result

    return results


def filter_unreleased(d: dict[str, dict[str, T]]) -> dict[str, dict[str, T]]:
    """Hide unfinished items for release and release-candidate builds"""
    return {k: v for k, v in d.items() if v.get('release', True)}


def braced_glob(p: Path | str) -> list[Path]:
    result = [
        Path(globbed)
        for expanded in bracex.iexpand(str(p), limit=0)
        for globbed in glob(expanded)
    ]

    return sorted(result)
