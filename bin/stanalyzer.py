import argparse
import os
from importlib import import_module
from typing import Optional

from stanalyzer.utils import read_json


def get_settings(analysis_name: Optional[str] = None) -> dict:
    parser = argparse.ArgumentParser(description='TODO')

    if analysis_name is None:
        parser.add_argument('analysis_name')

    parser.add_argument('setting=value', nargs='*',
                        help='can be used to overwrite settings.json OR to specify'
                        ' them on the command-line rather than GUI')
    args = parser.parse_args()

    if os.path.exists('settings.json'):
        settings = read_json('settings.json')

    settings.update(args)

    return settings


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = get_settings()

    analysis_name = settings.pop('analysis_name')
    analysis = import_module(f"stanalyzer.{analysis_name}")

    analysis.main(settings)


if __name__ == '__main__':
    main()
