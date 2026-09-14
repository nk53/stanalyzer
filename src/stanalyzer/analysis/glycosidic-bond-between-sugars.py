import argparse
from typing import Optional

import numpy as np

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda

ANALYSIS_NAME = 'glycosidic_bond_between_sugars'


def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time contacts"
    else:
        header_str = "#time contacts"

    print(header_str, file=outfile)

    return header_str


def _bond_label(o, c) -> str:
    """Label for a glycosidic bond, e.g. 'O3(CARA1BGALNA)-C1(CARA2BGAL)'.

    Lower residue number comes first for a stable column order."""
    first, second = (o, c) if (o.resid, o.name) < (c.resid, c.name) else (c, o)
    return (f'{first.name}({first.segid}{first.resid}{first.resname})-'
            f'{second.name}({second.segid}{second.resid}{second.resname})')


def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   out: sta.FileRef, interval: int = 1) -> None:
    """Writes per-frame glycosidic bond lengths to `out`.

    A glycosidic bond is any covalent bond between a carbon and an oxygen
    in *different* residues, with both endpoints inside `sel`. Endpoint
    identities come from the topology's bond list, not atom names, so the
    anomeric carbon (C1 in aldoses, C2 in ketoses such as sialic acid) and
    the acceptor oxygen (O1-O6, O11/12, OG1, ...) are found automatically.
    """
    u = mda.Universe(psf, traj)

    sel_atoms = u.select_atoms(sel or 'all')
    sel_ids = set(sel_atoms.ix)
    oxygens = sel_atoms.select_atoms('name O* and bonded name C*')

    bonds: dict[frozenset, tuple] = {}
    for o in oxygens:
        for c in o.bonded_atoms:
            if c.ix not in sel_ids or not c.name.startswith('C'):
                continue
            if c.resid == o.resid:
                continue
            bonds[frozenset({o.ix, c.ix})] = (o, c)

    bond_pairs = sorted(bonds.values(),
                        key=lambda pair: (pair[0].resid, pair[0].name,
                                          pair[1].resid, pair[1].name))

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        if bond_pairs:
            labels = '  '.join(_bond_label(o, c) for o, c in bond_pairs)
            print(f'# {labels}', file=outfile)

        for ts in u.trajectory[::interval]:
            values = ' '.join(
                f'{np.linalg.norm(o.position - c.position):.6f}'
                for o, c in bond_pairs)
            print(f'{ts.time:10.3f} {values}', file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for glycosidic_bond calculation")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()
