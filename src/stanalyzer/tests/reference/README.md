# Golden Reference Files

This directory stores **golden reference output files** for stanalyzer's
correctness tests. Each reference file captures the expected output of an
analysis run against the standard test inputs (`src/stanalyzer/tests/inputs/`).

## Layout

Reference files are organized by analysis name, mirroring the output directory
structure that stanalyzer produces under `<output_path>/<analysis_name>/`:

```
reference/
  <analysis_name>/
    <output_file>.dat
    ...
```

For example, the reference for the `rmsf` analysis lives at
`reference/rmsf/rmsf.dat`.

## How references are generated

References are produced by running the analysis with the standard test
arguments against the standard test inputs, then copying the resulting output
files into this directory. The exact command used to generate each reference
is recorded in the test that consumes it (see `test_cli.py`).

## How references are consumed

The helper `assert_output_matches_reference()` in `test_cli.py` compares an
analysis's actual output against the corresponding reference file using
`numpy.testing.assert_allclose` with a default relative tolerance of `1e-5`
and absolute tolerance of `1e-8`.

## Regenerating references

When an analysis's output format legitimately changes (e.g. a bug fix that
alters numerical output, or a new feature that adds columns), the reference
files must be regenerated. To regenerate:

1. Run the analysis with the standard test arguments against the standard
   test inputs.
2. Copy the resulting output file(s) into `reference/<analysis_name>/`.
3. Verify the new reference is correct (e.g. by inspecting the output or
   cross-checking against an independent implementation).
4. Commit the updated reference files alongside the code change.

> **Warning**: Do not regenerate references to silence a failing test. A
> failing comparison usually indicates a real regression. Only regenerate
> when the output change is intentional and verified.

### 2026-09 regeneration: contacts, pi_stacking, rmsf

These three references were regenerated on the merge of `test-updates` into
`master`, after master's analysis rewrite (PR #19) intentionally changed their
output. `contacts` and `pi_stacking` emit the same data in a new, deterministic
sorted row order; `rmsf` differs only by float32-vs-float64 alignment
accumulation (maximum relative difference ~2.9e-05). All three were verified
algorithm-equivalent before regeneration (identical sorted datasets; identical
rmsf index column, all values within 1e-4 relative), so this regeneration is
intentional and verified — not an attempt to silence a failing test. To
regenerate them again: `python generate_baselines.py --only contacts pi_stacking rmsf`.
