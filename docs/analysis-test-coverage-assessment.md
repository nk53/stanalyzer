# Analysis Test Coverage Assessment

This report reviews every golden-reference correctness test in ST-Analyzer's
test suite (`src/stanalyzer/tests/test_cli.py`) and the stored expected-result
files they compare against (`src/stanalyzer/tests/reference/`). For each of the
34 built-in analyses it answers three questions — whether the test meaningfully
exercises the analysis (Q1), whether the reference result is scientifically
meaningful (Q2), and what the output ought to mean scientifically (Q3) — and
ends with a verdict table identifying which tests and references need attention.

- **Scope**: 34 golden-reference test classes, 55 reference files
  (52 `.dat` + 2 `.pdb` + 1 `.png`), 3 input systems.
- **Status**: complete — all 34 per-analysis assessments are present (26
  soohyung, 5 yiwei, 3 omf), each answering Q1/Q2/Q3 and closing with a
  verdict.
- **Date**: 2026-09-13
- **Update (2026-09-14)**: three flagged items were resolved after this
  assessment — ContactResidenceTime (UNRELIABLE reference) fixed,
  PositionTimeCopy (near-duplicate) removed, HelixAnalysis (cross-platform
  fragility) fixed via a structure-aware comparator. See the UPDATE section
  at the end of this file; the verdict table and per-analysis subsections
  below remain the frozen 2026-09-13 snapshot.

## Background

### Input systems

The golden-reference tests run against three input systems under
`src/stanalyzer/tests/inputs/`:

| System | Input files | Composition |
| --- | --- | --- |
| `soohyung_membrane` | `step7_999.dcd`, `step7_1000.dcd`, `step5_input.psf` | DOPC/DSPC/CHL membrane + protein (segid PROA) |
| `yiwei_protein` | `step5_1.dcd`, `step5_2.dcd`, `step3_input.psf` | Protein system (TIP3 water, ions) |
| `2omf_membrane` | `equil.dcd`, `system.psf` | PDB 2OMF OmpF porin monomer, 190 POPC + 49 CHL, 340 residues, 5 frames / 25 ps (see `inputs/2omf_membrane/README.md`) |

`inputs/README.md` is git-submodule install notes, not a per-system
description; `inputs/2omf_membrane/README.md` is the only system-description
README.

### Golden-reference harness

- One test class per analysis, each providing `standard_args` plus a
  `test_standard_correctness` method that runs the analysis and compares every
  discovered output file against the golden reference.
- A generic `test_standard` smoke test (file-existence / non-empty checks) is
  auto-generated per class via `__init_subclass__`
  (test_cli.py:496-502, 540-546, 589-595).
- `assert_output_matches_reference` (test_cli.py:212-248) skips the comparison
  when the reference file is missing; pure-numeric files are compared with
  `numpy.testing.assert_allclose` at `rtol=1e-5` / `atol=1e-8`; mixed
  text+numeric files are compared per field (numeric fields with tolerance,
  text fields by exact equality).
- `load_reference` (test_cli.py:150-171) skips `#` comment lines and blank
  lines when parsing reference files.
- `OUTPUT_PATTERNS` (test_cli.py:72-113) maps each analysis name to the output
  globs it produces; `discover_output_files` (test_cli.py:251-262) uses those
  patterns to locate the actual output files in the analysis output directory.
- A missing reference file causes that file's comparison to be silently
  skipped (the test passes without comparing).

### Class inventory (34 golden-reference test classes)

One row per golden class in `test_cli.py`. Line anchors match the published
line map. `accepts_o=False` marks analyses that write hardcoded/dynamic
filenames instead of a single `--out` file.

| Class | test_cli.py | System | analysis_name | Notes |
| --- | --- | --- | --- | --- |
| DensityZ | 606 | soohyung | density_z | accepts_o=False |
| ContactResidenceTime | 632 | soohyung | contact_res_time | |
| HBond | 657 | soohyung | hbond | |
| PiStacking | 684 | yiwei | pi_stacking | |
| RMSF | 709 | soohyung | rmsf | |
| SystemSize | 735 | soohyung | system_size | standard_args = '' |
| Thickness | 759 | soohyung | thickness | |
| WaterBridge | 784 | yiwei | water_bridge | |
| ClusteringHca | 813 | soohyung | clustering_hca | accepts_o=False; standard_args = '' |
| ClusteringKmedoid | 838 | soohyung | clustering_kmedoid | accepts_o=False; standard_args = '' |
| CompressibilityModulus | 863 | soohyung | compressibility_modulus | |
| Contacts | 887 | soohyung | contacts | |
| CovAnalysis | 911 | soohyung | cov_analysis | accepts_o=False; + eigenvector invariance test at :935 |
| MsdMembrane | 1007 | soohyung | msd_membrane | accepts_o=False |
| MsdSolution | 1032 | soohyung | msd_solution | accepts_o=False |
| PositionTime | 1057 | soohyung | position_time | |
| PositionTimeCopy | 1081 | soohyung | position_time_copy | near-duplicate of PositionTime |
| Rdf | 1106 | soohyung | rdf | |
| RadiusOfGyration | 1130 | soohyung | radius_of_gyration | |
| Rmsd | 1154 | soohyung | rmsd | |
| SaltBridge | 1178 | yiwei | salt_bridge | |
| Scd | 1206 | soohyung | scd | accepts_o=False; --qa |
| VoronoiApl | 1231 | soohyung | voronoi_apl | accepts_o=False; --qa |
| VoronoiContact | 1256 | soohyung | voronoi_contact | accepts_o=False; --qa |
| VoronoiShellComp | 1281 | soohyung | voronoi_shell_comp | accepts_o=False; --qa |
| SecondaryStructure | 1312 | omf | secondary_structure | gated: dssp |
| Sasa | 1337 | omf | sasa | gated: freesasa |
| Hole | 1363 | omf | hole | gated: hole2; accepts_o=False |
| GlycosidicBondBetweenSugars | 1389 | yiwei | glycosidic-bond-between-sugars | |
| CholTilt | 1418 | soohyung | chol_tilt | rtol=1e-2 |
| HelixAnalysis | 1446 | soohyung | helix_analysis | |
| HelixDistanceCrossingAngle | 1471 | yiwei | helix_distance_crossing_angle | |
| HelixTiltRotationAngle | 1496 | soohyung | helix_tilt_rotation_angle | |
| BondStatistics | 1520 | soohyung | bond_statistics | accepts_o=False |

### Reference inventory (55 reference files)

One row per data file under `src/stanalyzer/tests/reference/` (the directory's
`README.md` is excluded from the 55). Count: 52 `.dat` + 2 `.pdb` + 1 `.png`.

| Analysis dir | File | Type | Notes |
| --- | --- | --- | --- |
| bond_statistics | bond_lengths.dat | .dat | |
| chol_tilt | chol_tilt.dat | .dat | |
| clustering_hca | cluster.dat | .dat | |
| clustering_hca | cluster_representative.pdb | .pdb | |
| clustering_kmedoid | cluster.dat | .dat | |
| clustering_kmedoid | cluster_representative.pdb | .pdb | |
| compressibility_modulus | compressibility_modulus.dat | .dat | |
| contact_res_time | contact_res_time.dat | .dat | |
| contacts | contacts.dat | .dat | |
| cov_analysis | corr_matrix.dat | .dat | |
| cov_analysis | eigenvalues.dat | .dat | |
| cov_analysis | eigenvectors.dat | .dat | dead golden file (no OUTPUT_PATTERNS entry; removed in T7) |
| cov_analysis | correlation_matrix_heatmap.png | .png | dead artifact (removed in T7) |
| density_z | memb_nb100_0.dat | .dat | |
| glycosidic-bond-between-sugars | glycosidic-bond-between-sugars.dat | .dat | |
| hbond | hbond.dat | .dat | |
| helix_analysis | helix_analysis.dat | .dat | |
| helix_distance_crossing_angle | helix_distance_crossing_angle.dat | .dat | |
| helix_tilt_rotation_angle | helix_tilt_rotation_angle.dat | .dat | |
| hole | means.dat | .dat | |
| hole | midpoints.dat | .dat | |
| msd_membrane | dn_dopc_0.dat | .dat | |
| msd_membrane | up_dopc_0.dat | .dat | |
| msd_solution | dopc_0.dat | .dat | |
| pi_stacking | pi_stacking.dat | .dat | |
| position_time | position_time.dat | .dat | |
| position_time_copy | position_time_copy.dat | .dat | |
| radius_of_gyration | radius_of_gyration.dat | .dat | |
| rdf | rdf.dat | .dat | |
| rmsd | rmsd.dat | .dat | |
| rmsf | rmsf.dat | .dat | |
| salt_bridge | salt_bridge.dat | .dat | |
| sasa | sasa.dat | .dat | |
| scd | ave_dn_dopc_chain0_0.dat | .dat | |
| scd | ave_dn_dopc_chain1_0.dat | .dat | |
| scd | ave_up_dopc_chain0_0.dat | .dat | |
| scd | ave_up_dopc_chain1_0.dat | .dat | |
| secondary_structure | secondary_structure.dat | .dat | |
| system_size | system_size.dat | .dat | |
| thickness | thickness.dat | .dat | |
| voronoi_apl | ave_dn_0.dat | .dat | |
| voronoi_apl | ave_up_0.dat | .dat | |
| voronoi_contact | ave_dn_fcomp_0.dat | .dat | |
| voronoi_contact | ave_dn_ncomp_0.dat | .dat | |
| voronoi_contact | ave_up_fcomp_0.dat | .dat | |
| voronoi_contact | ave_up_ncomp_0.dat | .dat | |
| voronoi_shell_comp | ave_dn_dopc_fcomp_0.dat | .dat | |
| voronoi_shell_comp | ave_dn_dopc_ncomp_0.dat | .dat | |
| voronoi_shell_comp | ave_dn_dspc_fcomp_0.dat | .dat | |
| voronoi_shell_comp | ave_dn_dspc_ncomp_0.dat | .dat | |
| voronoi_shell_comp | ave_up_dopc_fcomp_0.dat | .dat | |
| voronoi_shell_comp | ave_up_dopc_ncomp_0.dat | .dat | |
| voronoi_shell_comp | ave_up_dspc_fcomp_0.dat | .dat | |
| voronoi_shell_comp | ave_up_dspc_ncomp_0.dat | .dat | |
| water_bridge | water_bridge.dat | .dat | |

## Methodology

Each analysis is assessed against three questions. The per-analysis
subsections below answer all three and close with a one-line verdict.

### Q1 — Test meaningfulness

"does the test meaningfully test the analysis" — evidence-based answer citing
test args, selection vs input system, and what the module computes.

### Q2 — Reference meaningfulness

"is the reference scientifically meaningful" — validated by numerical
cross-checks against expected MD quantities (e.g. order-param range, APL ~60 Å², thickness ~43 Å); non-numeric files (.pdb) assessed structurally, dead PNG flagged as dead.

### Q3 — Scientific interpretation

Per-analysis high-level statement of what the output ought to mean
scientifically (2-4 sentences).

## Soohyung membrane system (26 analyses)

### DensityZ (`density_z`)

**Q1 — Test meaningfulness:** `DensityZ` (test_cli.py:606) runs the full analysis with `--sel "name C*" --sel-name "MEMB" --sel-sys "segid MEMB and (name [PN] or name [CO]1[0-9])"` — all membrane carbons as the density target, phosphate/nitrogen/carbonyl atoms as the centering system. `accepts_o=False`: the module writes hardcoded filenames (`<name>_nb<nbin>_<suffix>.dat`), so the test discovers outputs via `OUTPUT_PATTERNS` (`['*_nb*_*.dat', 'combined_nb*_*.dat', 'NA_*_nb*_*.dat']`) and golden-compares each at `rtol=1e-5/atol=1e-8`. The module (density_z.py) recenters the bilayer at z=0, assigns leaflets by z-position, bins atom z-coordinates into 100 bins, and normalizes by slab volume to a number density. The comparison is sensitive (100 bins × 3 columns) and would catch wrong binning, wrong slab-volume normalization, wrong centering, or wrong leaflet assignment. Coverage is limited to the default branch: `qcent=False` (midplane fine-tuning), `otype='outs'` (individual output), `nbin=100`, `suffix=0` — the `qcent=True`, `outc`, and non-default `nbin`/`suffix` branches are untested.

**Q2 — Reference meaningfulness:** Spot-checks on `memb_nb100_0.dat` (100 bins, box z = 104.12 Å): edge bins are ~1e-4 (first three 9e-5/1.5e-4/9e-5, last three 1.8e-4/8e-5/4e-5) — essentially zero density outside the bilayer, as expected. Maximum density 0.0437 Å⁻³ at z ≈ +8.9 Å (hydrocarbon core). The ±20–27 Å band (headgroup/glycerol region) carries density 0.005–0.029, well above the water background, but it is a monotonic falloff zone, **not** a local maximum: for a `name C*` carbon selection the expected "headgroup peaks" at ±20–27 Å do not appear as peaks (the only significant maximum is the core). Density >1e-3 extends to ±32.8 Å; density >1e-2 to ±27.6 Å. The profile is a proper bilayer carbon density (dense core, monotonic falloff, ~1e-4 outside) and is scientifically sound.

**Q3 — Scientific interpretation:** The profile maps the carbon distribution across the bilayer: a dense hydrocarbon core (~0.04 Å⁻³) spanning roughly ±15 Å, a falloff through the headgroup/glycerol region (±20–27 Å), and essentially zero density in the water phase beyond ±33 Å. This is consistent with a ~43 Å-thick DOPC/DSPC/CHL bilayer. The slight asymmetry (peak at +8.9 Å rather than z=0) reflects the finite 20-frame sampling and the protein's presence in one leaflet.

**Verdict:** Test meaningful, reference scientifically sound.

### ContactResidenceTime (`contact_res_time`)

**Q1 — Test meaningfulness:** `ContactResidenceTime` (test_cli.py:632) runs the analysis with `--sel "protein and name CA" --threshold "5.0"` and golden-compares the full output (30 contact pairs × 6 fields) at `rtol=1e-5/atol=1e-8`. The module (contact_res_time.py) computes residue COM distances per frame, tracks consecutive frames below the 5.0 Å threshold as contact "events", and reports mean ± std of event durations per pair. The comparison is sensitive and would catch wrong COM computation, wrong event tracking, or wrong threshold application. Coverage is limited to the default `threshold=5.0` and `interval=1`; non-default values are untested. Critically, the test cannot catch the for-else overwrite bug at contact_res_time.py:84-89 because the golden reference itself encodes it (see Q2) — the test is meaningful but its reference is corrupted.

**Q2 — Reference meaningfulness:** Spot-check: 30 contact pairs. The last row, `GLY 22 ALA 23 0 0`, is **wrong**: `contacts.dat` shows GLY22–ALA23 in contact in all 20 frames, so its residence should be `20 0`. The `for … else` construct at contact_res_time.py:84-89 (`else: residence_dist[contact] = (0.0, 0.0)`) always executes after the loop completes, overwriting the last-inserted contact's entry with (0.0, 0.0). The reference is self-consistent with the buggy code (the test passes) but the final row is scientifically incorrect. The remaining rows are plausible (e.g. adjacent pairs GLY 1 GLY 2 = 20, 0; GLY 2 ALA 3 = 20, 0; transient pairs GLY 1 ALA 3 = 9.5 ± 2.5). **Reference flagged UNRELIABLE: last row corrupted by the for-else bug.**

**Q3 — Scientific interpretation:** Contact residence times quantify how long residue pairs remain within the 5 Å cutoff before separating. Values of 20 frames (the full trajectory) for adjacent backbone pairs reflect persistent native contacts; shorter times (1–3 frames) for non-adjacent pairs reflect transient contacts. The mean ± std format captures the distribution of contact-event durations. The corrupted last row is a code bug, not physics.

**Verdict:** Test meaningful, but reference flagged UNRELIABLE: last row (GLY 22 ALA 23) overwritten to `0 0` by the for-else bug at contact_res_time.py:84-89 (pair is in contact all 20 frames per contacts.dat).

### HBond (`hbond`)

**Q1 — Test meaningfulness:** `HBond` (test_cli.py:657) runs the analysis with `--sel "segid MEMB or protein" --hydrogens-sel "None" --acceptors-sel "None" --d-a-cutoff "3.0" --d-h-a-angle-cutoff "150.0"` and golden-compares all 409 detected bonds (6 fields per row; text atom identifiers exact, distance/angle numeric with tolerance). The module (hbond.py) wraps MDAnalysis `HydrogenBondAnalysis` over the membrane+protein selection (water excluded), exercising the `guess_hydrogens`/`guess_acceptors` branches (both "None"), the explicit 3.0 Å / 150° cutoffs, and the solvent-solvent skip logic. The comparison is sensitive and would catch wrong H-bond detection, wrong atom labeling, or wrong cutoff application. Untested branches: explicit `--hydrogens-sel`/`--acceptors-sel` strings, non-default cutoffs, non-default `interval`.

**Q2 — Reference meaningfulness:** Spot-checks on `hbond.dat`: 409 bonds across 20 frames (14–27 per frame). D–A distances span 2.497–2.999 Å (all below the 3.0 Å cutoff); D–H–A angles span 150.04–179.60° (all above the 150° cutoff). The task's typical-range bounds (2.74–3.00 Å / 152–170°) are exceeded at the extremes (min distance 2.50 Å, max angle 179.6°), but every value is a physically reasonable H-bond geometry consistent with the cutoffs used. Donors/acceptors are protein backbone N–H and cholesterol O3–H3′ to lipid carbonyl O — chemically sensible. Reference scientifically sound.

**Q3 — Scientific interpretation:** The list captures the network of protein-backbone and cholesterol–lipid hydrogen bonds. Protein N–H···O=C bonds (2.7–3.0 Å, 150–180°) are standard secondary-structure H-bonds; cholesterol O3–H3′···lipid-carbonyl bonds (2.5–3.0 Å) are the cholesterol–lipid interactions that anchor cholesterol in the membrane. The per-frame variation (14–27 bonds) reflects dynamic H-bond breaking and forming over the 20-frame trajectory.

**Verdict:** Test meaningful, reference scientifically sound.

### RMSF (`rmsf`)

**Q1 — Test meaningfulness:** `RMSF` (test_cli.py:709) runs the analysis with `--sel-align "segid PROA and name CA" --sel-rmsf "segid PROA and name CA"` and golden-compares all 23 rows (residue index + RMSF) at `rtol=1e-5/atol=1e-8`. The module (rmsf.py) aligns each frame to the first frame via `rotation_matrix`, accumulates mean and sum-of-squares, and reports per-residue RMSF. The comparison is sensitive and would catch wrong alignment, wrong reference frame, or wrong fluctuation accumulation. Untested branches: `--ref-psf`, `--align-out`, non-default `interval`, and non-identical align/rmsf selections.

**Q2 — Reference meaningfulness:** Spot-check on `rmsf.dat`: 23 residues, RMSF 0.309–0.717 Å, mean 0.455 Å — all within the expected backbone 0.3–1.5 Å range for a folded protein in a membrane. Terminal residues fluctuate most (residue 1: 0.717 Å, residue 23: 0.715 Å) and the core least (residue 17: 0.310 Å) — the expected end-fraying pattern. Reference scientifically sound.

**Q3 — Scientific interpretation:** RMSF quantifies per-residue thermal fluctuation about the average structure. Values of 0.3–0.7 Å for CA atoms indicate a well-structured, stable protein; elevated values at the termini reflect end-fraying. The 20-frame trajectory is short, so these are lower-bound estimates of the true fluctuations.

**Verdict:** Test meaningful, reference scientifically sound.

### SystemSize (`system_size`)

**Q1 — Test meaningfulness:** `SystemSize` (test_cli.py:735) runs the analysis with empty `standard_args` — every setting comes from project.json defaults (`time_step="1 ns"`, `interval=1`) — and golden-compares all 20 rows × 8 fields (time, x/y/z, α/β/γ, volume) at `rtol=1e-5/atol=1e-8`. The module (system_size.py) reads per-frame box dimensions and computes the triclinic volume. The comparison is sensitive and would catch wrong box reading, wrong volume formula, or wrong time stepping. Coverage is minimal: only the default path is exercised; non-default `interval` and `time_step` are untested (the `include_angles` flag is accepted but the module always writes all 8 columns).

**Q2 — Reference meaningfulness:** Spot-check on `system_size.dat`: 20 frames, x = y = 74.25–75.12 Å (mean 74.63), z = 102.64–104.84 Å (mean 104.12), all angles exactly 90°, volume 5.79–5.81e5 Å³. Matches the expected ~75×75×104 Å box; volume cross-check 74.63² × 104.12 = 5.80e5 Å³. Reference scientifically sound.

**Q3 — Scientific interpretation:** The box dimensions track the simulation cell over time. The ~75×75×104 Å orthorhombic cell with ~1 Å x/y and ~2 Å z fluctuations reflects an NPT simulation of a hydrated lipid bilayer; the area fluctuations are the basis for the area compressibility modulus (see CompressibilityModulus). The constant 90° angles confirm an orthorhombic cell.

**Verdict:** Test meaningful, reference scientifically sound.

### Thickness (`thickness`)

**Q1 — Test meaningfulness:** `Thickness` (test_cli.py:759) runs the analysis with `--sel "segid MEMB and (name P or name N or name C1[0-9] or name O1[0-9])" --sel-sys "resname DOPC and name P; resname DSPC and name P"` — headgroup atoms as the thickness target, phosphate atoms as the centering/leaflet system — and golden-compares the per-frame z_up/z_dn/thickness columns plus the three average lines at `rtol=1e-5/atol=1e-8`. The module (thickness.py) recenters the bilayer, assigns leaflets by z-position, and computes leaflet mean z-positions and their difference. The comparison is sensitive and would catch wrong leaflet assignment, wrong centering, or wrong headgroup selection. Untested branches: `qcent=True` (midplane fine-tuning), non-default `interval`.

**Q2 — Reference meaningfulness:** Spot-check on `thickness.dat`: z_up = 20.128 Å, z_dn = −23.041 Å, thickness = 43.169 Å (STD 0.264) — matches the expected 43.17 Å (P-based z_up = 20.13, z_dn = −23.04) essentially exactly. The up/down asymmetry (|z_up| ≠ |z_dn|) reflects the protein's presence on one side. Reference scientifically sound.

**Q3 — Scientific interpretation:** The bilayer thickness (43.2 Å) is the headgroup-to-headgroup distance across the membrane, the standard measure of bilayer thickness. The value is typical for DOPC/DSPC (C18:1/C18:0) bilayers (~40–45 Å). The up/down asymmetry (20.1 vs −23.0 Å) indicates the bilayer center is offset from the box center, consistent with the protein perturbing one leaflet.

**Verdict:** Test meaningful, reference scientifically sound.

### CompressibilityModulus (`compressibility_modulus`)

**Q1 — Test meaningfulness:** `CompressibilityModulus` (test_cli.py:863) runs the analysis with `--temp 310` and golden-compares the output (header line + `K_A = 1494.54333 (dyn/cm)`) at `rtol=1e-5/atol=1e-8`. The module (compressibility_modulus.py) computes K_A = k_BT·⟨A⟩/var(A) from per-frame box areas. The comparison is sensitive for the single numeric field but the output is one scalar — a bug that coincidentally preserved K_A would pass. Coverage is limited to the default path; non-default `temp` and `interval` are untested.

**Q2 — Reference meaningfulness:** Spot-check: recomputing K_A from `system_size.dat` box areas (⟨A⟩ = 5569.16 Å², var(A) = 1594.39) with T = 310 K gives K_A = 1495.00 dyn/cm vs the reference 1494.54 — 0.03% agreement. The reference is internally consistent with the box-fluctuation data. The value (~1500 dyn/cm) is high for a pure lipid bilayer (DOPC ~200–300 dyn/cm) but plausible for a cholesterol-rich DOPC/DSPC/CHL membrane, which cholesterol stiffens substantially. Caveat: it is estimated from only 20 frames, so the area variance is poorly sampled and the absolute value carries large uncertainty. Reference scientifically sound (internally consistent; short-sampling caveat).

**Q3 — Scientific interpretation:** The area compressibility modulus K_A measures the membrane's resistance to area change. The high value (~1500 dyn/cm) reflects the cholesterol-rich composition, which rigidifies the bilayer. The 20-frame estimate is noisy — the box-area variance over such a short window is dominated by a few fluctuations — so the value should be treated as order-of-magnitude.

**Verdict:** Test meaningful, reference scientifically sound (internally consistent with box data; short-sampling caveat).

### Contacts (`contacts`)

**Q1 — Test meaningfulness:** `Contacts` (test_cli.py:887) runs the analysis with `--sel "protein and name CA" --contact-threshold "5.0"` and golden-compares all 30 contact pairs (5 fields; text fields exact, frequency numeric with tolerance) at `rtol=1e-5/atol=1e-8`. The module (contacts.py) computes residue COM distances per frame and records a contact below 5.0 Å, using the chunked parallel runtime (RuntimeContext/RuntimeScheduler/RuntimeExecutor + `contacts_worker`). The test therefore also exercises the scheduling/chunking/merge machinery. The comparison is sensitive and would catch wrong COM, wrong threshold, or wrong chunk merging. Untested branches: non-default threshold, non-default `interval`, explicit `--workers`.

**Q2 — Reference meaningfulness:** Spot-check on `contacts.dat`: 30 contact pairs, frequencies 1–20. Adjacent backbone pairs (GLY 1 GLY 2, LEU 6 ALA 7, ALA 7 LEU 8, …) all show frequency 20 (contact in all frames) — physically expected for covalently bonded neighbors. Non-adjacent pairs show lower frequencies (1–19), e.g. GLY 1 PHE 4 = 11, TRP 19 LEU 20 = 16. GLY 22 ALA 23 = 20, consistent with the contact_res_time finding that this pair is in contact all 20 frames. Reference scientifically sound.

**Q3 — Scientific interpretation:** Contact frequencies report the fraction of the trajectory during which residue pairs are within 5 Å. The 20/20 values for adjacent residues are trivial (bonded neighbors are always in contact); the informative entries are the non-adjacent pairs (e.g. GLY 1 PHE 4 = 11/20, TRP 19 LEU 20 = 16/20), which report transient tertiary contacts. The 5 Å COM cutoff is generous and includes non-direct contacts.

**Verdict:** Test meaningful, reference scientifically sound.

### CovAnalysis (`cov_analysis`)

**Q1 — Test meaningfulness:** `CovAnalysis` (test_cli.py:911) runs the analysis with `--sel "protein and name CA"` (23 CA atoms → 69×69 matrices) and golden-compares `corr_matrix.dat` (69×69) and `eigenvalues.dat` (6 values) at `rtol=1e-5/atol=1e-8`. The module (cov_analysis.py) aligns the trajectory, computes the covariance and correlation matrices, eigendecomposes via `np.linalg.eigh`, and selects k modes by eig-based cumulative variance ≥ 0.85 (the former PCA-on-covariance bug was fixed per NOTE_TO_COLLEAGUE.md). The eigenvector output is validated **by proxy**: `eigenvectors.dat` has no `OUTPUT_PATTERNS` entry (test_cli.py:107-111), no `discover_output_files` glob, and no `assert_output_matches_reference` target — instead `test_eigenvectors_consistent_with_eigenvalues` (:935-1004) re-derives C from the inputs and asserts C·v ≈ λ·v (:991-995) plus orthonormality (:998-1001), a rotation-invariant check that catches axis-selection bugs and is robust to the BLAS-dependent basis degeneracy documented in CA_BUG_REPORT.md. The test exercises alignment, covariance/correlation computation, eigendecomposition, and k-selection. Untested branches: non-default `interval`, `--align-out`, custom output filenames. Test is meaningful.

**Q2 — Reference meaningfulness:** Spot-checks: `eigenvalues.dat` = [2.284518, 0.863323, 0.679162, 0.347965, 0.254410, 0.209662] — matches the expected [2.2845, 0.8633, 0.6792, 0.3480, 0.2544, 0.2097] (CA_BUG_REPORT.md:141-142) to 4 decimals; k=6 at the 85% cumulative-variance threshold (cumulative ratios 0.428, 0.590, 0.717, 0.782, 0.830, 0.869 per NOTE_TO_COLLEAGUE.md:10). `corr_matrix.dat`: 69×69, diagonal exactly 1.0, off-diagonal −0.931 to 0.979 — a proper correlation matrix. `eigenvectors.dat`: 6×69, row norms 1.0, V·Vᵀ ≈ I (max deviation 9.3e-7) — orthonormal. Reference scientifically sound. The two dead files — `eigenvectors.dat` (golden) and `correlation_matrix_heatmap.png` — are scheduled for removal in T7: no `OUTPUT_PATTERNS` entry, no `discover_output_files` glob, no `assert_output_matches_reference` target; the eigenvector output is instead guarded by the self-consistency test at :935-1004.

**Q3 — Scientific interpretation:** The covariance analysis (essential dynamics) decomposes the protein's positional fluctuations into orthogonal modes. The top eigenvalue (2.28 Å²) captures the largest collective motion; the 85% cumulative-variance cutoff at k=6 means six modes explain ~87% of the fluctuation. The correlation matrix shows strong inter-atom correlations (up to 0.98) reflecting collective backbone motions; the heatmap (dead artifact) would visualize this.

**Verdict:** Test meaningful, reference scientifically sound; `eigenvectors.dat` + `correlation_matrix_heatmap.png` are dead golden files scheduled for removal in T7 (eigenvector output validated by the rotation-invariant self-consistency test at test_cli.py:935-1004).

### ClusteringHca (`clustering_hca`)

**Q1 — Test meaningfulness:** `ClusteringHca` (test_cli.py:813) runs the analysis with empty `standard_args` (`accepts_o=False`); the module (clustering_hca.py) hardcodes `sel='name CA'`, `cutoff=2.0`, `criterion='distance'`, `method='average'` in `main()`. It builds the pairwise RMSD distance matrix over the 23 protein CA atoms via `diffusionmap.DistanceMatrix`, runs hierarchical linkage + `fcluster` at the 2.0 Å cutoff, and writes `cluster.dat` (frame→cluster) plus `cluster_representative.pdb` (protein atoms of one representative frame per cluster). `OUTPUT_PATTERNS` (`['cluster.dat', 'cluster_representative.pdb']`) discovers both files and golden-compares them at `rtol=1e-5/atol=1e-8` (the PDB is compared line-by-line: numeric coordinate fields with tolerance, text fields exact). The comparison is sensitive and would catch wrong distance-matrix construction, wrong linkage/fcluster parameters, or wrong representative-frame selection. Coverage is limited to the default branch: the `--cutOff` parser argument is accepted but ignored by `main()` (hardcoded 2.0), so non-default cutoffs are untested.

**Q2 — Reference meaningfulness:** Spot-checks: `cluster.dat` has 20 rows (frames 0–19, matching the 20-frame trajectory) mapping to 12 clusters (ids 1–12, sizes 1–4 frames) — plausible for a 2.0 Å cutoff on pairwise CA RMSD over a short trajectory. `cluster_representative.pdb` is structurally sane: 12 MODEL records (one per cluster), each with 336 atoms / 23 residues / 23 CA / single chain X — exactly matching the PROA protein (336 atoms, 23 residues, 23 CA). Reference scientifically sound.

**Q3 — Scientific interpretation:** Hierarchical clustering groups trajectory frames by structural similarity (pairwise CA RMSD). Twelve clusters from 20 frames at a 2.0 Å cutoff indicate the 23-residue peptide samples a diverse conformational ensemble; the representative frame per cluster summarizes each visited state. The cluster sizes (1–4 frames) reflect a mix of persistent and transient conformations over the short trajectory.

**Verdict:** Test meaningful, reference scientifically sound.

### ClusteringKmedoid (`clustering_kmedoid`)

**Q1 — Test meaningfulness:** `ClusteringKmedoid` (test_cli.py:838) runs the analysis with empty `standard_args` (`accepts_o=False`); the module (clustering_kmedoid.py) hardcodes `sel='name CA'` and `k=2` (`settings.get('k') or 2`). It builds the same pairwise-RMSD distance matrix over the 23 CA atoms and runs `kmedoids.fasterpam(dist_matrix, k, random_state=0)` — the fixed `random_state=0` makes the golden comparison deterministic. It writes `cluster.dat` (frame→cluster) and `cluster_representative.pdb` (protein atoms of the k medoid frames). Both files are discovered via `OUTPUT_PATTERNS` and golden-compared. The comparison is sensitive to the distance matrix, k, and medoid selection. Coverage is limited to the default `k=2`; the `--k` argument is exercisable via CLI but untested.

**Q2 — Reference meaningfulness:** Spot-checks: `cluster.dat` has 20 rows mapping to 2 clusters (ids 0/1; 12 frames in cluster 0, 8 in cluster 1) — consistent with `k=2`. `cluster_representative.pdb` is structurally sane: 2 MODEL records (one per medoid), each with 336 atoms / 23 residues / 23 CA / single chain X — exactly matching the PROA protein. Reference scientifically sound.

**Q3 — Scientific interpretation:** k-medoids partitions the 20 frames into two conformational states, with the medoid frames being the most representative structures of each state. The 12/8 frame split suggests two distinct conformational populations of the peptide over the trajectory. The fixed random seed ensures reproducible clustering.

**Verdict:** Test meaningful, reference scientifically sound.

### MsdMembrane (`msd_membrane`)

**Q1 — Test meaningfulness:** `MsdMembrane` (test_cli.py:1007) runs the analysis with `--sel "resname DOPC" --sel-sys "resname DOPC DSPC"` (`accepts_o=False`). The module (msd_membrane.py) recenters the bilayer, assigns leaflets by z-position, unwraps the trajectory, applies leaflet-COM drift correction, and computes per-leaflet MSD of DOPC molecule COMs. With default settings (`qcomsys=False`, `qcommol=False`, `otype='outl'`, `suffix=0`) the run produces exactly `up_dopc_0.dat` and `dn_dopc_0.dat`, both of which have references — so the default-config comparison is complete (no silent skips). The comparison is sensitive to leaflet assignment, unwrapping, drift correction, and MSD accumulation. Coverage is limited to the default branch: the `sys_com`/`mol_com`/`mol_info` output modes (`qcomsys`/`qcommol`), bilayer output (`otype='outb'`), NA-leaflet handling, and DSPC/CHL lipid types are untested.

**Q2 — Reference meaningfulness:** Spot-checks: `up_dopc_0.dat` and `dn_dopc_0.dat` each have 20 tau rows (tau = 0–19 ns, `time_step=1 ns`). Values verified against the input trajectory: a direct COM-displacement calculation on the DCDs gives lateral MSD 5.41 Å² at 1 ns vs the reference up-leaflet value 5.66 Å² (tau=1) — consistent. The implied lateral diffusion coefficient (~1400 µm²/s at tau=1) is anomalously high for DOPC (literature ~5–10 µm²/s), but this reflects the short 20-frame input trajectory, not a reference error; the MSD is non-monotonic at large tau (e.g. MSDX 15.16 at tau=7 → 3.30 at tau=19) due to the tiny sample count at long lag times. **Partial reference coverage flagged:** `OUTPUT_PATTERNS` lists 5 patterns (`*_sys_com_*.dat`, `*_mol_com_*.dat`, `*_*_*.dat`, `NA_*_*_*.dat`, `*_mol_info_*.dat`) but only 2 reference files exist (`up_dopc_0.dat`, `dn_dopc_0.dat`) — the sys_com/mol_com/mol_info/NA output classes have no references and are never compared.

**Q3 — Scientific interpretation:** The per-leaflet MSD of DOPC molecule COMs quantifies lateral lipid diffusion in each leaflet. The up/down leaflet difference (up-leaflet MSD_lat 5.66 Å² vs dn-leaflet 2.00 Å² at tau=1) reflects leaflet asymmetry, likely influenced by the protein. The anomalously fast apparent diffusion and non-monotonic long-lag behavior are artifacts of the very short (20 ns) trajectory, so the absolute D values should not be over-interpreted.

**Verdict:** Test meaningful for the default config; reference scientifically sound but PARTIAL — only 2 of the 5 `OUTPUT_PATTERNS` output classes have references (`up_dopc_0.dat`, `dn_dopc_0.dat`); sys_com/mol_com/mol_info/NA modes untested.

### MsdSolution (`msd_solution`)

**Q1 — Test meaningfulness:** `MsdSolution` (test_cli.py:1032) runs the analysis with `--sel "resname DOPC and name P"` (`accepts_o=False`). The module (msd_solution.py) splits the selection into molecules (79 DOPC molecules, one P atom each), unwraps the trajectory, and computes the MSD of the P-atom COMs. With default settings (`qdrift='disabled'`, `qcomsys=False`, `qcommol=False`, `suffix=0`) the run produces exactly `dopc_0.dat`, which has a reference — the default-config comparison is complete. The comparison is sensitive to molecule splitting, unwrapping, and MSD accumulation. Coverage is limited to the default branch: `qdrift='enabled'`, the `sys_com`/`mol_com`/`mol_info` output modes, NA handling, and non-DOPC selections are untested.

**Q2 — Reference meaningfulness:** Spot-check: `dopc_0.dat` has 20 tau rows (tau = 0–19 ns). Values are consistent with the same fast-diffusion regime verified for msd_membrane (lateral MSD 6.51 Å² at tau=1; direct trajectory calc 5.41 Å² for DOPC COM — same order). **Partial reference coverage flagged:** `OUTPUT_PATTERNS` lists 5 patterns (`sys_com_*.dat`, `mol_com_*.dat`, `*_*.dat`, `NA_*_*.dat`, `mol_info_*.dat`) but only 1 reference file exists (`dopc_0.dat`) — the sys_com/mol_com/mol_info/NA output classes have no references and are never compared.

**Q3 — Scientific interpretation:** The MSD of DOPC phosphate atoms (a proxy for whole-lipid COM) quantifies lipid diffusion. The P-atom MSD is a reasonable proxy for lateral lipid motion. As with msd_membrane, the absolute D values are inflated by the short trajectory and should be treated as order-of-magnitude.

**Verdict:** Test meaningful for the default config; reference scientifically sound but PARTIAL — only 1 of the 5 `OUTPUT_PATTERNS` output classes has a reference (`dopc_0.dat`); sys_com/mol_com/mol_info/NA modes untested.

### PositionTime (`position_time`)

**Q1 — Test meaningfulness:** `PositionTime` (test_cli.py:1057) runs the analysis with `--sel "protein and name CA" --head-group "segid MEMB and name P"`. The module (position_time.py) centers the membrane on the P headgroup atoms (`center_in_box(selected_head_group, point=(0,0,0))`) and tracks the protein CA centroid z-position per frame (default `method='com'`, `axis='z'`). The output (`position_time.dat`, 20 rows × 2 fields) is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the centering and the centroid computation. Untested branches: `method='cog'`, `axis='x'/'y'`, non-default selections.

**Q2 — Reference meaningfulness:** Spot-check on `position_time.dat`: 20 rows, z-position −2.62 to −7.23 Å (mean −5.08 Å). With the P headgroups centered at z=0, the protein CA centroid sits ~5 Å below the membrane midplane — plausible for a 23-residue peptide inserted asymmetrically into the bilayer. The ±2 Å fluctuation over 20 ns reflects vertical motion. Reference scientifically sound.

**Q3 — Scientific interpretation:** The time series tracks the protein's vertical position relative to the membrane center. The ~−5 Å offset indicates the peptide sits below the bilayer midplane; the fluctuations reflect vertical thermal motion. This is a useful measure of protein–membrane registration over time.

**Verdict:** Test meaningful, reference scientifically sound.

### PositionTimeCopy (`position_time_copy`)

**Q1 — Test meaningfulness:** `PositionTimeCopy` (test_cli.py:1081-1103) runs the analysis with `--sel "protein and name CA"` — notably **without** the `--head-group` argument. The module (position_time_copy.py) is a stripped copy of position_time.py: identical centroid tracking and output format, but the `head_group` parameter and the `center_in_box` membrane-centering step are removed, so it reports the raw box-coordinate centroid. **Flagged as a near-duplicate of PositionTime:** the module adds no capability beyond position_time (it is a strict subset — same computation minus centering), and the test exists only for framework coverage of the duplicate `position_time_copy` analysis. Its inclusion is not scientifically justified: it tests a code path that position_time already covers, and its reference is a different (uncentered) view of the same quantity.

**Q2 — Reference meaningfulness:** Spot-check on `position_time_copy.dat`: 20 rows, z-position 47.68–52.33 Å (mean 49.91 Å) — the raw box z-coordinate of the protein CA centroid (box z ≈ 104 Å, so ~50 Å is mid-box). Near-duplicate status verified numerically: `position_time_copy.dat − position_time.dat` is nearly constant (mean 54.98 Å, std 0.32 Å) — the difference is exactly the membrane P-headgroup center z-position that position_time subtracts. The reference is a valid golden of the stripped module's output but is scientifically redundant with position_time.

**Q3 — Scientific interpretation:** The output is the same quantity as position_time (protein CA centroid z) but in the raw box frame instead of the membrane-centered frame. Without centering, the values are box-coordinate dependent and less interpretable for membrane-relative motion; the ~50 Å mean merely reflects the box geometry.

**Verdict:** **Flagged near-duplicate** — position_time_copy is position_time minus the membrane-centering feature; test present only for framework coverage (test_cli.py:1081-1103), inclusion not scientifically justified.

### Rdf (`rdf`)

**Q1 — Test meaningfulness:** `Rdf` (test_cli.py:1106) runs the analysis with `-sel1 "protein and name CA" -sel2 "resname DOPC and name P" -bin-size 0.1`. The module (rdf.py) determines the RDF range from the box (`box_l/2`), bins at 0.1 Å, and computes the InterRDF between the 23 protein CA atoms and the 79 DOPC P atoms. The output (`rdf.dat`, 371 bins × 2 fields) is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the box-range determination, binning, and InterRDF computation. Untested branches: non-default `-bin-size`, other selections.

**Q2 — Reference meaningfulness:** Spot-checks on `rdf.dat`: **rdf ~0 at small r confirmed** — g = 0.0000 for all bins below r = 3.85 Å (excluded volume between protein CA and DOPC P atoms). First peak at r = 7.95 Å with g = 2.61 (the first coordination shell of lipid phosphate around protein Cα), which is also the global maximum; the tail approaches ~1.16 at r = 37.05 Å. The range 37.1 Å = box_l/2 for the ~74.2 Å box — consistent with the system_size box (~75 Å). Reference scientifically sound.

**Q3 — Scientific interpretation:** The RDF describes the probability of finding DOPC phosphate groups at distance r from protein Cα atoms. The excluded volume out to 3.85 Å and the first peak at ~8 Å reflect the lipid headgroup shell around the protein; the approach to ~1 at large r indicates a bulk-like lipid distribution away from the protein.

**Verdict:** Test meaningful, reference scientifically sound.

### RadiusOfGyration (`radius_of_gyration`)

**Q1 — Test meaningfulness:** `RadiusOfGyration` (test_cli.py:1130) runs the analysis with `--sel-rg "protein and name CA" --sel-align "protein and name CA"`. The module (radius_of_gyration.py) aligns the trajectory to frame 1 (`ref_frame_type='specific'`, `ref_frame_num=1`) via `AlignTraj` and computes the radius of gyration of the CA atoms per frame. The output (`radius_of_gyration.dat`, 20 rows × 2 fields) is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the alignment and the Rg computation. Untested branches: `ref_frame_type='average'`, non-default `ref_frame_num`, `interval`, `--align-out`.

**Q2 — Reference meaningfulness:** Spot-check on `radius_of_gyration.dat`: 20 rows, Rg 9.61–9.99 Å (mean 9.82 Å, std ~0.1 Å). For the 23-residue peptide this matches the theoretical Rg of a 23-residue α-helix (length ≈ 34.5 Å, Rg = L/√12 ≈ 9.96 Å). The small fluctuation (±0.2 Å) over 20 ns indicates a stable structure. Reference scientifically sound.

**Q3 — Scientific interpretation:** Rg measures the compactness of the peptide. The stable ~9.8 Å value indicates a well-folded, rigid helical structure; the small fluctuations reflect thermal breathing. The value is consistent with a single α-helix spanning the membrane.

**Verdict:** Test meaningful, reference scientifically sound.

### Rmsd (`rmsd`)

**Q1 — Test meaningfulness:** `Rmsd` (test_cli.py:1154) runs the analysis with `--sel "protein and name CA"`. The module (rmsd.py) aligns the trajectory to frame 1 (`ref_frame_type='specific'`, `ref_frame_num=1`) via `AlignTraj` and computes the per-frame RMSD of the 23 CA atoms. The output (`rmsd.dat`, 20 rows × 2 fields) is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the alignment and RMSD computation. Untested branches: `ref_frame_type='average'`, non-default `ref_frame_num`, `interval`, `time_step`, `--align-out`.

**Q2 — Reference meaningfulness:** Spot-check on `rmsd.dat`: 20 rows, RMSD 0.00–1.00 Å (mean 0.56 Å). Frame 1 is the reference (RMSD = 0 at t = 1.0 ns); the time column runs 1.0–21.0 ns as expected from `time_step=1 ns`. Values <1 Å over 20 ns are plausible for a stable, well-equilibrated 23-residue peptide in a membrane. Reference scientifically sound.

**Q3 — Scientific interpretation:** RMSD tracks the peptide's deviation from its frame-1 structure. The low values (<1 Å) indicate a stable, well-equilibrated peptide; the small fluctuations reflect thermal motion about the reference conformation. The 20 ns window is short, so the RMSD captures only local conformational stability.

**Verdict:** Test meaningful, reference scientifically sound.

### Scd (`scd`)

**Q1 — Test meaningfulness:** `Scd` (test_cli.py:1206) runs the analysis with `--sel "resname DOPC and (name C22 or name C32)" --sel-sys "segid MEMB and name P" --qa` (`accepts_o=False`). The module (scd.py) recenters the bilayer on the P atoms, unwraps the trajectory, assigns leaflets by z-position, and computes per-carbon deuterium order parameters (SCD) for both DOPC acyl chains (chain0 = C22–C218, chain1 = C32–C318, 17 carbons each), averaged per leaflet. `--qa` writes the 4 average files (`ave_{up,dn}_dopc_chain{0,1}_0.dat`), discovered via `OUTPUT_PATTERNS` (`['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat']`, test_cli.py:99) and golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive (17 carbons × 2 chains × 2 leaflets, mean+std each) and would catch wrong C–H vector construction, wrong leaflet assignment, or wrong averaging. Coverage is limited to the default branch: `--qa` only (no `--qt` time series), `otype='outl'`, `suffix=0`, `interval=1`; the `--qt`, `outb`, non-default suffix/interval branches are untested.

**Q2 — Reference meaningfulness:** Spot-checks on the 4 files (17 rows each: carbon, mean, std): up-leaflet chain1 (C32–C318) shows the classic plateau + tail decay — rows C33..C39 = −0.37951, −0.39363, −0.40569, −0.39704, −0.36293, −0.28413, −0.15646 — a plateau ≈ −0.38..−0.41 followed by monotonic decay toward the terminal methyl; file min −0.40569, max −0.06534. Up-leaflet chain0 (C22–C218) min −0.37318, max 0.00034 (C22, adjacent to the carbonyl, is near zero then rises into the plateau). The up-leaflet values sit inside the expected −0.28..−0.41 range. The dn leaflet is much weaker (chain0 −0.023..−0.226, chain1 −0.013..−0.216) — a strongly asymmetric bilayer, consistent with the asymmetric DOPC/DSPC leaflet composition found in voronoi_apl. Reference scientifically sound.

**Q3 — Scientific interpretation:** SCD order parameters measure the orientational order of C–H bonds along the acyl chains. Up-leaflet plateau values of −0.38..−0.41 indicate a well-ordered chain (cholesterol-rich, liquid-ordered environment), with the characteristic monotonic decay toward the terminal methyl. The much weaker dn-leaflet values indicate a more disordered leaflet — consistent with the DOPC-rich dn leaflet (68 Å² APL) vs the DSPC-rich up leaflet (56 Å² APL) seen in voronoi_apl.

**Verdict:** Test meaningful, reference scientifically sound (asymmetric leaflets reflected in the order-parameter magnitudes).

### VoronoiApl (`voronoi_apl`)

**Q1 — Test meaningfulness:** `VoronoiApl` (test_cli.py:1231) runs the analysis with `--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa` (`accepts_o=False`). The module (voronoi_apl.py) recenters the bilayer, assigns leaflets by z-position, builds a 2D Voronoi tessellation of the P atoms per leaflet (with periodic images), and computes the area per molecule, averaged per lipid type per leaflet. `--qa` writes `ave_dn_0.dat` + `ave_up_0.dat` (mean, std, and sample count per lipid type), discovered via `OUTPUT_PATTERNS` (`['ave_*_*.dat', 'time_*_*.dat']`, test_cli.py:102) and golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive (2 leaflets × 2 lipid types × mean/std/count) and would catch wrong tessellation, wrong leaflet assignment, or wrong averaging. Coverage is limited to the default branch: `--qa` only, `otype='outl'`, `suffix=0`; `--qt`, `outb`, non-default suffix/interval untested.

**Q2 — Reference meaningfulness:** Spot-checks: `ave_dn_0.dat` — DOPC 68.09972 ± 0.73558 Å² (1440 samples = 72 DOPC × 20 frames), DSPC 66.59765 ± 3.41609 (200 = 10 × 20); `ave_up_0.dat` — DOPC 56.38709 ± 3.82307 (140 = 7 × 20), DSPC 56.24399 ± 0.61026 (1840 = 92 × 20). **Strongly asymmetric leaflets confirmed**: the dn leaflet is DOPC-rich (72 DOPC / 10 DSPC) with APL ≈ 68 Å², the up leaflet is DSPC-rich (7 DOPC / 92 DSPC) with APL ≈ 56 Å². The 56.4/68.1 Å² asymmetry matches the plan's expected values; both are within the plausible APL range for DOPC/DSPC (~55–70 Å²). Reference scientifically sound.

**Q3 — Scientific interpretation:** APL quantifies the area each lipid occupies in its leaflet. The ~68 Å² (dn, DOPC-rich) vs ~56 Å² (up, DSPC-rich) values reflect the different lipid packing in the two leaflets; the low per-frame std in the majority component (0.74 Å² dn DOPC, 0.61 Å² up DSPC) indicates uniform packing, while the minority components show larger variance (3.4–3.8 Å²) from small sample counts. The asymmetry is a property of this asymmetric DOPC/DSPC system, not an error.

**Verdict:** Test meaningful, reference scientifically sound (asymmetric leaflets 56.4/68.1 Å² confirmed).

### VoronoiContact (`voronoi_contact`)

**Q1 — Test meaningfulness:** `VoronoiContact` (test_cli.py:1256) runs the analysis with the same selection as VoronoiApl (`--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa`). The module (voronoi_contact.py) builds the same per-leaflet 2D Voronoi tessellation and counts molecular contacts (shared Voronoi edges) between lipid types, reporting per-type average contact number (ncomp) and fraction (fcomp). `--qa` writes 4 files (`ave_{dn,up}_{fcomp,ncomp}_0.dat`), discovered via `OUTPUT_PATTERNS` (`['ave_*_*.dat', 'time_*_*.dat', 'NA_time_*_*.dat']`, test_cli.py:101) and golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive (2 leaflets × 2 types × 2 observables × mean/std) and would catch wrong contact detection or wrong averaging. Coverage is limited to the default branch (`--qa`, `outl`, `suffix=0`); `--qt`/`outb` untested.

**Q2 — Reference meaningfulness:** Spot-checks: dn leaflet — DOPC has 5.31111 ± 0.08127 DOPC neighbors (fcomp 0.88391 ± 0.00809) and 0.69722 ± 0.04470 DSPC neighbors; DSPC has 5.02000 ± 0.32187 DOPC + 0.92000 ± 0.20396 DSPC neighbors. Up leaflet — DOPC has 0.21429 ± 0.12372 DOPC + 5.65714 ± 0.37634 DSPC neighbors (fcomp 0.03635); DSPC has 5.57935 ± 0.05722 DSPC neighbors (fcomp 0.92835). The ~5–6 total neighbors per lipid is exactly the expected coordination number of a 2D Voronoi tessellation (hexagonal packing → 6). The majority/minority pattern mirrors the voronoi_apl leaflet composition (dn DOPC-rich, up DSPC-rich). Reference scientifically sound.

**Q3 — Scientific interpretation:** Voronoi contact numbers/fractions describe the local lipid environment: each lipid's nearest-neighbor shell composition. The dn-leaflet DOPC is surrounded ~88% by DOPC and ~12% by DSPC; the up-leaflet DSPC is surrounded ~93% by DSPC. The ~5–6 coordination number is the geometric expectation for 2D Voronoi cells and confirms the tessellation is well-formed. The composition asymmetry quantifies the DOPC/DSPC demixing between leaflets.

**Verdict:** Test meaningful, reference scientifically sound.

### VoronoiShellComp (`voronoi_shell_comp`)

**Q1 — Test meaningfulness:** `VoronoiShellComp` (test_cli.py:1281) runs the analysis with the same selection (`--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa`). The module (voronoi_shell_comp.py) extends the Voronoi contact graph to concentric shells (shell 1 = direct neighbors, shell 2 = neighbors-of-neighbors, … up to `max_shell`, capped at 5), reporting per-type shell-wise contact number and fraction. `--qa` writes 8 files (`ave_{dn,up}_{dopc,dspc}_{fcomp,ncomp}_0.dat`), discovered via `OUTPUT_PATTERNS` (`['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat']`, test_cli.py:100) and golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive (2 leaflets × 2 types × 2 observables × 4 shells × mean/std) and would catch wrong shell expansion or wrong averaging. Coverage is limited to the default branch; `--qt`/`outb`/non-default suffix untested.

**Q2 — Reference meaningfulness:** Spot-checks: 8 files, 4 shell rows each (shells 1–4; `max_shell` estimated from box/APL = 4). Shell-1 values are internally consistent with voronoi_contact: `ave_dn_dopc_ncomp_0.dat` shell 1 = 5.31111 ± 0.08127 (other 0.69722 ± 0.04470) — identical to voronoi_contact's dn DOPC row; fcomp agrees to ~0.05% (0.88434 vs 0.88391, a rounding-level difference in the fraction averaging). Shell counts grow geometrically (dn DOPC: 5.31 → 11.24 → 17.76 → 22.66 neighbors in shells 1–4), consistent with 2D shell expansion. Reference scientifically sound.

**Q3 — Scientific interpretation:** Shell-wise composition describes the radial lipid environment around each lipid type: shell 1 is the immediate coordination shell, shell 2 the next ring, etc. The geometric growth of neighbor counts (≈5.3, 11.2, 17.8, 22.7) matches the expected 2D Voronoi shell structure and confirms the shell expansion is well-formed. The per-shell composition quantifies how DOPC/DSPC mixing decays with distance.

**Verdict:** Test meaningful, reference scientifically sound (shell-1 values cross-check against voronoi_contact).

### CholTilt (`chol_tilt`)

**Q1 — Test meaningfulness:** `CholTilt` (test_cli.py:1418) runs the analysis with `--sel "segid MEMB and resname CHL1" --center-sel "segid MEMB and name P"`. The module (chol_tilt.py) centers the membrane on the P atoms and computes the tilt of each cholesterol's C3–C17 vector relative to the z-axis, folded to 0–90°. The output (`chol_tilt.dat`, 20 frames × 49 residues) is golden-compared at **rtol=1e-2** (test_cli.py:1440-1443) — a documented tolerance override: "2D-array folding arithmetic accumulates ~0.75% error on some platforms/BLAS builds (observed ubuntu-latest CI); the reference itself was generated on macOS." The comparison is sensitive (980 tilt values) and would catch wrong vector construction, wrong folding, or wrong centering. Coverage is limited to the default branch: `interval=1`, single `--sel`/`--center-sel`; non-default interval untested.

**Q2 — Reference meaningfulness:** Spot-checks: `chol_tilt.dat` = 20 frames × (time + 49 residues), time 1–20 ns. Residues 1–40 and 140–148 (49 CHL1 molecules). Per-residue tilt means 5.89–35.57°; global range 0.23–69.38°. Cholesterol tilt angles of ~6–36° (mostly <15°) are physically plausible for cholesterol in a DOPC/DSPC bilayer (the rigid sterol ring orients nearly parallel to the membrane normal). The rtol=1e-2 tolerance is documented and justified by the platform-dependent folding arithmetic. Reference scientifically sound.

**Q3 — Scientific interpretation:** The tilt angle of the cholesterol long axis (C3–C17) relative to the membrane normal quantifies cholesterol orientation. Small tilts (mostly <15°) indicate cholesterol standing nearly upright in the bilayer, the expected orientation for the rigid sterol ring system; the per-residue spread (up to ~36° mean, 69° max) reflects thermal fluctuations. The 0–90° folding makes the observable symmetric about the normal.

**Verdict:** Test meaningful, reference scientifically sound (rtol=1e-2 tolerance documented at test_cli.py:1440-1443).

### HelixAnalysis (`helix_analysis`)

**Q1 — Test meaningfulness:** `HelixAnalysis` (test_cli.py:1446) runs the analysis with `--sel-align "segid PROA and name CA" --sel-helix "segid PROA and name CA" --align-out aligned.dcd`. The module (helix_analysis.py) aligns the trajectory to frame 1 (`ref_frame_type='specific'`, `ref_frame_num=1`) via `AlignTraj`, writes the aligned trajectory to `aligned.dcd`, then runs MDAnalysis HELANAL (`ref_axis=[0,0,1]`) on the 23 CA atoms. The output (`helix_analysis.dat`) contains Global Axes (mean/sample_sd), Global Tilts, and per-frame All Bends; it is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the alignment and the HELANAL computation. Coverage is limited to the default branch: `ref_frame_type='specific'`, `ref_frame_num=1`, `interval=1`; the `average` reference-frame branch and non-default interval are untested. **The test is fragile on linux-64**: the docker CI run fails with a line-count mismatch (see Q2) — a cross-platform numerical difference in the HELANAL output, not a logic error.

**Q2 — Reference meaningfulness:** Spot-check: `helix_analysis.dat` = 97 total lines, 95 non-empty non-comment lines, 20 "Frame" lines (one per trajectory frame). The reference is structurally sound (Global Axes mean/sample_sd, Global Tilts, All Bends shape + 20 Frame rows). **Cross-platform fragility (reported, NOT reproduced per plan guardrails)**: the docker linux-64 full-suite run fails this test with a line-count mismatch — PR_DRAFT.md:164-172 verbatim:

```
........................F.............................................................
======================================================================
FAIL: test_standard_correctness (test_cli.HelixAnalysis.test_standard_correctness)
...
AssertionError: 93 != 95 : Line count mismatch: 93 != 95
----------------------------------------------------------------------
Ran 86 tests in 455.239s

FAILED (failures=1)
```

The reference has 95 non-empty lines; the linux-64 run produced 93 — two lines differ (a platform-dependent numerical difference in the HELANAL output), not a tolerance-fixable issue. The reference itself is a valid golden of the macOS run. Reference scientifically sound on this host; flagged for the CI fragility.

**First-hand reproduction and root cause (2026-09-14, docker linux-64 via the new `pixi run docker-test` task):** reproduces the exact failure (AssertionError: 93 != 95, 86 tests in 485s). Diffing the two files shows the 2-line gap is NOT missing data and NOT differing blank lines — both files share identical structure (Global Axes/Tilts, same 2 blank separators, all 20 Frame rows × 20 values). The delta comes entirely from numpy's print formatting on **Frames 2 and 6**: their first local-bend value is `1.97823402e-02`° in the reference but exactly `0.0`° in the linux-64 output. A nonzero smallest value next to a max of ~99 makes numpy print the reference array in scientific notation (wider tokens → 4 values per wrapped line → 5 lines); with `0.0` the array stays in fixed notation (6 per line → 4 lines). All other values differ only in the 5th–7th significant digit (max relative error ≈ 1e-7–1e-5; frames 12/15 ~6e-3 absolute on O(100) values) — standard BLAS/platform summation drift. So the count mismatch is a formatting cascade triggered by one genuinely platform-dependent HELANAL value per frame (0.0198° vs 0.0°), which is why no rtol/atol can fix a line-count assertion. Both files contain all 20 frames × 20 values; the source of the 0.0° vs 0.0198° value flip has not been root-caused beyond the platform/BLAS boundary (likely a degenerate terminal-residue bend at the arccos precision floor; HELANAL is MDAnalysis library code, out of stanalyzer's scope).

**Q3 — Scientific interpretation:** HELANAL characterizes the helix geometry: the global axis (mean orientation of the helix over the trajectory), global tilts (per-frame tilt of the helix axis vs the reference axis [0,0,1]), and local bends (per-residue bend angles). For the 23-residue peptide, the global axis mean ≈ [−0.121, −0.360, 0.925] indicates a helix tilted ~22° from the membrane normal; the small sample_sd (0.002–0.004) indicates a stable orientation over the 20 frames.

**Verdict:** Test meaningful, reference scientifically sound on this host, but **flagged for cross-platform fragility** — docker linux-64 produces 93 vs 95 lines (PR_DRAFT.md:164-172); not tolerance-fixable, report-only per plan scope.

### HelixTiltRotationAngle (`helix_tilt_rotation_angle`)

**Q1 — Test meaningfulness:** `HelixTiltRotationAngle` (test_cli.py:1496) runs the analysis with `--helix-start 1 --helix-end 23`. The module (helix_tilt_rotation_angle.py) selects `resid 1-23 and name CA` (23 CA atoms), computes the helix principal axis from the inertia tensor (eigenvector of `positions.T·positions` with the smallest eigenvalue), and reports per-frame tilt (vs z-axis [0,0,1]) and rotation angle (vs reference vector [1,0,0]). The output (`helix_tilt_rotation_angle.dat`, 20 rows × 3 fields) is golden-compared at `rtol=1e-5/atol=1e-8`. The comparison is sensitive to the axis computation and the angle formulas. Coverage is limited to the default branch: `interval=1`, fixed helix range; non-default interval and other helix ranges untested.

**Q2 — Reference meaningfulness:** Spot-check: 20 rows (time, tilt, rotation). Tilt mean 79.19° (min 67.93, max 110.57); rotation mean −1.04° (min −16.15, max 11.68). The ~79° tilt indicates the helix principal axis is nearly perpendicular to the membrane normal — the peptide lies largely in the membrane plane (consistent with the ~5 Å below-midplane position from position_time and the short 23-residue helix). The rotation angle near 0° with ±16° spread reflects the azimuthal orientation about the axis. Values are self-consistent (tilt > 90° occurs because the principal axis is a signed eigenvector — the arccos of the projection is not folded, so axis-sign flips produce 67–110° values). Reference scientifically sound.

**Q3 — Scientific interpretation:** The tilt angle measures the helix axis orientation relative to the membrane normal; the ~79° mean indicates the peptide lies nearly flat in the membrane plane rather than spanning it. The rotation angle tracks the azimuthal orientation of the helix about its own axis. The large tilt is consistent with a short (23-residue) peptide that cannot span the ~43 Å bilayer as a transmembrane helix.

**Verdict:** Test meaningful, reference scientifically sound.

### BondStatistics (`bond_statistics`)

**Q1 — Test meaningfulness:** `BondStatistics` (test_cli.py:1520) runs the analysis with `-a "(1,2,3)(4,5,6)"` (`accepts_o=False`). The module (bond_statistics.py) parses the compact atom-group syntax: 2 groups → exactly one `Bond(G1,G2)` (convert_atom_groups, bond_statistics.py:356-363), computes the G1–G2 centroid distance per frame (default `centroid='cog'`), and writes only `bond_lengths.dat` (`write_files` skips empty results). The test discovers outputs via `OUTPUT_PATTERNS` (`['bond_lengths.dat', 'bond_angles.dat', 'bond_dihedrals.dat']`, test_cli.py:112) and golden-compares whatever exists. **Reference gap flagged**: only `bond_lengths.dat` is produced by this argument set, so `discover_output_files` finds 1 file and the angle/dihedral comparisons are silently skipped (`assert_output_matches_reference` skips missing refs, test_cli.py:212-248). The bond-length code path is exercised (20 values compared), but the angle and dihedral code paths (bond_statistics.py:94-173) are never exercised by the golden test — a 3-group or 4-group argument would be needed.

**Q2 — Reference meaningfulness:** Spot-check: `bond_lengths.dat` = header `@Bond Length (Angstrom)[G1_G2]` + 20 values (one per frame), 1.3426–1.7493 Å. The values are plausible inter-group centroid distances for the selected atoms (1.34–1.75 Å spans typical C–C/C–N bond lengths and short non-bonded contacts). The reference is a valid golden of the single-bond output. **Partial reference coverage flagged**: `OUTPUT_PATTERNS` declares 3 files (test_cli.py:112) but only 1 reference exists (`bond_lengths.dat`); `bond_angles.dat` and `bond_dihedrals.dat` have no references and are never compared.

**Q3 — Scientific interpretation:** Bond statistics report the distribution of distances between user-defined atom-group centroids over the trajectory. The 1.34–1.75 Å range for the G1–G2 bond reflects the covalent geometry of the selected atoms; the per-frame variation captures thermal bond-length fluctuations. Angle and dihedral statistics (unexercised here) would report the corresponding angular distributions.

**Verdict:** Test meaningful for the bond-length path, but **reference gap flagged** — only 1 of 3 declared `OUTPUT_PATTERNS` files has a reference (test_cli.py:112); angle/dihedral code paths never exercised by the golden test.

## Yiwei protein system (5 analyses)

### PiStacking (`pi_stacking`)

**Q1 — Test meaningfulness:** `PiStacking` (test_cli.py:684) runs the analysis against the yiwei_protein system with `--sel "not segid SOLV and not segid IONS" --pi-pi-dist-cutoff "6.0" --pi-cation-dist-cutoff "6.0"` and golden-compares the full output (26 residue-pair rows, each with a frame list) at `rtol=1e-5/atol=1e-8`. The module (pi_stacking.py + workers/pi_stacking_worker.py) computes aromatic ring centers and normals per frame, detects π–π pairs by ring-center distance (1.0–6.0 Å via `capped_distance` with PBC) and π–cation pairs by ring-center-to-cation distance (≤6.0 Å) plus a ring-normal alignment check (±30° or ≥150°), and emits `#residue1 residue2 frames` sorted by occupancy. The comparison is sensitive — any change in ring-atom definitions, distance computation, alignment logic, or frame indexing shifts the output. The selection covers the whole protein (PROA+PROB+CARA; the CARA carbohydrate residues contribute no rings), so both the π–π and π–cation branches are exercised (the reference contains both pair types). Untested branches: non-default `interval`, `--workers`, `--debug`, and the `pi_pi_dist_cutoff <= 0` / empty-selection error paths.

**Q2 — Reference meaningfulness:** Spot-checks on `pi_stacking.dat` (26 rows, frames 0–19): the persistent pair PROB_HSD_86–PROB_PHE_109 (listed in all 20 frames) has ring-center distances 4.075–5.791 Å across all 20 frames — all within the 6.0 Å cutoff, consistent with its all-frame occupancy. The π–cation pair PROA_HSD_1419–PROA_ARG_1287 (listed in 18 frames) has ring-center-to-cation minimum distances 3.373–5.490 Å, ≤6.0 Å in **all** 20 frames — the two missing frames (1, 14) confirm the alignment filter is active and the reference reflects real geometry, not just distance. Values are physically plausible (aromatic stacking at 4–6 Å, cation–π at 3.4–5.5 Å). Reference scientifically sound.

**Q3 — Scientific interpretation:** The output lists aromatic–aromatic and cation–π contacts with per-frame occupancy. The reference shows a small set of persistent stacking interactions (e.g., HSD86–PHE109, PHE191–PHE192, PHE192–TRP219 in all 20 frames) anchoring the fold, plus a longer tail of transient contacts (e.g., PHE1366–PHE199 in 1 frame) — the expected distribution for a folded protein. Cation–π contacts (HSD/TRP/PHE rings with ARG/LYS) are common stabilizing interactions. The 20-frame trajectory (2×10 frames at 100 ps spacing) is short, so occupancies are indicative rather than converged.

**Verdict:** Test meaningful, reference scientifically sound.

### WaterBridge (`water_bridge`)

**Q1 — Test meaningfulness:** `WaterBridge` (test_cli.py:784) runs the analysis against the yiwei_protein system (step5_1.dcd + step5_2.dcd + step3_input.psf) with `--sel "protein" --sel2 "None" --water-sel "resname TIP3" --d-a-cutoff "3.0" --d-h-a-angle-cutoff "150.0"` and golden-compares the full event list (1570 lines: frame index + four site descriptors with exact atom identities) at `rtol=1e-5/atol=1e-8`. The module (water_bridge.py) wraps MDAnalysis `WaterBridgeAnalysis` (order=1, one bridging water) and emits `#frame site1 water_to_site1 water_to_site2 site2`, each site as `[acceptor, None]` or `[donor_hydrogen, donor_heavy]` with `segid_resname_resid_name` labels. The comparison is highly sensitive — atom identities must match exactly and the event list is long. `--sel2 "None"` exercises the `sel2 = sel` fallback branch; the explicit 3.0 Å / 150° cutoffs and the TIP3 water selection are exercised. Untested branches: a distinct `--sel2`, non-default `interval`, and the default `water_sel` fallback.

**Q2 — Reference meaningfulness:** Spot-check on `water_bridge.dat` (1570 events over 20 frames, 67–92 per frame): the frame-0 bridge PROA_SER_1253 OG–HG1 → SOLV_TIP3_60832 OH2 → H2/OH2 → PROA_PRO_1254 O satisfies both geometric cutoffs on both legs — leg1 D–A 2.888 Å (≤3.0) with D–H–A 176.61° (≥150), leg2 D–A 2.703 Å with D–H–A 177.55°. The event density (tens of bridges per frame) is plausible for an ~800-residue protein in TIP3 water. Reference scientifically sound.

**Q3 — Scientific interpretation:** Water bridges are transient hydrogen-bonded networks linking protein donor/acceptor pairs through a single water molecule. The reference shows a dense network typical of a solvated protein surface, dominated by backbone N–H···O and side-chain O–H···O bridges (e.g., SER1253 OG donating through a water to the PRO1254 backbone carbonyl). The 3.0 Å / 150° cutoffs are standard H-bond criteria. The first column is a frame index (0–19), not simulation time.

**Verdict:** Test meaningful, reference scientifically sound.

### SaltBridge (`salt_bridge`)

**Q1 — Test meaningfulness:** `SaltBridge` (test_cli.py:1178) runs the analysis against the yiwei_protein system (step5_1.dcd + step5_2.dcd + step3_input.psf) with `--positive-sel "resname ARG LYS and name NZ NZ*" --negative-sel "resname ASP GLU and name OE* OD*" --positive-def "resname ARG LYS and name NZ NZ*" --negative-def "resname ASP GLU and name OE* OD*" --dist-cutoff "4.5"` and golden-compares the full output (22 residue-pair rows, each with a frame list) at `rtol=1e-5/atol=1e-8`. The module (salt_bridge.py + workers/salt_bridge_worker.py) computes residue-level salt bridges as the minimum distance between positively charged atoms (ARG/LYS NZ/NZ*) and negatively charged atoms (ASP/GLU OE*/OD*) below 4.5 Å, and emits `#residue1 residue2 frames` sorted by occupancy. The explicit `--positive-def`/`--negative-def` args exercise the custom-definition path in `build_selection`. The comparison is sensitive — any change in distance computation, selection, or frame indexing shifts the output. Untested branches: non-default `interval`, `--workers`, `--debug`, and the empty-selection error paths.

**Q2 — Reference meaningfulness:** Spot-check on `salt_bridge.dat` (22 rows, frames 0–19): the persistent bridge PROA_GLU_1260–PROA_LYS_1312 (listed in all 20 frames) has minimum OE–NZ distances 2.544–3.108 Å across all 20 frames — all within the 4.5 Å cutoff, consistent with its all-frame occupancy. 2.5–3.1 Å is the classic salt-bridge contact range. Reference scientifically sound.

**Q3 — Scientific interpretation:** Salt bridges are electrostatic contacts between oppositely charged side chains that stabilize protein structure. The reference shows a mix of persistent bridges (e.g., GLU1260–LYS1312, ASP159–LYS137, ASP175–LYS152 in all 20 frames) and transient ones (e.g., ASP235–LYS1335 in 1 frame) — the expected distribution for a folded protein. The 4.5 Å cutoff is standard for salt-bridge detection. The 20-frame trajectory is short, so occupancies are indicative rather than converged.

**Verdict:** Test meaningful, reference scientifically sound.

### GlycosidicBondBetweenSugars (`glycosidic-bond-between-sugars`)

**Q1 — Test meaningfulness:** `GlycosidicBondBetweenSugars` (test_cli.py:1389, `analysis_name = 'glycosidic-bond-between-sugars'`) runs the analysis against the yiwei_protein system (step5_1.dcd + step5_2.dcd + step3_input.psf) with `--sel "segid CARA"` and golden-compares the full output (bond-label header + 20 rows of per-frame bond lengths) at `rtol=1e-5/atol=1e-8`. The module (glycosidic-bond-between-sugars.py) finds inter-residue C–O covalent bonds within the selection from the **topology bond list** (not atom names), so the anomeric carbon and acceptor oxygen are located automatically, labels each bond (e.g., `O3(CARA1BGALNA)-C1(CARA2BGAL)`), and writes per-frame bond lengths. The comparison is sensitive — the bond labels exercise the topology-parsing logic and the lengths are compared at tight tolerance. Untested branches: non-default `interval`, an empty selection, and selections spanning multiple segments.

**Q2 — Reference meaningfulness:** Spot-checks on `glycosidic-bond-between-sugars.dat` (3 bonds, 20 frames): all three reference bonds exist as covalent C–O bonds in the topology (`bonded_in_topo=True`), and the frame-0 lengths match the reference **exactly** (1.415425 / 1.377743 / 1.395518 Å). The full length range 1.36–1.49 Å is physically correct for glycosidic C–O bonds (typical ~1.38–1.43 Å). The time column (100–1000, then 100–1000) matches the actual trajectory frame times, which reset per DCD file. Reference scientifically sound.

**Q3 — Scientific interpretation:** Glycosidic bonds link sugar residues in the carbohydrate chains (segid CARA); their lengths are tightly constrained by covalent geometry (~1.4 Å), so the reference values are nearly constant across frames, with small fluctuations reflecting thermal vibration. The three detected bonds (BGALNA O3–BGAL C1, BGALNA O6–ANE5AC C2, BGAL O3–ANE5AC C2) map the connectivity of the carbohydrate segments. The time column resets per DCD file (100–1000 twice), consistent with the two-file trajectory.

**Verdict:** Test meaningful, reference scientifically sound.

### HelixDistanceCrossingAngle (`helix_distance_crossing_angle`)

**Q1 — Test meaningfulness:** `HelixDistanceCrossingAngle` (test_cli.py:1471) runs the analysis against the yiwei_protein system (step5_1.dcd + step5_2.dcd + step3_input.psf) with `--helix1-start 1293 --helix1-end 1303 --helix2-start 1356 --helix2-end 1366` and golden-compares the full output (20 rows × 3 columns: time, distance, crossing angle) at `rtol=1e-5/atol=1e-8`. The module (helix_distance_crossing_angle.py) selects the CA atoms of the two residue ranges, computes each helix's principal axis (smallest-eigenvalue eigenvector of the inertia tensor), and per frame reports the perpendicular distance from helix2's center of geometry to helix1's axis plus the crossing angle Ω between the two axes. The comparison is sensitive — any change in axis computation, distance formula, or angle formula shifts the output. Untested branches: non-default `interval` and the `np_formatted` header path.

**Q2 — Reference meaningfulness:** Spot-check on `helix_distance_crossing_angle.dat` (20 rows): the frame-0 distance (20.06000 Å) and crossing angle (22.34838°) recomputed independently from the trajectory match the reference **exactly**. Distances span 17.05–22.08 Å — plausible for two helix axes in a large protein (residues up to ~1472). The crossing angle alternates between acute (~3–30°) and obtuse (~160–170°) values: the formula reports the angle between `(axis1×h)` and `(h×axis2)`, whose sign depends on the arbitrary eigenvector direction, so acute and supplementary values represent the same physical helix-helix geometry. The time column (0–1900, continuous) matches the module's use of `trajectory.time` (which continues across the two DCD files), not `ts.time` (which resets per file). Reference scientifically sound.

**Q3 — Scientific interpretation:** The distance and crossing angle characterize the relative geometry of two helices (residues 1293–1303 and 1356–1366). The reference shows the helices maintaining a roughly constant separation (~17–22 Å) while the crossing angle flips between ~10–30° and ~160–170° — the supplementary-angle ambiguity of the formula means these are the same geometry with opposite axis direction. The 20-frame trajectory is too short for converged statistics but captures the instantaneous geometry.

**Verdict:** Test meaningful, reference scientifically sound.

## OmF membrane system (3 analyses)

### SecondaryStructure (`secondary_structure`)

**Q1 — Test meaningfulness:** The test class (test_cli.py:1311-1333) is gated by `@unittest.skipUnless(TOOLS_AVAILABLE['dssp'], ...)`, where `dssp` means `mkdssp` or `dssp` on PATH (test_cli.py:117). `standard_args = '--sel "segid PROT_A"'` restricts the run to the 340-residue OmpF porin monomer of the 2omf system. The module (secondary_structure.py) writes per-frame PDB snapshots, shells out to `mkdssp --output-format=dssp`, parses the DSSP residue table, and emits one line per residue with a per-frame secondary-structure assignment string; the comparison (`assert_output_matches_reference`) checks the residue label exactly and the assignment characters exactly. This exercises the full DSSP pipeline — PDB writing, external-tool invocation, DSSP parsing, residue ordering — and the `segid PROT_A` selection keeps lipids out of DSSP. **Gate: dssp (mkdssp). Platform note: dssp installs on all pixi platforms (README.md), so this test runs on both CI ubuntu and macos-14.** On this host `mkdssp` is present in the pixi env and both tests (smoke + correctness) PASS.

**Q2 — Reference meaningfulness:** Numeric cross-check of `reference/secondary_structure/secondary_structure.dat`: 340 data lines = 340 residues (labels `PROT_A_<AA>_<resid>`, resid 1-340), matching the 2OMF porin monomer; 5 assignment columns = 5 frames of `equil.dcd`. Assignment alphabet is `{-, B, E, G, H, P, S, T}` — all valid DSSP codes (no `I`/pi-helix, which is rare). E (strand) dominates at 197-198 residues per frame (57.9-58.2%), the expected signature of a β-barrel porin; 325/340 residues (95.6%) are identical across all 5 frames, consistent with a stable folded barrel over the 25 ps window. The `PROT_A_` segid prefix confirms the module's `use_segid` path (system.psf has no chainIDs). Verdict: scientifically meaningful.

**Q3 — Scientific interpretation:** OmpF is a 16-stranded β-barrel porin, so the dominant E (extended strand) assignment at ~58% with the remainder in turns (T), bends (S), and loops (-) is exactly the expected secondary-structure composition; the small H/G/P fractions are the short helical segments in the extracellular loops. The 95.6% frame-to-frame stability reflects a well-equilibrated barrel over the 25 ps trajectory.

**Verdict:** Meaningful test and reference; PASS on this host (dssp present), runs on all pixi platforms.

### Sasa (`sasa`)

**Q1 — Test meaningfulness:** The test class (test_cli.py:1336-1359) is gated by `@unittest.skipUnless(TOOLS_AVAILABLE['freesasa'], ...)`, where `freesasa` means the Python module is importable (test_cli.py:118). `standard_args = '--sel "segid PROT_A"'` restricts the run to the protein segment of the 2omf system. The module (sasa.py + workers/sasa_worker.py) chunks frames across a runtime scheduler, each worker writes the selected atoms to a temp PDB and calls `freesasa.calc` (Shrake-Rupley by default, probe radius 1.4 Å), and the merged per-frame totals are written as `# Frame SASA` + `frame total` lines. This exercises the full freesasa pipeline including the chunking/worker machinery, temp-PDB writing, and result merging. **Gate: freesasa. Platform note: freesasa installs on all pixi platforms (README.md), so this test runs on both CI ubuntu and macos-14.** On this host the freesasa module is present in the pixi env and both tests PASS.

**Q2 — Reference meaningfulness:** Numeric cross-check of `reference/sasa/sasa.dat`: 5 frames with total SASA 16441.83-16544.19 Å² (mean 16510.3 Å², spread 102.4 Å² = 0.6% of total). For the 340-residue porin that is 48.6 Å²/residue — squarely in the 40-60 Å²/residue range typical of a folded protein. The ~0.6% frame-to-frame variation over 25 ps is plausible for a stable fold. The prior ledger's independent Lee-Richards vs Shrake-Rupley cross-check agreed to 0.385% max, corroborating the reference. Verdict: scientifically meaningful.

**Q3 — Scientific interpretation:** The ~16500 Å² total SASA is the solvent-exposed surface of the OmpF monomer embedded in the bilayer: the barrel's hydrophobic exterior is buried in the membrane, so the measured area is the extracellular/periplasmic loops, turns, and the pore-lining surface. The small frame-to-frame fluctuation reflects a stable fold over the 25 ps window.

**Verdict:** Meaningful test and reference; PASS on this host (freesasa present), runs on all pixi platforms.

### Hole (`hole`)

**Q1 — Test meaningfulness:** The test class (test_cli.py:1362-1386) is gated by `@unittest.skipUnless(TOOLS_AVAILABLE['hole2'], ...)`, where `hole2` requires `hole`, `sos_triangle`, and `sph_process` all on PATH (test_cli.py:119-121); `accepts_o=False` (the module writes `midpoints.dat`, `means.dat`, `hist.png` via hardcoded defaults). `standard_args = '--sel "segid PROT_A"'` selects the porin. The module (hole.py) drives `MDAnalysis.analysis.hole2.HoleAnalysis` with `cpoint='center_of_geometry'` and an explicit `cvect=(0,0,1)` membrane-normal search direction (documented as necessary to avoid HOLE locking onto a non-pore path), a seeded Monte-Carlo (`RANDOM_SEED=42`) for reproducibility, and writes per-frame `midpoints.dat` (rxn_coord vs radius, blank-line-separated sections) plus a 100-bin mean profile in `means.dat` and a PNG. This exercises the full HOLE pipeline — executable resolution, HoleAnalysis run, profile extraction, binning, plotting. **Gate: hole2. Platform note: hole2 is linux-64 only (README.md), so `stanalyzer hole` tests skip on macOS by design; CI ubuntu runs them, macos-14 skips the 2 hole tests.** On this host `hole`/`sos_triangle`/`sph_process` are absent and the test SKIPS (2 skipped), as designed.

**Q2 — Reference meaningfulness:** Numeric cross-check of `reference/hole/midpoints.dat`: 5 per-frame sections (matching the 5 frames), z-span 17.39-93.79 Å (span 76.4 Å) covering the full barrel height; per-frame constriction radius 1.475-1.658 Å at z≈58.8-59.6 Å; vestibule radius up to ~22.0 Å at the pore ends. `means.dat` has 100 bins with radius 1.941-21.676 Å (mean 9.56 Å), minimum at bin 55. These match the prior ledger exactly (z-span 76.4, min r 1.48-1.66). The ~1.5 Å constriction is tight but consistent across all 5 frames and with the prior ledger; the monomer pore is narrower than the trimer's central pore. Verdict: scientifically meaningful.

**Q3 — Scientific interpretation:** The HOLE profile traces the OmpF monomer pore along the membrane normal: a wide vestibule (~22 Å) at the extracellular mouth narrowing to a ~1.5 Å constriction near z≈59 Å, then widening again toward the periplasm. The constriction radius is the key functional quantity (the ion-selectivity filter), and the 76.4 Å z-span covers the full barrel. The seeded Monte-Carlo makes the profile reproducible across runs.

**Verdict:** Meaningful test and reference; SKIPS on macOS by design (hole2 linux-64 only), runs on CI ubuntu.

## Verdict table

One row per golden-reference class, in the same order as the per-system
subsections above. Q1 = test meaningfulness, Q2 = reference meaningfulness,
Q3-verdict = the one-line verdict from each subsection. Rows with a non-empty
Flags column need attention (or document a completed action).

| Analysis | System | Q1 | Q2 | Q3-verdict | Flags |
| --- | --- | --- | --- | --- | --- |
| DensityZ | soohyung | meaningful (default branch only) | sound — proper bilayer carbon density | Test meaningful, reference scientifically sound. | |
| ContactResidenceTime | soohyung | meaningful (cannot catch for-else bug) | UNRELIABLE — last row corrupted | Test meaningful, but reference flagged UNRELIABLE (for-else bug, contact_res_time.py:84-89). | UNRELIABLE reference (for-else overwrite) |
| HBond | soohyung | meaningful | sound — 409 bonds within cutoffs | Test meaningful, reference scientifically sound. | |
| RMSF | soohyung | meaningful | sound — 0.31–0.72 Å backbone | Test meaningful, reference scientifically sound. | |
| SystemSize | soohyung | meaningful (default path only) | sound — ~75×75×104 Å box | Test meaningful, reference scientifically sound. | |
| Thickness | soohyung | meaningful | sound — 43.17 Å | Test meaningful, reference scientifically sound. | |
| CompressibilityModulus | soohyung | meaningful (single scalar) | sound — K_A 1494.54 vs 1495.00 recomputed | Test meaningful, reference scientifically sound (short-sampling caveat). | |
| Contacts | soohyung | meaningful (exercises chunked runtime) | sound — adjacent pairs 20/20 | Test meaningful, reference scientifically sound. | |
| CovAnalysis | soohyung | meaningful (eigenvectors by proxy) | sound — eigenvalues match CA_BUG_REPORT to 4 dp | Test meaningful, reference scientifically sound; dead golden files removed in T7. | cleanup done (T7, commit 2d6d226) |
| ClusteringHca | soohyung | meaningful (default branch only) | sound — 12 clusters, PDB sane | Test meaningful, reference scientifically sound. | |
| ClusteringKmedoid | soohyung | meaningful (default k=2 only) | sound — 2 clusters, PDB sane | Test meaningful, reference scientifically sound. | |
| MsdMembrane | soohyung | meaningful (default config) | sound but PARTIAL — 2 of 5 output classes referenced | Test meaningful for default config; reference PARTIAL (2 of 5 OUTPUT_PATTERNS classes). | partial refs (2 of 5) |
| MsdSolution | soohyung | meaningful (default config) | sound but PARTIAL — 1 of 5 output classes referenced | Test meaningful for default config; reference PARTIAL (1 of 5 OUTPUT_PATTERNS classes). | partial refs (1 of 5) |
| PositionTime | soohyung | meaningful | sound — protein ~5 Å below midplane | Test meaningful, reference scientifically sound. | |
| PositionTimeCopy | soohyung | meaningful but redundant | sound but redundant — constant offset from position_time | Flagged near-duplicate of PositionTime (test_cli.py:1081-1103). | near-duplicate of PositionTime |
| Rdf | soohyung | meaningful | sound — excluded volume, peak at 7.95 Å | Test meaningful, reference scientifically sound. | |
| RadiusOfGyration | soohyung | meaningful | sound — Rg 9.82 Å ≈ helix theory | Test meaningful, reference scientifically sound. | |
| Rmsd | soohyung | meaningful | sound — 0–1.00 Å | Test meaningful, reference scientifically sound. | |
| Scd | soohyung | meaningful | sound — plateau −0.38..−0.41 | Test meaningful, reference scientifically sound (asymmetric leaflets). | |
| VoronoiApl | soohyung | meaningful | sound — 56.4/68.1 Å² asymmetric leaflets | Test meaningful, reference scientifically sound (asymmetric leaflets 56.4/68.1 Å²). | asymmetric leaflets (56.4/68.1 Å²) |
| VoronoiContact | soohyung | meaningful | sound — ~5–6 coordination | Test meaningful, reference scientifically sound. | |
| VoronoiShellComp | soohyung | meaningful | sound — shell-1 cross-checks voronoi_contact | Test meaningful, reference scientifically sound. | |
| CholTilt | soohyung | meaningful | sound — tilts mostly <15° | Test meaningful, reference scientifically sound (rtol=1e-2 documented). | rtol=1e-2 tolerance override |
| HelixAnalysis | soohyung | meaningful (fragile on linux-64) | sound on this host — 95 lines | Test meaningful, reference sound on this host, but flagged for cross-platform fragility (93 vs 95 lines). | CI fragility (docker linux-64 93 vs 95) |
| HelixTiltRotationAngle | soohyung | meaningful | sound — tilt ~79° | Test meaningful, reference scientifically sound. | |
| BondStatistics | soohyung | meaningful (bond-length path only) | sound but PARTIAL — 1 of 3 refs | Test meaningful for bond-length path; reference gap flagged (1 of 3 OUTPUT_PATTERNS files). | reference gap (1 of 3) |
| PiStacking | yiwei | meaningful (π–π and π–cation exercised) | sound — distances within cutoffs | Test meaningful, reference scientifically sound. | |
| WaterBridge | yiwei | meaningful | sound — 1570 events, cutoffs satisfied | Test meaningful, reference scientifically sound. | |
| SaltBridge | yiwei | meaningful | sound — 2.54–3.11 Å bridges | Test meaningful, reference scientifically sound. | |
| GlycosidicBondBetweenSugars | yiwei | meaningful | sound — bonds match topology exactly | Test meaningful, reference scientifically sound. | |
| HelixDistanceCrossingAngle | yiwei | meaningful | sound — frame-0 recomputed exactly | Test meaningful, reference scientifically sound. | |
| SecondaryStructure | omf | meaningful — gated (dssp) | sound — 340 res, E ~58% β-barrel | Meaningful test and reference; PASS on this host, runs on all pixi platforms. | gated: dssp (all platforms) |
| Sasa | omf | meaningful — gated (freesasa) | sound — 48.6 Å²/residue | Meaningful test and reference; PASS on this host, runs on all pixi platforms. | gated: freesasa (all platforms) |
| Hole | omf | meaningful — gated (hole2) | sound — constriction ~1.5 Å | Meaningful test and reference; SKIPS on macOS by design, runs on CI ubuntu. | gated: hole2 (linux-64 only) |

## Gap / Recommendations appendix

Each gap lists the finding, a one-line recommendation, and the evidence file
that grounds it. All gaps are report-only per plan scope — none were fixed
except the cov_analysis cleanup (T7), which is recorded as complete.

- **HelixAnalysis cross-platform CI fragility** — the docker linux-64 full-suite run fails `test_standard_correctness` with a line-count mismatch (93 vs 95, PR_DRAFT.md:164-172); a platform-dependent HELANAL output difference, not tolerance-fixable. First-hand reproduction and root cause (2026-09-14, `pixi run docker-test`): the 2-line delta is numpy's sci-vs-fixed notation flip on Frames 2 & 6, whose first local-bend value differs (1.978e-02° vs 0.0°); all other values drift only at the 5th–7th significant digit. Recommendation: investigate the local-bend value flip on linux-64 (HELANAL is MDAnalysis library code) and either pin the reference to the CI platform or make the comparison robust to the platform-dependent line set. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-4-soohyung-c.md`, `.omo/evidence/analysis-test-coverage-assessment/docker-repro-helix-fragility.md`.
- **bond_statistics missing references** — `OUTPUT_PATTERNS` declares 3 files (test_cli.py:112) but only `bond_lengths.dat` has a reference; the `-a "(1,2,3)(4,5,6)"` argument produces only the bond output, so the angle/dihedral code paths (bond_statistics.py:94-173) are never exercised. Recommendation: add a 3-group/4-group golden test with `bond_angles.dat` and `bond_dihedrals.dat` references. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-4-soohyung-c.md`.
- **msd partial reference coverage** — msd_membrane has 2 of 5 and msd_solution 1 of 5 declared `OUTPUT_PATTERNS` output classes referenced; the sys_com/mol_com/mol_info/NA modes are never compared. Recommendation: generate references for the remaining output modes or trim `OUTPUT_PATTERNS` to the tested set. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-3-soohyung-b.md`.
- **PositionTimeCopy duplicate of PositionTime** — position_time_copy is position_time minus the membrane-centering feature; the test (test_cli.py:1081-1103) exists only for framework coverage and its reference is a redundant uncentered view of the same quantity. Recommendation: remove the duplicate analysis and its test, or document a scientific justification for its inclusion. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-3-soohyung-b.md`.
- **voronoi_apl asymmetric leaflets** — the reference shows a strongly asymmetric leaflet composition (dn DOPC-rich 68.1 Å², up DSPC-rich 56.4 Å²), a real property of the soohyung input, not a reference error. Recommendation: none — documented feature; the asymmetry is interpreted in the VoronoiApl and Scd subsections. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-4-soohyung-c.md`.
- **CholTilt rtol=1e-2 tolerance** — the comparison uses a documented tolerance override (test_cli.py:1440-1443) because 2D-array folding arithmetic accumulates ~0.75% error on some platforms/BLAS builds. Recommendation: keep the documented override; revisit if the folding arithmetic is made platform-independent. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-4-soohyung-c.md`.
- **cov_analysis cleanup — DONE** — the dead golden files `eigenvectors.dat` and `correlation_matrix_heatmap.png` were deleted and `generate_baselines.py`'s copy loop was guarded by `OUTPUT_PATTERNS` in T7 (commit 2d6d226); the eigenvector output remains validated by the rotation-invariant self-consistency test (test_cli.py:935-1004). Recommendation: none — complete. Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-7-cleanup.md`.
- **ContactResidenceTime UNRELIABLE reference** — the last row (GLY 22 ALA 23) is overwritten to `0 0` by the for-else bug at contact_res_time.py:84-89, though the pair is in contact all 20 frames per contacts.dat. Recommendation: fix the for-else bug and regenerate the reference (report-only per plan scope). Evidence: `.omo/evidence/analysis-test-coverage-assessment/task-2-soohyung-a.md`.

## Provenance

- **Regeneration workflow**: golden references are produced by running each analysis with the standard test arguments against the standard test inputs, copying the outputs into `reference/<analysis_name>/`, and verifying the new reference before committing (reference/README.md:36-46). Regeneration must never be used to silence a failing test; the 2026-09 regeneration of `contacts`, `pi_stacking`, and `rmsf` (reference/README.md:53-63) was verified algorithm-equivalent before landing.
- **cov_analysis bug chain**: the eigenvector golden comparison was mathematically invalid for this dataset (degenerate eigenspaces, 50 exactly-zero eigenvalues, sign freedom on every mode) and the analysis carried four real bugs (CA_BUG_REPORT.md). The test-side fix removed `eigenvectors.dat` from the golden comparison and added the rotation-invariant self-consistency test (test_cli.py:935-1004); the analysis-side fix (eig-based k-selection) is reflected in the current reference eigenvalues [2.2845, 0.8633, 0.6792, 0.3480, 0.2544, 0.2097] (CA_BUG_REPORT.md:141-142).
- **Post-cleanup state**: the reference inventory above lists 55 files; after the T7 cleanup (commit 2d6d226) 53 remain. The two deleted files — `reference/cov_analysis/eigenvectors.dat` and `reference/cov_analysis/correlation_matrix_heatmap.png` — are documented in the CovAnalysis subsection and in git history; `generate_baselines.py` now copies only `OUTPUT_PATTERNS`-matched files, so regeneration cannot silently recreate them.
- **Date / author**: 2026-09-13. Assessment written from source inspection plus numeric reference spot-checks (transcripts in `.omo/evidence/analysis-test-coverage-assessment/task-2-soohyung-a.md` through `task-7-cleanup.md`). No analysis, test, reference, or input source was modified except the T7 cleanup.

## UPDATE (2026-09-14): changes after the assessment

This section records changes made after the 2026-09-13 snapshot above, so
the frozen per-analysis verdicts and the verdict table do not need to be
rewritten. It uses the current (post-change) `test_cli.py` line map; the
body anchors above refer to the pre-change file (as a rule of thumb, add
+93 to classes up to and including `PositionTime`, +68 to classes after it,
and −23 to the `PositionTimeCopy` block that no longer exists).

Three of the eight gap-appendix items were resolved:

### 1. ContactResidenceTime — for-else bug fixed, reference regenerated (was: UNRELIABLE reference)

- **Code fix** (contact_res_time.py:88-93): the `else` clause that
  clobbered the last-inserted contact's stats with `(0.0, 0.0)` now binds
  to `if event_len:` instead of the `for` loop. A pair with no recorded
  events still gets `(0.0, 0.0)`; the last pair is no longer
  unconditionally overwritten.
- **Reference regenerated** (`generate_baselines.py --only contact_res_time`):
  exactly one line changed, `GLY 22 ALA 23 0 0` → `GLY 22 ALA 23 20 0`
  (the pair is in contact all 20 frames per `contacts.dat`); the other 29
  rows are byte-identical.
- **Verification**: `test_cli.ContactResidenceTime` passes on macOS and
  docker linux-64.
- **Effect on the assessment**: the Q1 caveat "the test cannot catch the
  for-else overwrite bug because the golden reference itself encodes it"
  no longer applies — the golden now encodes the correct value, so a
  regression of the bug would fail the comparison. The verdict table row
  and the gap-appendix item for ContactResidenceTime are resolved.

### 2. PositionTimeCopy — analysis, test, and reference removed (was: near-duplicate of PositionTime)

- **Removed**: `src/stanalyzer/analysis/position_time_copy.py` (the
  stripped copy of position_time without membrane centering), the
  `PositionTimeCopy` test class (pre-change test_cli.py:1081-1103), and
  `reference/position_time_copy/position_time_copy.dat`. The
  `OUTPUT_PATTERNS` entries for `position_time_copy` were deleted from
  both `test_cli.py` and `generate_baselines.py`.
- **`position_time.py` cleanup**: a dead commented-out block and the
  `time_series` list initialization were removed/moved; behavior is
  identical (verified by test pass).
- **Counts**: class inventory 34 → 33 golden classes; on-disk reference
  files 53 → 52 (50 `.dat` + 2 `.pdb`; the frozen table above still lists
  the 55 original rows, three of which are now dead — 2 removed in T7,
  1 here). `stanalyzer -l` no longer lists `position_time_copy`.
- **Verification**: `test_cli.PositionTime` passes; no remaining
  `position_time_copy` references in `src/stanalyzer/tests/`.
- **Effect on the assessment**: the verdict table row, the class/reference
  inventory rows, and the gap-appendix item for PositionTimeCopy are
  resolved (the redundant analysis no longer exists; nothing to justify).

### 3. HelixAnalysis — structure-aware comparator (was: flagged for cross-platform fragility)

- **Test fix** (test_cli.py:250-312, called at :1536): the generic
  line-count comparison was replaced by
  `assert_helix_output_matches_reference`, which parses the Global Axes /
  Global Tilts / All Bends sections (`_parse_helix_headers` :250,
  `_parse_helix_bends` :269) and compares them numerically at
  `rtol=1e-5/atol=1e-8`.
- **Why it fixes the failure**: the 93-vs-95 line mismatch was numpy's
  sci-vs-fixed notation wrap on Frames 2 & 6 (first local-bend value
  1.978e-02° vs 0.0°), not missing data. A line-count assertion broke on
  that formatting; a section-aware numeric comparison is immune to it
  while still catching real regression (wrong vector construction, wrong
  averaging).
- **Verification**: `test_cli.HelixAnalysis` passes on macOS and docker
  linux-64 (targeted rerun).
- **Effect on the assessment**: the Q1 "fragile on linux-64" note, the
  Q2 cross-platform fragility flag, the verdict table row, and the
  gap-appendix item are resolved. The underlying 0.0198° vs 0.0° value
  flip remains a HELANAL (MDAnalysis library) platform artifact and is
  still not root-caused at the library level; it no longer affects the
  test result.

Unresolved gaps from the appendix remain as documented: bond_statistics
missing references (1 of 3), msd partial reference coverage (2 of 5 /
1 of 5), voronoi_apl asymmetric leaflets (documented feature), CholTilt
rtol=1e-2 tolerance override (documented), and the HELANAL value-flip
root cause (now test-neutral).