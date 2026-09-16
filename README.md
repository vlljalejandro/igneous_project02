# Project: Igneous Geochemistry — Exercise Data Processing (EPMA + LA-ICP-MS)

## Context
KAUST Igneous Geochemistry course (F. van der Zwan). **Set 4 — Paul & Alejandro V**
(paired work; agree early who drafts which report section and keep one shared notebook set).
Deliverable: the "Exercise data processing written reports — 20%" report, due Thu 10/08.
Two EPMA datasets and one LA-ICP-MS dataset, raw instrument output → publication-ready
tables and figures + a 5½-page report.

## Data — everything in `data/`
```
data/
  Information regarding your files.docx  <- READ FIRST: sample list, which analyses are
                                            standards vs unknowns, which set is Set 4
  Standards.xlsx                         <- certified/reference values (GeoReM-style)
  bgls_20122016_all.txt                  EPMA #1 — basaltic glasses, session 20/12/2016
  bgls_20122016_oxide.txt                  _all   = counts, net cps, SD%, detection limits
  bgls_20122016_std-cnd.txt                _oxide = oxide wt% + totals
  ol_01032017_all.txt                    EPMA #2 — olivines, session 01/03/2017
  ol_01032017_oxide.txt                    _mole  = cation/mole proportions
  ol_01032017_mole.txt                     _ratios= e.g. Fo, Mg#
  ol_01032017_ratios.txt                   _std-cnd = standards + analytical conditions
  ol_01032017_std-cnd.txt                            (kV, nA, beam dia, count times)
  170919s_9_global_tephra 1_Froukje2.csv LA-ICP-MS — GLITTER-style export, 19/09/2017
```
File-content assumptions above are inferred from filenames — verify against the .docx in
notebook 00 and correct this file if wrong. Confirm in the .docx that these three files
are in fact Set 4's; if the sets are split differently, fix the list before processing.
`*_std-cnd.txt` is the primary source for the report's methodology section; do not
paraphrase conditions from memory.

## Rules
- Never modify files in `data/`. Outputs go to `outputs/` and `figures/`.
- Hardcoded paths in a config cell at the top of each notebook. No argparse, no CLI.
- Prefer short, direct pandas/numpy. Factor into `src/geohelpers.py` only if used by ≥2
  notebooks (candidates: standards loader, Q-test, drift fit, 2SD/2SE).
- Every filtering or correction decision gets a one-line justification in a markdown cell —
  the methodology section is built from these, and "WHY" is explicitly graded.
- Figures: matplotlib, one plot per figure, PNG 300 dpi + kept inline. Error bars are 2SE
  unless stated. Units (wt% vs ppm) on every axis.
- Data below detection limit: flag, never silently drop or set to zero.

## Notebook sequence
**00_inspect.ipynb** — Open every file raw (`head`, no parser). Determine delimiter, header
depth, comment/sample-ID column, encoding. JEOL exports often have multi-line or fixed-width
headers; the filename with a space needs quoting/`Path`. Write one loader per file type,
confirm row counts against the .docx, save tidy dataframes to `outputs/interim/`. Nothing
else happens here.

**01_epma_glass.ipynb** and **02_epma_olivine.ipynb** — same skeleton:
1. Split into calibration standards, control standards run as unknowns (e.g. VG2, VGA99,
   San Carlos Ol, KL2), and samples, using the Comment/Sample-ID field.
2. Precision: per standard per oxide → mean, SD, 2SD, 2RSD%, n.
3. Accuracy: deviation of mean from certified value in `Standards.xlsx` (absolute and %).
   Flag any oxide where |accuracy error| > precision — the problem criterion from the lecture.
4. Drift: plot each standard oxide vs analysis order. If a trend exists, derive a correction
   factor (certified / measured mean) for the session or for batches between standard
   brackets, apply to samples, record the factors.
5. Filter samples: totals outside 98–102% (state whether volatiles or oxidation state could
   legitimately explain low totals); element-specific contamination (olivine analyses with
   anomalous Al₂O₃/CaO/K₂O = beam overlap onto glass or inclusions; glass analyses with Ol
   contamination). Report how many points removed and why.
6. Outliers within a sample's replicate spots: Q-test (Q_exp = gap/range, lecture tables)
   before exclusion. Never drop a point just because it looks odd.
7. Averages: minimum 3 points; report n, mean, 2SD, 2SE per sample. Check for two
   populations (different glass chips, zoned olivine) rather than averaging across them.
8. Compare low-concentration elements (Cl, P, Ni, Cr) to detection limits in `_all`.
9. Olivine only: recompute Fo and Mg# from the oxides rather than trusting `_ratios`, and
   confirm they agree.

**03_laicpms.ipynb**
1. Identify calibration standard (likely NIST SRM 610/612), secondary/control standards
   (KL2, ATHO, BCR2G, GOR128-G, GOR132-G, BM90 — confirm from the .docx), and unknowns.
2. Confirm which element is the internal standard (constant within a material, e.g. Ca43 or
   Si29) and which EPMA values it was normalized to.
3. Per standard: mean, SD, %RSD, 2RSD, n → precision; deviation from reference → accuracy.
   Produce the accuracy/precision table for every reported element.
4. Drift: k = reference/measured per standard analysis, plot vs analysis number, fit a line,
   decide session-wide vs batch-wise correction, apply.
5. Oxide/molecular interference where relevant — ⁴⁵Sc on olivine is partly ²⁹Si¹⁶O: get the
   Si* oxide-production factor from a Si-only material (quartz) and correct. Check REE for
   LREE-oxide overlaps.
6. Decide element by element what is reportable. 2RSD > ~20% or accuracy worse than ~15% →
   exclude or report with an explicit caveat. Justify in writing.
7. Optionally normalize majors to 100% and compare against EPMA for the same samples as an
   independent consistency check.

**04_figures_and_tables.ipynb** — final figures and the paper-style table.

## Figures (pick the subset the report can carry)
- Standard time series with certified value ± 2SD band (drift evidence)
- Accuracy vs precision per element, per technique
- Histogram of analytical totals with the 98–102% window marked
- Harker diagrams (MgO vs major oxides) for the glasses
- Fo vs NiO / CaO / MnO for the olivines
- Primitive-mantle-normalized multi-element and chondrite-normalized REE patterns
- EPMA vs LA-ICP-MS 1:1 plot for shared elements

## Report (5½ pages)
| Section | Pages | Content |
|---|---|---|
| Introduction | ½ | What the samples are, what the analytical question is |
| Detailed methodology | 2 | Every processing step **and why**, with example figures |
| Assessment of data quality | 1 | Precision, accuracy, detection limits, what is and isn't reportable |
| Methodology for manuscript | ½ | The same work as a journal methods paragraph: instruments, conditions, standards, their reproducibility |
| Short presentation of data | 1 | The geochemistry itself, 2–3 figures |
| Conclusion | ½ | What was learned, interesting features of the data |

Optional appendix: Excel table formatted as for a paper — sample, n, mean ± 2SE, with
standards and reference values in the same table.

## Working style
Concise and direct. Minimal, targeted code over architecture. Ask before inventing data
semantics — if a column's meaning is ambiguous, check the .docx or flag it rather than guess.
