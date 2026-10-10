# Igneous Geochemistry: exercise data processing (Set 4)

Processing of two electron microprobe (EPMA) sessions, basaltic glasses and olivines, and one LA-ICP-MS session of melt inclusions, from raw instrument exports to corrected, quality-checked data, figures and paper-style tables.

## Repository layout

```
igneous_project02/
├── data/                      raw instrument files (read only, never modified)
├── src/geohelpers.py          file readers and reference values shared by the notebooks
├── 01_olivine.ipynb           EPMA olivine, session 01/03/2017
├── 02_glass.ipynb             EPMA basaltic glasses, session 20/12/2016
├── 03_laicpms.ipynb           LA-ICP-MS melt inclusions, session 19/09/2017
├── 04_figures_tables.ipynb    report figures, host-olivine check and paper-style tables
├── report.tex                 the report (two columns, LaTeX); figures from figures/report/
├── report.pdf                 the compiled report
├── outputs/                   Excel files written by the notebooks
└── figures/                   PNG figures (300 dpi): olivine/, glass/, laicpms/, report/
```

## Requirements

Python 3 with pandas, numpy, matplotlib, openpyxl and Jupyter. Tested with Python 3.13, pandas 3.0, numpy 2.5, matplotlib 3.11 and openpyxl 3.1 (01 last run with Python 3.12, numpy 2.4 and matplotlib 3.10).

```
pip install pandas numpy matplotlib openpyxl jupyter
```

## How to run

Start Jupyter in the repository root. Every path is relative to the root and set in the config cell at the top of each notebook (no command-line arguments). Run the notebooks in order:

1. `01_olivine.ipynb`
2. `02_glass.ipynb`
3. `03_laicpms.ipynb`
4. `04_figures_tables.ipynb` (reads the outputs of 01 to 03; the KL2 example in figure 1 also reads the GLITTER file and `standard_icpms.csv`)

Notebooks 01 to 03 are independent of each other.

The report compiles from the repository root, so that the figure paths resolve: `pdflatex report.tex` (run twice for the references). The appendix tables in `report.tex` are the tables of `outputs/paper_tables.xlsx`.

## Input data

| File | Content | Used by |
|---|---|---|
| `ol_01032017_oxide.txt` | olivine oxides (wt%) per analysis | 01 |
| `ol_01032017_ratios.txt` | Fo and Fa (mol%) per analysis | 01 |
| `ol_01032017_all.txt` | full record per analysis; detection limits | 01 |
| `ol_01032017_std-cnd.txt` | analytical conditions, calibration standards and their values | 01, methods |
| `ol_01032017_mole.txt` | cation proportions | not used |
| `bgls_20122016_oxide.txt` | glass oxides (wt%) per analysis | 02 |
| `bgls_20122016_all.txt` | full record per analysis; detection limits | 02 |
| `bgls_20122016_std-cnd.txt` | analytical conditions, calibration standards and their values | 02, methods |
| `170919s_9_global_tephra 1_Froukje2.csv` | GLITTER export: concentrations, 1σ errors, detection limits | 03, 04 |
| `standard_icpms.csv` | KL2-G reference values for the LA-ICP-MS correction | 03, 04 |
| `Standards.xlsx` | reference values: sheet *major elements* (EPMA), sheet *Trace elements LA sorted* (LA-ICP-MS) | 01, 02, 03 |
| `Information regarding your files.docx` | sample list, naming rules, LA measurement schedule | read by hand |

## Processing

### EPMA: 01 olivine and 02 glass (steps 1 to 10)

1. Join the oxide table and, for olivine, the Fo/Fa table on the analysis number.
2. Put the Comment column first, then the oxides.
3. Plot every standard element against analysis number and remove obvious outliers by eye.
4. Compute the mean, SD and correction factor CF = mean / reference. Every analysis is divided by its CF and the totals are recomputed.
   - Olivine: SC-Olivine. Reference values from `Standards.xlsx`, with the more precise `_std-cnd` values for SiO₂, MgO and FeO.
   - Glass: VG-2 only, the one standard measured in the session with reference values for all 13 oxides (an arbitrary choice, reported as such). Reference values from `Standards.xlsx` (S converted to SO₃). VGA99, Scapolite and KE3 are not used.
5. Flag correction factors outside 0.95 to 1.05.
6. Estimate accuracy from held-out standards: the CF is recomputed from part of the standard analyses and tested on the rest (olivine: the third analysis of each bracket; glass: VG-2 No. 56 and 101, trained on No. 3, 4, 57 and 100).
7. Olivine only: recompute Fo from MgO and FeO and compare it with the instrument's `Fo+M`.
8. Compute precision (2SD, 2RSD) and accuracy of the corrected standards.
9. Screen the samples:
   - Totals outside 98 to 102 wt% are flagged.
   - Obviously wrong analyses, or single wrong oxides, are excluded from the averages and figures; they keep their row with a note.
   - Values below the detection limit of their own analysis are flagged, never removed.
10. Save the final corrected data and sample or grain means (n, mean, 2SD).

### LA-ICP-MS: 03 (steps 1 to 7)

1. Transpose the GLITTER export (one row per analysis). Reference values for all standards come from `Standards.xlsx`.
2. Use GLITTER's MDL-filtered concentrations to build the correction, and flag values below the detection limit (BDL). Ignore NIST 612: it does not count as time. Time is the position in the run sequence, starting with 004KL2 at time 1. KL2 is split into train (runs 5, 24, 25, 42, 44, 45, 64, 65, 66, 81, 82), test (runs 6, 26, 46, 67, 83) and discarded (run 4).
3. For every train analysis and element, CF = measured / reference (`standard_icpms.csv`); fit a line of CF against time. Ca43, the internal standard, is not corrected.
4. Correct every analysis, using the concentrations not filtered for MDL: corrected = original / (slope × time + intercept).
5. For KL2 (train, test and both) and BM90/21-G, GOR132-G, GOR128-G and BCR-2G: average, STDev, 2STDev, RSTD = 2STDev / average, and accuracy. Ca43 is listed for completeness; as the internal standard it is fixed in GLITTER, so it is left out of the medians.
6. Exclude BDL values and BCR-2G run 70 (normalised to the wrong Ca). Flag an element when its RSTD on KL2 exceeds 20 % or its accuracy on BCR-2G is worse than ±15 %.
7. Normalise the REE to CI chondrite and the multi-element pattern to primitive mantle (Sun & McDonough, 1989). Flag inclusions with Pb/Pb* above 5 as probable Pb contamination.

### Figures and tables: 04

Five report figures in `figures/report/` (300 dpi):

1. `fig1_methodology_examples.png`: a standard with its removed outlier in each EPMA session, the KL2 boron correction factor against time with its fitted line, and the EPMA totals.
2. `fig2_data_quality.png`: EPMA correction factors per oxide (±5 % band), EPMA precision against hold-out accuracy, and LA-ICP-MS reportability (KL2 precision against BCR-2G accuracy).
3. `fig3_glass_harker.png`: glass Harker diagrams against MgO (spots, and sample means ± 2SD).
4. `fig4_olivine.png`: olivine NiO, MnO and CaO against Fo, with the Fo in equilibrium with the C21 glass.
5. `fig5_melt_inclusions.png`: CI-normalised REE and primitive-mantle-normalised multi-element patterns per suite.

Host-olivine check: the LA-ICP-MS major elements of each inclusion analysis are converted to oxides and summed (67 to 416 wt%). The sum follows Si, Mg, Fe and Mn (r = 0.96 to 0.99) and partly Sc (0.69) and Li (0.39), so the laser also ablated host olivine. Elements with r > 0.5 (Si, Mg, Fe, Mn, Sc) are marked h in the inclusion table and not interpreted; all other elements have |r| ≤ 0.27.

`paper_tables.xlsx` has one sheet per dataset with n and mean ± 2SD. The EPMA standards are given as uncorrected mean, reference value, correction factor, precision and hold-out accuracy. The LA-ICP-MS standards are BCR-2G and the KL2 test split with their accuracy, plus the KL2 precision. Footnotes are generated from the outputs of 01 to 03.

## Outputs

| File | Content |
|---|---|
| `outputs/olivine_step4_standards.xlsx` | SC-Olivine analyses, correction factors, all analyses raw and corrected |
| `outputs/olivine_final.xlsx` | corrected samples with flags, grain means, standard QC, hold-out accuracy, detection limits |
| `outputs/glass_step4_standards.xlsx` | VG-2 analyses, correction factors, all analyses raw and corrected |
| `outputs/glass_final.xlsx` | corrected samples with flags, sample means, standard QC, hold-out accuracy, detection limits |
| `outputs/laicpms_standards_QC.xlsx` | statistics per standard (KL2 train, test and both; BCR-2G, GOR132-G, GOR128-G, BM90/21-G), element flags, KL2 fit |
| `outputs/laicpms_final.xlsx` | corrected analyses with flags, inclusion means, element flags, run list with time stamps and splits |
| `outputs/paper_tables.xlsx` | glass, olivine and melt inclusion tables formatted for a paper |

## Conventions

- CF = measured / reference; corrected value = raw value / CF (EPMA: one factor per oxide; LA-ICP-MS: a linear fit of CF against time).
- Uncertainties are 2SD (2σ) unless stated. Precision is 2RSD in % (called RSTD in 03); accuracy is the deviation from the reference value in %.
- Time: EPMA uses the analysis number; LA-ICP-MS uses the position in the run sequence without the NIST 612 analyses (004KL2 = 1). Instrument timestamps are not used.
- Flagging never deletes data: excluded analyses keep their row and a note.
- Files in `data/` are never modified.

## Changing a decision

Every judgement call is a short, commented variable in the notebook that applies it, with the reason in the markdown cell above it.

| Decision | Where |
|---|---|
| Outliers removed from the standards | `DROP` (01, 02); `EXCLUDE_RUN` (03, step 6) |
| Sample analyses excluded from the figures | `BAD` (01, 02); `BAD_EL` (02) |
| Standard used for the glass correction | `STD` (02, config) |
| Hold-out split for accuracy | `PER_BRACKET` (01); `TRAIN`, `TEST` (02, 03) |
| KL2 analysis discarded | `DISCARD` (03, config) |
| Element reportability thresholds | 03, step 6 |
| Pb contamination threshold | 03, step 7 |
| Kd and Fe²⁺/ΣFe for the olivine-melt equilibrium band | `KD`, `FE2` (04, config) |
| Elements marked as host olivine in the inclusions (r > 0.5 with the oxide sum) | `HOST` (04, host-olivine check) |

## References

- Jochum, K. P. et al. (2005). GeoReM: a new geochemical database for reference materials and isotopic standards. *Geostandards and Geoanalytical Research* 29, 333–338. (reference values in `Standards.xlsx`)
- Roeder, P. L. & Emslie, R. F. (1970). Olivine-liquid equilibrium. *Contributions to Mineralogy and Petrology* 29, 275–289. (olivine-melt Fe-Mg Kd)
- Sun, S.-s. & McDonough, W. F. (1989). Chemical and isotopic systematics of oceanic basalts: implications for mantle composition and processes. *Geological Society, London, Special Publications* 42, 313–345. (CI chondrite and primitive mantle values)
- Williams, M. J. et al. (2020). pyrolite: Python for geochemistry. *Journal of Open Source Software* 5(50), 2314. (source of the Sun & McDonough values in `src/geohelpers.py`)
