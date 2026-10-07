# Igneous Geochemistry: exercise data processing (Set 4)

Processing of two electron microprobe (EPMA) sessions, basaltic glasses and olivines, and one LA-ICP-MS session of melt inclusions, from raw instrument exports to corrected, quality-checked data, figures and paper-style tables.

KAUST Igneous Geochemistry course (F. van der Zwan), data processing exercise, Set 4 (Paul & Alejandro V).

## Repository layout

```
igneous_project02/
├── data/                      raw instrument files (read only, never modified)
├── src/geohelpers.py          file readers and reference values shared by the notebooks
├── 01_olivine.ipynb           EPMA olivine, session 01/03/2017
├── 02_glass.ipynb             EPMA basaltic glasses, session 20/12/2016
├── 03_laicpms.ipynb           LA-ICP-MS melt inclusions, session 19/09/2017
├── 04_figures_tables.ipynb    figures for the data section and paper-style tables
├── outputs/                   Excel files written by the notebooks
└── figures/                   PNG figures (300 dpi): olivine/, glass/, laicpms/, report/
```

## Requirements

Python 3 with pandas, numpy, matplotlib, openpyxl and Jupyter. Tested with Python 3.12, pandas 3.0, numpy 2.4, matplotlib 3.10 and openpyxl 3.1.

```
pip install pandas numpy matplotlib openpyxl jupyter
```

## How to run

Start Jupyter in the repository root. Every path is relative to the root and set in the config cell at the top of each notebook (no command-line arguments). Run the notebooks in order:

1. `01_olivine.ipynb`
2. `02_glass.ipynb`
3. `03_laicpms.ipynb`
4. `04_figures_tables.ipynb` (reads only the outputs of 01 to 03)

Notebooks 01 to 03 are independent of each other.

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
| `170919s_9_global_tephra 1_Froukje2.csv` | GLITTER export: concentrations, 1σ errors, detection limits | 03 |
| `Standards.xlsx` | reference values: sheet *major elements* (EPMA), sheet *Trace elements LA sorted* (LA-ICP-MS) | 01, 02, 03 |
| `Information regarding your files.docx` | sample list, naming rules, LA measurement schedule | read by hand |

## Processing

### EPMA: 01 olivine and 02 glass (steps 1 to 10)

1. Join the oxide table and, for olivine, the Fo/Fa table on the analysis number.
2. Put the Comment column first, then the oxides.
3. Plot every standard element against analysis number and remove obvious outliers by eye.
4. Compute the mean, SD and correction factor CF = mean / reference. Every analysis is divided by its CF and the totals are recomputed.
   - Reference values come from `Standards.xlsx`, replaced by the more precise `_std-cnd` values for each calibration standard.
   - Glass: each oxide is corrected with its own calibration standard (VGA99; Scapolite for SO₃ and Cl; KE3 for F). MnO uses VG-2, because its calibration standard (rhodonite) was not measured in the session.
5. Flag correction factors outside 0.95 to 1.05.
6. Estimate accuracy from held-out standards: the CF is computed from part of each standard bracket and tested on the rest. For glass, VG-2 is also used, since it is not a calibration standard.
7. Olivine only: recompute Fo from MgO and FeO and compare it with the instrument's `Fo+M`.
8. Compute precision (2SD, 2RSD) and accuracy of the corrected standards.
9. Screen the samples:
   - Totals outside 98 to 102 wt% are flagged.
   - Obviously wrong analyses are excluded from the figures.
   - Values below the detection limit of their own analysis are flagged, never removed.
10. Save the final corrected data and sample or grain means (n, mean, 2SD).

### LA-ICP-MS: 03

1. Transpose the GLITTER export and assign roles from the measurement schedule in the .docx:
   - KL2 bracketing triplets: calibration.
   - KL2 between brackets, BCR-2G, GOR128-G, GOR132-G and BM90/21-G: controls.
   - NIST 612: not used.
2. Keep the concentrations not filtered for the detection limit; flag values below each analysis's own detection limit (BDL).
3. Remove obvious outliers among the standards. Average the KL2 calibration analyses per bracket and fit a line against run number.
4. Divide every analysis by CF(run) = fitted KL2(run) / KL2 reference. Run number is used as time. Ca43, the internal standard, is not corrected.
5. For every standard type, report the average, deviation from the reference, SD, SD %, 2SD and 2RSD %.
6. (Step 5b) Exclude BDL values. Flag an element when its 2RSD on KL2 exceeds 20 % or its deviation on BCR-2G exceeds 15 %.
7. (Step 6) Normalise the REE to CI chondrite and the multi-element pattern to primitive mantle (Sun & McDonough, 1989). Flag inclusions with Pb/Pb* above 5 as probable Pb contamination.

### Figures and tables: 04

- Glass Harker diagram (MgO against CaO).
- Olivine Fo against NiO, with the Fo that would be in equilibrium with the C21 glass.
- Paper-style tables: mean ± 2SD with n per sample, grain or inclusion, and the standards next to their reference values.

## Outputs

| File | Content |
|---|---|
| `outputs/olivine_step4_standards.xlsx` | SC-Olivine analyses, correction factors, all analyses raw and corrected |
| `outputs/olivine_final.xlsx` | corrected samples with flags, grain means, standard QC, hold-out accuracy, detection limits |
| `outputs/glass_step4_standards.xlsx` | standard analyses and statistics, correction factors, all analyses raw and corrected |
| `outputs/glass_final.xlsx` | corrected samples with flags, sample means, standard QC, hold-out and VG-2 accuracy, detection limits |
| `outputs/laicpms_standards_QC.xlsx` | accuracy and precision per standard type, element flags, KL2 trend |
| `outputs/laicpms_final.xlsx` | corrected analyses with flags, inclusion means, element flags, run list |
| `outputs/paper_tables.xlsx` | glass, olivine and melt inclusion tables formatted for a paper |

## Conventions

- CF = mean / reference; corrected value = raw value / CF.
- Uncertainties are 2SD (2σ) unless stated. Precision is 2RSD in %; accuracy is the deviation from the reference value in %.
- The analysis or run number is used as time; timestamps are not used.
- Flagging never deletes data: excluded analyses keep their row and a note.
- Files in `data/` are never modified.

## Changing a decision

Every judgement call is a short, commented variable in the notebook that applies it, with the reason in the markdown cell above it.

| Decision | Where |
|---|---|
| Outliers removed from the standards | `DROP` (01, 02); `DROP_RUN`, `DROP_EL` (03) |
| Sample analyses excluded from the figures | `BAD` (01, 02); `BAD_EL` (02) |
| Standard used for the glass MnO correction | `CF_STD["MnO"]` (02, step 4) |
| KL2 analyses used for calibration | `KL2_CAL` (03, config) |
| Element reportability thresholds | 03, step 5b |
| Pb contamination threshold | 03, step 6 |
| Kd and Fe²⁺/ΣFe for the olivine-melt equilibrium band | `KD`, `FE2` (04, config) |

## References

- Jochum, K. P. et al. (2005). GeoReM: a new geochemical database for reference materials and isotopic standards. *Geostandards and Geoanalytical Research* 29, 333–338. (reference values in `Standards.xlsx`)
- Roeder, P. L. & Emslie, R. F. (1970). Olivine-liquid equilibrium. *Contributions to Mineralogy and Petrology* 29, 275–289. (olivine-melt Fe-Mg Kd)
- Sun, S.-s. & McDonough, W. F. (1989). Chemical and isotopic systematics of oceanic basalts: implications for mantle composition and processes. *Geological Society, London, Special Publications* 42, 313–345. (CI chondrite and primitive mantle values)
- Williams, M. J. et al. (2020). pyrolite: Python for geochemistry. *Journal of Open Source Software* 5(50), 2314. (source of the Sun & McDonough values in `src/geohelpers.py`)
