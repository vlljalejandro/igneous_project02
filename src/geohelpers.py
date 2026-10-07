"""Readers for the JEOL EPMA exports, shared by 01_olivine and 02_glass."""
import numpy as np, pandas as pd

# element (ppm) -> oxide (wt%) conversion, for detection limits
OXF = {"Si": ("SiO2", 2.1393), "Ti": ("TiO2", 1.6681), "Al": ("Al2O3", 1.8895), "Cr": ("Cr2O3", 1.4615),
       "Fe": ("FeO", 1.2865), "Mn": ("MnO", 1.2912), "Mg": ("MgO", 1.6583), "Ca": ("CaO", 1.3992),
       "Ni": ("NiO", 1.2725), "Na": ("Na2O", 1.3480), "K": ("K2O", 1.2046), "P": ("P2O5", 2.2914),
       "S": ("SO3", 2.4972), "Cl": ("Cl", 1.0), "F": ("F", 1.0)}


def read_jeol(path):
    """JEOL summary table (_oxide, _ratios): tab-separated, 3 title lines, statistics footer -> one row per analysis."""
    d = pd.read_csv(path, sep="\t", skiprows=3, skipinitialspace=True, encoding="latin-1")
    d.columns = [c.strip() for c in d.columns]
    d = d[pd.to_numeric(d["No."], errors="coerce").notna()].copy()          # drops Minimum/Maximum/Average/Sigma rows
    for c in d.columns.drop("Comment"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["No."], d["Comment"] = d["No."].astype(int), d["Comment"].str.strip()
    return d.dropna(axis=1, how="all").reset_index(drop=True)


def read_dl(path):
    """Detection limits per analysis from a JEOL _all file, converted to oxide wt% -> rows = No., columns = oxide."""
    rows, no, block = {}, None, False
    for line in open(path, encoding="latin-1"):
        if "UNK No." in line:
            no = int(line.split(":")[1].split()[0])
        elif line.lstrip().startswith("Element Peak"):
            block = True
        elif block and line.strip():
            t = line.split()
            ox, f = OXF[t[1]]
            rows.setdefault(no, {})[ox] = float(t[-1]) * f / 1e4
        else:
            block = False
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("No.")


def read_cal(path):
    """Calibration standard and assumed value per oxide ('Standard Data' block of a _std-cnd file)."""
    rows, block = [], False
    for line in open(path, encoding="latin-1"):
        t = line.split()
        if line.startswith("Standard Data"):
            block = True
        elif block and t and t[0].isdigit():
            rows.append(dict(oxide=t[1], standard=t[2], value=float(t[3])))
        elif block and rows:
            break
    return pd.DataFrame(rows)


def read_ref(path, rows):
    """Reference values from the 'major elements' sheet of Standards.xlsx.
    rows: {name in the xlsx: name used in the session}. S is converted to SO3."""
    xl = pd.read_excel(path, sheet_name="major elements", header=None).fillna("").astype(str).apply(lambda c: c.str.strip())
    is_hdr = (xl == "SiO2").any(axis=1)                       # the sheet has several header rows
    out = {}
    for name, alias in rows.items():
        r = xl.index[xl[0] == name][0]
        hdr = xl.loc[is_hdr[:r].index[is_hdr[:r]][-1]]        # nearest header row above
        s = pd.Series(xl.loc[r].values, index=hdr).drop(["Comment", "Total", ""], errors="ignore")
        out[alias] = pd.to_numeric(s, errors="coerce")
    ref = pd.DataFrame(out).T.rename_axis(columns=None)
    if "S" in ref:
        ref["SO3"] = ref.pop("S") * 2.4972
    return ref


def mean_2sd(df, by, cols):
    """n, mean and 2SD per group, columns: <oxide> and <oxide>_2SD."""
    g = df.groupby(by)[cols]
    return pd.concat([g.size().rename("n"), g.mean(), (2 * g.std()).add_suffix("_2SD")], axis=1)[
        ["n"] + [c for ox in cols for c in (ox, f"{ox}_2SD")]]


# ---------------- LA-ICP-MS (GLITTER export) ----------------

def read_glitter(path, block):
    """One block of a GLITTER csv export, transposed: rows = analyses, columns = isotopes (Li7, ...).
    block: start of the block title, e.g. 'Trace Element Concentrations, Not filtered for MDL'."""
    import io
    lines = open(path, encoding="latin-1").read().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("GLITTER!:")]
    i = next(i for i in starts if lines[i].startswith("GLITTER!: " + block))
    j = next((k for k in starts if k > i), len(lines))
    d = pd.read_csv(io.StringIO("\n".join(l for l in lines[i + 1:j] if l.strip())), index_col=0)
    return d.apply(pd.to_numeric, errors="coerce").T.rename_axis(index="label", columns=None)


def read_ref_la(path):
    """GeoReM values (µg/g) and their 95 % uncertainty from the 'Trace elements LA sorted' sheet of Standards.xlsx
    -> (values, uncertainty), rows = standard, columns = element symbol."""
    x = pd.read_excel(path, sheet_name="Trace elements LA sorted", header=None)
    cols = x.iloc[0, 1:].tolist()
    stop = cols.index("H2O")                                         # trace elements end where the oxides start
    rows = x.index[x[0].notna() & (x[0] != "uncertainty")]
    val = pd.DataFrame(x.loc[rows, 1:stop].values, index=x.loc[rows, 0], columns=cols[:stop]).apply(pd.to_numeric, errors="coerce")
    unc = pd.DataFrame(x.loc[rows + 1, 1:stop].values, index=x.loc[rows, 0], columns=cols[:stop]).apply(pd.to_numeric, errors="coerce")
    return val.rename_axis(None), unc.rename_axis(None)


# Sun & McDonough (1989), Geol. Soc. Spec. Publ. 42, 313-345: CI chondrite and primitive mantle, ppm
# (values as compiled in pyrolite, Williams et al. 2020, files CH_SunMcDonough1989 / PM_SunMcDonough1989)
CI_SM89 = {"La": 0.237, "Ce": 0.612, "Pr": 0.095, "Nd": 0.467, "Sm": 0.153, "Eu": 0.058, "Gd": 0.2055, "Tb": 0.0374,
           "Dy": 0.254, "Ho": 0.0566, "Er": 0.1655, "Tm": 0.0255, "Yb": 0.17, "Lu": 0.0254}
PM_SM89 = {"Cs": 0.032, "Rb": 0.635, "Ba": 6.989, "Th": 0.085, "U": 0.021, "Nb": 0.713, "Ta": 0.041, "K": 250,
           "La": 0.687, "Ce": 1.775, "Pb": 0.185, "Pr": 0.276, "Sr": 21.1, "P": 95, "Nd": 1.354, "Zr": 11.2, "Hf": 0.309,
           "Sm": 0.444, "Eu": 0.168, "Ti": 1300, "Gd": 0.596, "Tb": 0.108, "Dy": 0.737, "Y": 4.55, "Ho": 0.164,
           "Er": 0.48, "Tm": 0.074, "Yb": 0.493, "Lu": 0.074, "Li": 1.6}
