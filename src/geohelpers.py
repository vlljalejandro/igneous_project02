import re, numpy as np, pandas as pd
from pathlib import Path

def load_jeol_table(path):
    """_oxide / _mole / _ratios: tab-separated, 1 header row, stats footer."""
    df = pd.read_csv(path, sep="\t", skiprows=3, engine="python", quoting=3)
    df.columns = [c.strip() for c in df.columns]
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    df["No."] = df["No."].astype(str).str.strip()
    df = df[df["No."].str.fullmatch(r"\d+")].copy()      # drops Minimum/Maximum/Average/Sigma/No. of data footer
    df["No."] = df["No."].astype(int)
    df["Comment"] = df["Comment"].astype(str).str.strip()
    for c in df.columns.drop(["No.","Comment"]):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.reset_index(drop=True)

BLK = re.compile(r"Unknown Specimen No\.\s+(\d+)")
def load_jeol_all(path):
    """_all: per-analysis blocks -> (meta df, long df of net cps / SD% / DL)."""
    txt = Path(path).read_text().replace("\r","")
    blocks = re.split(r"(?=Unknown Specimen No\.)", txt)[1:]
    meta, rows = [], []
    for b in blocks:
        no  = int(BLK.search(b).group(1))
        m   = {"No.": no,
               "Comment": re.search(r"Comment\s*:\s*(.*)", b).group(1).strip(),
               "date":    re.search(r"Dated on (.*)", b).group(1).strip(),
               "kV":   float(re.search(r"Acc\. Voltage :\s*([\d.]+)", b).group(1)),
               "dia":  float(re.search(r"Probe Dia\. :\s*(\d+)", b).group(1)),
               "curr_A": float(re.search(r"Curr\.\(A\) :\s*([\dE.+-]+)", b).group(1))}
        s = re.search(r"Stage\s*:\s*X=\s*([\d.-]+)\s*Y=\s*([\d.-]+)\s*Z=\s*([\d.-]+)", b)
        m.update(zip(("X","Y","Z"), map(float, s.groups())))
        meta.append(m)
        tbl = b.split("D.L.(ppm)")[1].split("\n\n")[0]
        for ln in tbl.strip().split("\n"):
            p = ln.split()
            if len(p) >= 7 and p[0].isdigit():
                rows.append({"No.": no, "element": p[1], "peak_mm": float(p[2]),
                             "net_cps": float(p[3]), "bg_minus": float(p[4]),
                             "bg_plus": float(p[5]), "sd_pct": float(p[6].rstrip("?")),
                             "dl_ppm": float(p[-1].rstrip("?")), "flagged": "?" in ln})
    meta = pd.DataFrame(meta)
    meta["datetime"] = pd.to_datetime(meta["date"], format="%b %d %H:%M %Y")
    return meta, pd.DataFrame(rows)

def load_std_cnd(path):
    """_std-cnd: calibration standard per oxide + count times + WDS setup."""
    txt = Path(path).read_text().replace("\r","")
    std = pd.DataFrame(
        [l.split("\t")[1:4] for l in txt.split("Standard Data")[1].split("Standard Intensity")[0]
         .strip().split("\n")[1:] if l.count("\t") >= 3],
        columns=["oxide","standard","mass_pct"])
    std = std.apply(lambda s: s.str.strip())
    std["mass_pct"] = std["mass_pct"].astype(float)
    cnd = pd.DataFrame(
        [l.split("\t")[1:5] for l in txt.split("Element\tPeak\tBack\tPksk")[1]
         .split("\n \n")[0].strip().split("\n") if l.count("\t") >= 4],
        columns=["element","peak_s","bg_s","_"]).drop(columns="_")
    cnd = cnd.apply(lambda s: s.str.strip())
    cnd = cnd[cnd["peak_s"].str.match(r"[\d.]+$")].reset_index(drop=True)
    cnd[["peak_s","bg_s"]] = cnd[["peak_s","bg_s"]].astype(float)
    return std, cnd

def load_glitter(path):
    """GLITTER export -> {block name: (values df, below-MDL mask df)}; rows=analyses."""
    lines = Path(path).read_text().replace("\r","").split("\n")
    starts = [i for i,l in enumerate(lines) if l.startswith("GLITTER!:")]
    out = {}
    for i in starts:
        name = lines[i].split("GLITTER!:")[1].strip().rstrip(".")
        cols = lines[i+1].split(",")[1:]
        cols = [c for c in cols if c]
        body = []
        for l in lines[i+2:]:
            if not l.strip() or l.startswith("GLITTER!:"): break
            body.append(l.split(","))
        df = pd.DataFrame({r[0]: r[1:len(cols)+1] for r in body if r[0]}, index=cols)
        bdl = df.apply(lambda s: s.astype(str).str.strip().str.startswith("<"))
        val = df.apply(lambda s: pd.to_numeric(s.astype(str).str.lstrip("<"), errors="coerce"))
        out[name] = (val, bdl)
    return out

def parse_id(s):
    """'UnkA-C21MI_9_3' -> ('UnkA-C21MI_9', 3); last number = measurement/spot."""
    m = re.match(r"^(.*)_(\d+)$", s.strip())
    return (m.group(1), int(m.group(2))) if m else (s.strip(), np.nan)


# ---------- statistics shared by 01 / 02 ----------

Q_CRIT = {          # Dean & Dixon Q, one gap, from the lecture tables
    90: {3:.941, 4:.765, 5:.642, 6:.560, 7:.507, 8:.468, 9:.437, 10:.412,
         11:.392, 12:.376, 13:.361, 14:.349, 15:.338, 16:.329, 17:.320, 18:.313},
    95: {3:.970, 4:.829, 5:.710, 6:.625, 7:.568, 8:.526, 9:.493, 10:.466,
         11:.444, 12:.426, 13:.410, 14:.396, 15:.384, 16:.374, 17:.365, 18:.356}}

def q_test(values, conf=90):
    """Dean-Dixon Q on the min and max of a replicate set.
    Returns (index_to_reject | None, Q_exp, Q_crit). Never rejects more than one point."""
    v = pd.Series(values).dropna().sort_values()
    n = len(v)
    if n < 3 or n > max(Q_CRIT[conf]) or v.iloc[-1] == v.iloc[0]:
        return None, np.nan, np.nan
    rng = v.iloc[-1] - v.iloc[0]
    q_lo = (v.iloc[1] - v.iloc[0]) / rng
    q_hi = (v.iloc[-1] - v.iloc[-2]) / rng
    q_exp, idx = (q_hi, v.index[-1]) if q_hi >= q_lo else (q_lo, v.index[0])
    qc = Q_CRIT[conf][n]
    return (idx if q_exp > qc else None), q_exp, qc

def summarize(df, cols, by):
    """n, mean, 2SD, 2RSD%, 2SE per group per oxide, long format."""
    g = df.groupby(by)[cols]
    n, mu, sd = g.count(), g.mean(), g.std()
    out = pd.concat({"n": n, "mean": mu, "2SD": 2*sd, "2RSD%": 200*sd/mu,
                     "2SE": 2*sd/np.sqrt(n)}, axis=1).stack(future_stack=True)
    out.index.names = list(np.atleast_1d(by)) + ["oxide"]
    return out.reset_index()


def load_ref_majors(path, sheet="major elements"):
    """Standards.xlsx 'major elements' -> wt% oxide table indexed by standard name.
    The sheet reports sulphur as element S; converted to SO3 (x 2.4972)."""
    r = pd.read_excel(path, sheet)
    r.columns = [str(c).strip() for c in r.columns]
    r = r.set_index(r.columns[0])
    r.index = r.index.astype(str).str.strip()
    r = r[r.index.str.startswith("Std")].apply(pd.to_numeric, errors="coerce")
    if "S" in r.columns:
        r["SO3"] = r["S"] * 2.4972
    r.index = r.index.str.replace("^Std ", "", regex=True)
    return r


# ---------- olivine-specific (02) ----------

OXIDE_MW = {"SiO2": 60.084, "TiO2": 79.866, "Al2O3": 101.961, "FeO": 71.844, "MnO": 70.937,
            "MgO": 40.304, "CaO": 56.077, "NiO": 74.692, "Cr2O3": 151.990, "Na2O": 61.979,
            "K2O": 94.196, "P2O5": 141.945}

def load_ref_olivine(path, sheet="major elements"):
    """The SC_Olivine reference row of Standards.xlsx, which sits in a second table
    further down the sheet with its own header row. Returns a wt% oxide Series."""
    raw = pd.read_excel(path, sheet, header=None)
    i = raw.index[raw[0].astype(str).str.strip().str.lower() == "sc_olivine"][0]
    hdr = raw.loc[i - 1].astype(str).str.strip()
    s = pd.Series(raw.loc[i].values, index=hdr).apply(pd.to_numeric, errors="coerce")
    return s[[c for c in s.index if c in OXIDE_MW]].dropna()

def cation_ratio(df, si="SiO2", m=("MgO", "FeO", "MnO", "CaO", "NiO")):
    """Divalent cations per Si. Must be 2.000 for stoichiometric olivine, whatever the Fo.
    Independent of any uniform scale error on the analysis, since it is a ratio."""
    M = sum(df[o] / OXIDE_MW[o] for o in m if o in df)
    return M / (df[si] / OXIDE_MW["SiO2"])

def forsterite(df):
    """Fo = 100 Mg/(Mg+Fe), mol%. For olivine this equals Mg# (all Fe as Fe2+)."""
    mg, fe = df.MgO / OXIDE_MW["MgO"], df.FeO / OXIDE_MW["FeO"]
    return 100 * mg / (mg + fe)
