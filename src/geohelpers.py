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
