import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
OUT_TXT = os.path.join(HERE, "used_capital_table.txt")

V = {
    2005: dict(used="q34b_05", scale=100, estab="q6ab_05", emp="Eq1v04", assets="Eq1p04"),
    2007: dict(used="q30b_07", scale=1, estab="q6a_07", emp="Eq1x06", assets="Eq1q06"),
    2009: dict(used="q30ab_09", scale=1, estab="q6a_09", emp="EAq1x08", assets="EAq1q08"),
    2011: dict(used="q30ab_11", scale=100, estab="q6a_11", emp="EAq1x10", assets="EAq1q10"),
    2013: dict(used="q30ab_13", scale=100, estab="q6a_13", emp="EAq1x12", assets="EAq1q12"),
    2015: dict(used="q31ab_15", scale=100, estab="q6a_15", emp="EAq1k14", assets="EAq1h14"),
}

frames = []
for yr, m in V.items():
    f = os.path.join(DATA, "enterprise %d.dta" % yr)
    d = pd.read_stata(f, convert_categoricals=False,
                      columns=["id", m["used"], m["estab"], m["emp"], m["assets"]])
    d = d[~d["id"].isna()].drop_duplicates(subset="id", keep="first")
    w = pd.DataFrame({"fid": d["id"].astype(int), "wave": yr})
    w["used_share"] = (pd.to_numeric(d[m["used"]], errors="coerce")
                       / m["scale"]).clip(0, 1).values
    est = pd.to_numeric(d[m["estab"]], errors="coerce")
    w["birthyear"] = est.where((est >= 1900) & (est <= yr)).values
    emp = pd.to_numeric(d[m["emp"]], errors="coerce")
    w["size"] = pd.cut(emp.where(emp > 0), [0, 9.5, 49.5, np.inf],
                       labels=["micro", "small", "medium+"]).values
    a = pd.to_numeric(d[m["assets"]], errors="coerce")
    w["assets"] = a.where(a > 0).values
    frames.append(w)

df = pd.concat(frames, ignore_index=True)

bmode = (df.dropna(subset=["birthyear"])
           .groupby(["fid", "birthyear"]).size().rename("n").reset_index()
           .sort_values(["fid", "n", "birthyear"], ascending=[True, False, True])
           .drop_duplicates("fid").set_index("fid")["birthyear"])
age = df["wave"] - df["fid"].map(bmode)
df["age"] = age.where((age >= 0) & (age <= 80))
df["agebin"] = pd.cut(df["age"], [-0.5, 5.5, 10.5, 200],
                      labels=["0-5", "6-10", "11+"])


def firm_stats(y):
    y = np.asarray(y, float)
    n = y.size
    return 100 * y.mean(), 100 * y.std(ddof=1) / np.sqrt(n), n


def asset_stats(y, w):
    y = np.asarray(y, float)
    w = np.asarray(w, float)
    sw = w.sum()
    m = np.average(y, weights=w)
    v = np.sum((w ** 2) * (y - m) ** 2) / (sw ** 2)
    return 100 * m, 100 * np.sqrt(v), y.size


def block(frame, key, order, labmap, lead):
    b = frame.dropna(subset=["used_share", key])
    bw = b.dropna(subset=["assets"])
    head = "Firm-wtd" if lead == "firm" else "Asset-wtd"
    L = ["  %-24s %10s %8s %20s %9s" % ("", head, "s.e.", "95% CI", "N"),
         "  " + "-" * 74]
    memo = []
    for k in order:
        g = b[b[key] == k]
        gw = bw[bw[key] == k]
        if lead == "firm":
            m, se, n = firm_stats(g["used_share"])
            m2, _, _ = asset_stats(gw["used_share"], gw["assets"])
        else:
            m, se, n = asset_stats(gw["used_share"], gw["assets"])
            m2, _, _ = firm_stats(g["used_share"])
        L.append("  %-24s %10.1f %8.2f   [%6.1f, %6.1f] %9s"
                 % (labmap[k], m, se, m - 1.96 * se, m + 1.96 * se,
                    "{:,}".format(int(n))))
        memo.append("%s %.1f" % (labmap[k].split(" (")[0].strip(), m2))
    other = "asset-weighted" if lead == "firm" else "firm-weighted"
    L.append("")
    L.append("  Memo, %s: %s" % (other, "   ".join(memo)))
    return L


AGE_LAB = {"0-5": "Young (0-5 years)", "6-10": "6-10 years", "11+": "11+ years"}
SZ_LAB = {"micro": "Micro (1-9 workers)", "small": "Small (10-49 workers)",
          "medium+": "Medium+ (50+ workers)"}

O = []
O.append("=" * 78)
O.append("  USED CAPITAL IN VIETNAMESE SMEs")
O.append("  Manufacturing SMEs, Viet Nam, pooled waves 2005-2015")
O.append("=" * 78)
O.append("")
O.append("  Used capital as a share of the firm's equipment stock, percent.")
O.append("  All firms answering the question, zeros included.")
O.append("")
O.append("")
O.append("  TABLE 1.  By firm age")
O.append("")
O += block(df, "agebin", ["0-5", "6-10", "11+"], AGE_LAB, lead="firm")
O.append("")
O.append("")
O.append("  TABLE 2.  By employment size, asset-weighted")
O.append("")
O += block(df, "size", ["micro", "small", "medium+"], SZ_LAB, lead="asset")
O.append("")
O.append("")
O.append("-" * 78)
O.append("  Notes")
O.append("-" * 78)
O.append("  Sample: Survey of Small and Medium Scale MANUFACTURING Enterprises in")
O.append("  Viet Nam (CIEM, ILSSA, UCPH and UNU-WIDER), waves 2005 to 2015, pooled.")
O.append("  The used share is the respondent's own percentage, censored to the")
O.append("  unit interval. Table 1 is Figure 1 in numbers.")
O.append("")
O.append("  Firm age is the survey year minus the establishment year, with the")
O.append("  establishment year pinned at the firm's modal report across waves;")
O.append("  the self-reported year moves between waves for about a third of")
O.append("  firms. Age 10 belongs to the 6-10 group, so the top group is 11+.")
O.append("")
O.append("  Table 1 leads with the firm-weighted mean because the age gradient")
O.append("  is a statement about firms. Table 2 leads with the asset-weighted")
O.append("  mean because the model's new-to-used capital ratio is a ratio of")
O.append("  aggregate capital stocks, not an average across firms: two thirds of")
O.append("  these firms are micro but hold about a fifth of the capital. Each")
O.append("  table gives the other weighting as a memo.")
O.append("")
O.append("  Asset-weighted standard errors use the linearised variance of a")
O.append("  weighted mean. The asset distribution is concentrated -- the top")
O.append("  percentile of firms holds 29 percent of assets and the top decile 67")
O.append("  -- so a weighted mean rests on far fewer effective observations than")
O.append("  N suggests. That is why the medium+ interval is wide enough to")
O.append("  contain the micro point estimate: the asset-weighted size ordering")
O.append("  runs in the direction the model predicts, but no pairwise difference")
O.append("  between size classes is significant. Report the ordering, not an")
O.append("  estimated gradient.")
O.append("")
O.append("  Firm-weighted standard errors are unweighted and are not clustered by")
O.append("  firm; the within-firm regressions reported separately cluster.")
O.append("=" * 78)

txt = "\n".join(O)
print(txt)
with open(OUT_TXT, "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
print("\n[saved] %s" % OUT_TXT)
