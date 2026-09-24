import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
OUT_TXT = os.path.join(HERE, "by_wave_table.txt")
OUT_TEX = os.path.join(HERE, "by_wave_table.tex")

V = {
    2005: dict(new="q34a_05",  used="q34b_05",  self="q34c_05",  mach="q32_05",
               wt="Eq1s04",    scale=100, mkt=False, item="q34"),
    2007: dict(new="q30a_07",  used="q30b_07",  self="q30c_07",  mach="q28_07",
               wt="Eq1t06",    scale=1,   mkt=False, item="q30"),
    2009: dict(new="q30aa_09", used="q30ab_09", self="q30ac_09", mach="q28_09",
               wt="EAq1t08",   scale=1,   mkt=True,  item="q30a"),
    2011: dict(new="q30aa_11", used="q30ab_11", self="q30ac_11", mach="q28_11",
               wt="EAq1t10",   scale=100, mkt=True,  item="q30a"),
    2013: dict(new="q30aa_13", used="q30ab_13", self="q30ac_13", mach="q28_13",
               wt="EAq1t12",   scale=100, mkt=True,  item="q30a"),
    2015: dict(new="q31aa_15", used="q31ab_15", self="q31ac_15", mach="q29_15",
               wt="q65ac_15",  scale=100, mkt=True,  item="q31a"),
}

rows = []
pool = []
for yr, m in sorted(V.items()):
    path = os.path.join(DATA, "enterprise %d.dta" % yr)
    want = ["id", m["new"], m["used"], m["self"], m["mach"], m["wt"]]
    d = pd.read_stata(path, convert_categoricals=False, columns=want)
    d = d[~d["id"].isna()].drop_duplicates(subset="id", keep="first")
    firms = len(d)

    num = lambda c: pd.to_numeric(d[c], errors="coerce")
    n, u, s = num(m["new"]), num(m["used"]), num(m["self"])
    mach = num(m["mach"])
    sc = float(m["scale"])

    raw_reporting = int(u.notna().sum())
    user = mach.between(2, 5)
    tot = pd.concat([n, u, s], axis=1).sum(axis=1, min_count=3)
    complete = n.notna() & u.notna() & s.notna() & tot.between(0.99 * sc, 1.01 * sc)
    keep = user & complete
    reporting = int(keep.sum())

    us = (u[keep] / sc).clip(0, 1)
    any_used = 100.0 * (us > 0).mean()
    users_only = 100.0 * us[us > 0].mean()
    all_firms = 100.0 * us.mean()

    w = pd.to_numeric(d[m["wt"]], errors="coerce")
    g = pd.DataFrame({"n": n, "u": u, "w": w})[keep]
    g = g[(g["w"] > 0) & g["w"].notna() & (g["n"] + g["u"] > 0)]
    tn = float((g["n"] * g["w"]).sum())
    tu = float((g["u"] * g["w"]).sum())
    agg = 100.0 * tu / (tn + tu) if (tn + tu) > 0 else np.nan
    knku = tn / tu if tu > 0 else np.nan

    rows.append(dict(year=yr, item=m["item"], firms=firms, raw=raw_reporting,
                     reporting=reporting, any_used=any_used,
                     users_only=users_only, all_firms=all_firms,
                     agg=agg, knku=knku, nw=len(g), mkt=m["mkt"]))
    pool.append(us)

allrep = pd.concat(pool, ignore_index=True)
rows.append(dict(year="Pooled", item="", firms=sum(r["firms"] for r in rows),
                 raw=sum(r["raw"] for r in rows), reporting=len(allrep),
                 any_used=100.0 * (allrep > 0).mean(),
                 users_only=100.0 * allrep[allrep > 0].mean(),
                 all_firms=100.0 * allrep.mean(),
                 agg=np.nan, knku=np.nan, nw=np.nan, mkt=True))


def fmt(x, nd=1):
    return "---" if x is None or (isinstance(x, float) and np.isnan(x)) else "%.*f" % (nd, x)


def grp(n):
    return "---" if (isinstance(n, float) and np.isnan(n)) else "{:,}".format(int(n))


T = ["=" * 92,
     "  THE USED-CAPITAL ITEM, BY SURVEY ROUND",
     "  Survey of Small and Medium Scale Manufacturing Enterprises, Viet Nam",
     "=" * 92, "",
     "  %-7s %-6s %8s %9s %9s %9s %10s %9s %10s %8s"
     % ("round", "item", "firms", "raw item", "in block", "any used",
        "users only", "all firms", "aggregate", "Kn/Ku"),
     "  " + "-" * 90]
for r in rows:
    T.append("  %-7s %-6s %8s %9s %9s %9s %10s %9s %10s %8s"
             % (r["year"], r["item"], grp(r["firms"]), grp(r["raw"]),
                grp(r["reporting"]), fmt(r["any_used"]), fmt(r["users_only"]),
                fmt(r["all_firms"]), fmt(r["agg"], 2) + ("" if r["mkt"] else "*"),
                fmt(r["knku"], 4)))
T += ["",
      "  firms      = firms interviewed in the round",
      "  raw item   = firms with a non-missing used-capital item in the raw file",
      "  in block   = the universe used here: machinery user with a complete",
      "               new/used/self triple.  2005 loses 295 non-machinery firms",
      "               that the raw file nonetheless gives a vintage triple.",
      "  any used   = % of in-block firms with a strictly positive used share",
      "  users only = mean used share among those firms, %",
      "  all firms  = mean used share across in-block firms, %",
      "  aggregate  = sum(w u) / sum(w (n+u)), self-constructed dropped, w the",
      "               firm's machinery and equipment stock.  Kn/Ku is its",
      "               complement; 2013 is the KnKu_ratio calibration target.",
      "  *          = weight at accounting value, not market price (2005, 2007);",
      "               not comparable in level with 2009-2015.",
      "=" * 92]
txt = "\n".join(T)

L = ["% rows for tab:by_wave -- generated by by_wave_table.py",
     "% columns: round & interviewed & in block & any used & users only & aggregate"]
for r in rows:
    if r["year"] == "Pooled":
        L.append("\\midrule")
        L.append("Pooled & $%s$ & $%s$ & $%s$ & $%s$ & --- \\\\"
                 % (grp(r["firms"]), grp(r["reporting"]),
                    fmt(r["any_used"]), fmt(r["users_only"])))
    else:
        L.append("%d & $%s$ & $%s$ & $%s$ & $%s$ & $%s$%s \\\\"
                 % (r["year"], grp(r["firms"]), grp(r["reporting"]),
                    fmt(r["any_used"]), fmt(r["users_only"]), fmt(r["agg"]),
                    "" if r["mkt"] else "$^{*}$"))
tex = "\n".join(L)

print(txt)
print()
print(tex)
with open(OUT_TXT, "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
with open(OUT_TEX, "w", encoding="utf-8") as fh:
    fh.write(tex + "\n")
print("\n[saved] %s\n[saved] %s" % (OUT_TXT, OUT_TEX))
