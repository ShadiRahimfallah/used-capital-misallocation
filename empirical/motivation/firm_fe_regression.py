import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
OUT_TXT = os.path.join(HERE, "firm_fe_results.txt")
OUT_TEX = os.path.join(HERE, "firm_fe_table.tex")


class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, s):
        self.f.write(s)
        self.f.flush()
        self.stdout.write(s)

    def flush(self):
        self.f.flush()
        self.stdout.flush()

    def close(self):
        self.f.close()


_tee = Tee(OUT_TXT)
sys.stdout = _tee

print("=" * 78)
print("FIRM FIXED-EFFECTS REGRESSION - used capital, Vietnam SME panel 2005-2015")
print("=" * 78)

V = {
    2005: dict(used="q34b_05", scale=100, estab="q6ab_05", emp="Eq1v04",
               debt="Eq1u04", assets="Eq1p04", imp=None, sector="q13_05",
               prov="q3be_05"),
    2007: dict(used="q30b_07", scale=1, estab="q6a_07", emp="Eq1x06",
               debt="Eq1v06", assets="Eq1q06", imp=None, sector="q17a_07",
               prov="province"),
    2009: dict(used="q30ab_09", scale=1, estab="q6a_09", emp="EAq1x08",
               debt="EAq1v08", assets="EAq1q08", imp="EAq1y08", sector="q17a_09",
               prov="q3be_09"),
    2011: dict(used="q30ab_11", scale=100, estab="q6a_11", emp="EAq1x10",
               debt="EAq1v10", assets="EAq1q10", imp="EAq1y10", sector="q17a_11",
               prov="q3be_11"),
    2013: dict(used="q30ab_13", scale=100, estab="q6a_13", emp="EAq1x12",
               debt="EAq1v12", assets="EAq1q12", imp="EAq1y12", sector="q17a_13",
               prov="q3be_13"),
    2015: dict(used="q31ab_15", scale=100, estab="q6a_15", emp="EAq1k14",
               debt="EAq1i14", assets="EAq1h14", imp=None, sector="q17a_15",
               prov="q3ce1_15"),
}

frames = []
print("\n--- 1. BUILD PANEL ---------------------------------------------------")
for yr, m in V.items():
    f = os.path.join(DATA, "enterprise %d.dta" % yr)
    cols = ["id"] + [v for k, v in m.items() if k != "scale" and v is not None]
    d = pd.read_stata(f, convert_categoricals=False, columns=list(dict.fromkeys(cols)))
    d = d[~d["id"].isna()].drop_duplicates(subset="id", keep="first")
    w = pd.DataFrame({"fid": d["id"].astype(int), "wave": yr})

    w["used_share"] = (pd.to_numeric(d[m["used"]], errors="coerce") / m["scale"]).clip(0, 1).values

    est = pd.to_numeric(d[m["estab"]], errors="coerce")
    est = est.where((est >= 1900) & (est <= yr))
    w["birthyear"] = est.values
    age = yr - est
    w["age"] = age.where((age >= 0) & (age <= 80)).values
    w["young"] = (w["age"] <= 5).astype(float).where(w["age"].notna())

    emp = pd.to_numeric(d[m["emp"]], errors="coerce")
    emp = emp.where(emp > 0)
    w["emp"] = emp.values
    w["size"] = pd.cut(emp, [0, 9.5, 49.5, np.inf],
                       labels=["micro", "small", "medium+"]).values

    debt = pd.to_numeric(d[m["debt"]], errors="coerce")
    debt = debt.where(debt >= 0)
    w["zerodebt"] = (debt == 0).astype(float).where(debt.notna()).values

    a = pd.to_numeric(d[m["assets"]], errors="coerce")
    w["lnA"] = np.log(a.where(a > 0)).values

    if m["imp"]:
        w["importer"] = (pd.to_numeric(d[m["imp"]], errors="coerce") > 0).astype(float).values
    else:
        w["importer"] = np.nan

    w["sector"] = pd.to_numeric(d[m["sector"]], errors="coerce").values
    w["prov"] = pd.to_numeric(d[m["prov"]], errors="coerce").values
    frames.append(w)

df = pd.concat(frames, ignore_index=True)

df["age_raw"] = df["age"]
df["young_raw"] = df["young"]
bmode = (df.dropna(subset=["birthyear"])
           .groupby(["fid", "birthyear"]).size().rename("n").reset_index()
           .sort_values(["fid", "n", "birthyear"], ascending=[True, False, True])
           .drop_duplicates("fid").set_index("fid")["birthyear"])
df["birthyear_fix"] = df["fid"].map(bmode)
age = df["wave"] - df["birthyear_fix"]
df["age"] = age.where((age >= 0) & (age <= 80))
df["young"] = (df["age"] <= 5).astype(float).where(df["age"].notna())

df["used_share"] = 100.0 * df["used_share"]

df["small"] = (df["size"] == "small").astype(float).where(df["size"].notna())
df["medplus"] = (df["size"] == "medium+").astype(float).where(df["size"].notna())
df["zd_young"] = df["zerodebt"] * df["young"]
for lo, hi, nm in [(-0.5, 2.5, "a02"), (2.5, 5.5, "a35"), (5.5, 10.5, "a610")]:
    df[nm] = ((df["age"] > lo) & (df["age"] <= hi)).astype(float).where(df["age"].notna())

print("  pooled firm-wave observations: %d   unique firm ids: %d"
      % (len(df), df.fid.nunique()))
nw = df.groupby("fid").size().value_counts().sort_index()
print("  firms by number of waves observed:")
for k, v in nw.items():
    print("     %d wave(s): %5d firms  (%6d obs)" % (k, v, k * v))

print("\n--- 2. PANEL-ID VALIDATION -------------------------------------------")
chk = df.dropna(subset=["birthyear"]).groupby("fid")["birthyear"].nunique()
multi = df[df.groupby("fid")["fid"].transform("size") >= 2].dropna(subset=["birthyear"])
chk2 = multi.groupby("fid")["birthyear"].nunique()
spread = (multi.groupby("fid")["birthyear"].max()
          - multi.groupby("fid")["birthyear"].min())
print("  firms whose reported birth year is constant across all their waves: "
      "%d of %d (%.1f%%)" % ((chk == 1).sum(), len(chk), 100 * (chk == 1).mean()))
print("  among multi-wave firms only:                                        "
      "%d of %d (%.1f%%) constant" % ((chk2 == 1).sum(), len(chk2),
                                      100 * (chk2 == 1).mean()))
print("  disagreement in years (max-min) among the inconsistent firms: "
      "median %.0f, p75 %.0f, p90 %.0f, max %.0f"
      % tuple(spread[spread > 0].quantile([.5, .75, .9, 1.0])))
print("  -> birth year pinned at the firm's modal report; Age and Young rebuilt "
      "from it.\n     Self-reported version retained as age_raw/young_raw "
      "(robustness column).")
pv = df.dropna(subset=["prov"]).groupby("fid")["prov"].nunique()
print("  firms whose province is constant across all their waves:            "
      "%d of %d (%.1f%%)" % ((pv == 1).sum(), len(pv), 100 * (pv == 1).mean()))
sc = df.dropna(subset=["sector"]).groupby("fid")["sector"].nunique()
print("  firms whose sector   is constant across all their waves:            "
      "%d of %d (%.1f%%)" % ((sc == 1).sum(), len(sc), 100 * (sc == 1).mean()))
print("  -> province/sector are (nearly) firm-invariant, hence absorbed by alpha_i.")

print("\n--- 3. IDENTIFYING VARIATION FOR THE WITHIN ESTIMATOR ----------------")
needF = ["used_share", "zerodebt", "young", "lnA", "small", "medplus", "wave", "fid"]
dF = df.dropna(subset=needF).copy()
dF = dF[dF.groupby("fid")["fid"].transform("size") >= 2].copy()
print("  estimation sample (non-missing, firm seen >=2 waves): %d obs, %d firms"
      % (len(dF), dF.fid.nunique()))


def switchers(d, col):
    g = d.groupby("fid")[col]
    v = g.nunique()
    return int((v > 1).sum()), int(v.size)


for c, lab in [("young", "Young (age<=5) changes"),
               ("zerodebt", "ZeroDebt changes"),
               ("zd_young", "ZeroDebt x Young changes"),
               ("small", "size class 'small' changes")]:
    s, t = switchers(dF, c)
    print("  firms for which %-28s : %5d of %5d (%5.1f%%)  <- these firms carry the coefficient"
          % (lab, s, t, 100 * s / t))

wsd = dF.groupby("fid")["used_share"].std()
print("  within-firm sd of the used share: mean %.1f pp, median %.1f pp; "
      "%d firms (%.1f%%) never change it"
      % (wsd.mean(), wsd.median(), (wsd == 0).sum(), 100 * (wsd == 0).mean()))
print("  transitions of Young among its switchers: %d firm-pairs go 1->0 "
      "(cross age 5), %d go 0->1 (should be ~0)"
      % (((dF.sort_values(["fid", "wave"]).groupby("fid")["young"].diff()) == -1).sum(),
         ((dF.sort_values(["fid", "wave"]).groupby("fid")["young"].diff()) == 1).sum()))

texnum = lambda v: format(int(v), ",").replace(",", "{,}")
STAR = lambda p: "***" if p < .01 else ("**" if p < .05 else ("*" if p < .1 else ""))
RES = {}


def pack(names, b, se, p, n, r2, G, note, r2lab="R2"):
    o = {k: (b[i], se[i], p[i]) for i, k in enumerate(names)}
    o.update(_N=int(n), _R2=r2, _G=G, _note=note, _R2lab=r2lab)
    return o


def pooled(d, formula, keep, note):
    m = smf.ols(formula, d).fit(cov_type="cluster", cov_kwds={"groups": d["fid"]})
    keep = [k for k in keep if k in m.params.index]
    return m, pack(keep, [m.params[k] for k in keep], [m.bse[k] for k in keep],
                   [m.pvalues[k] for k in keep], m.nobs, m.rsquared,
                   d.fid.nunique(), note)


def firm_fe(d, xcols, note, ycol="used_share"):
    d = d.copy()
    yrd = pd.get_dummies(d["wave"], prefix="yr", drop_first=True).astype(float)
    X = pd.concat([d[xcols].astype(float).reset_index(drop=True),
                   yrd.reset_index(drop=True)], axis=1)
    y = d[ycol].astype(float).reset_index(drop=True)
    g = d["fid"].reset_index(drop=True).to_numpy()

    Xd = X - X.groupby(g).transform("mean")
    yd = y - y.groupby(g).transform("mean")
    keepc = [c for c in Xd.columns if Xd[c].abs().max() > 1e-10]
    Xd = Xd[keepc]

    m = sm.OLS(yd, Xd).fit(cov_type="cluster", cov_kwds={"groups": g})
    n, k, G = int(m.nobs), Xd.shape[1], len(np.unique(g))
    nabs = G
    infl = np.sqrt((n - k) / float(n - k - nabs))
    bse = m.bse * infl
    dfree = G - 1
    pval = 2 * stats.t.sf(np.abs(m.params / bse), dfree)
    within_r2 = 1 - m.ssr / float(np.sum(yd ** 2))
    names = [c for c in xcols if c in Xd.columns]
    idx = [list(Xd.columns).index(c) for c in names]
    note += ("; SE clustered by firm and d.o.f.-corrected for the %d absorbed "
             "firm effects (SEs are %.0f%% larger than the uncorrected ones)"
             % (nabs, 100 * (infl - 1)))
    out = pack(names, [m.params.iloc[i] for i in idx], [bse.iloc[i] for i in idx],
               [pval[i] for i in idx], n, within_r2, G, note, r2lab="within-R2")
    out["_V"] = np.asarray(m.cov_params()) * (infl ** 2)
    out["_b"] = np.asarray(m.params)
    out["_cols"] = list(Xd.columns)
    out["_dfree"] = dfree
    return m, out


def lincom(res, weights):
    c = np.zeros(len(res["_cols"]))
    for k, w in weights.items():
        c[res["_cols"].index(k)] = w
    b = float(c @ res["_b"])
    se = float(np.sqrt(c @ res["_V"] @ c))
    return b, se, 2 * stats.t.sf(abs(b / se), res["_dfree"])


def waldtest(res, names):
    R = np.zeros((len(names), len(res["_cols"])))
    for i, k in enumerate(names):
        R[i, res["_cols"].index(k)] = 1.0
    Rb = R @ res["_b"]
    RVR = R @ res["_V"] @ R.T
    W = float(Rb @ np.linalg.solve(RVR, Rb)) / len(names)
    return W, len(names), res["_dfree"], stats.f.sf(W, len(names), res["_dfree"])


print("\n--- 4. ESTIMATES ------------------------------------------------------")

need1 = ["used_share", "zerodebt", "young", "age", "lnA", "small", "medplus",
         "wave", "sector", "prov", "fid"]
d1 = df.dropna(subset=need1)
f1 = ("used_share ~ zerodebt + young + zerodebt:young + lnA + age + small + "
      "medplus + C(wave) + C(sector) + C(prov)")
m1, RES["(1) pooled"] = pooled(
    d1, f1, ["zerodebt", "young", "zerodebt:young", "lnA", "age", "small", "medplus"],
    "BETWEEN-firm: sector+year+province FE, NO firm FE (the paper's eq. 1)")

m2, RES["(2) firm FE"] = firm_fe(
    dF, ["zerodebt", "young", "zd_young", "lnA", "small", "medplus"],
    "WITHIN-firm: firm FE + year FE; Age drops (= t - birthyear, collinear); "
    "sector/province absorbed")

m3, RES["(3) firm FE, lean"] = firm_fe(
    dF, ["zerodebt", "young", "zd_young"],
    "WITHIN-firm, lean: firm FE + year FE only; drops ln(Assets) and size, which "
    "are themselves outcomes of firm growth (bad controls)")

needA = needF + ["a02", "a35", "a610"]
dA = df.dropna(subset=needA).copy()
dA = dA[dA.groupby("fid")["fid"].transform("size") >= 2]
m4, RES["(4) firm FE, age bins"] = firm_fe(
    dA, ["a02", "a35", "a610", "zerodebt", "lnA", "small", "medplus"],
    "WITHIN-firm age profile, omitted bin = 10+ years; identified only off the "
    "NONLINEARITY in age, since linear age is absorbed by firm+year FE")

d5 = dF[dF.groupby("fid")["fid"].transform("size") >= 3]
m5, RES["(5) firm FE, >=3 waves"] = firm_fe(
    d5, ["zerodebt", "young", "zd_young", "lnA", "small", "medplus"],
    "WITHIN-firm, firms observed in >=3 waves only")

d6 = df.dropna(subset=needF + ["importer"]).copy()
d6 = d6[d6.groupby("fid")["fid"].transform("size") >= 2]
m6, RES["(6) firm FE + import"] = firm_fe(
    d6, ["zerodebt", "young", "zd_young", "lnA", "importer", "small", "medplus"],
    "WITHIN-firm, waves 2009/2011/2013 where imports are recorded")

d7 = df.dropna(subset=["used_share", "zerodebt", "young_raw", "lnA", "small",
                       "medplus", "wave", "fid"]).copy()
d7 = d7[d7.groupby("fid")["fid"].transform("size") >= 2]
d7["zd_young"] = d7["zerodebt"] * d7["young_raw"]
d7["young"] = d7["young_raw"]
m7, RES["(7) firm FE, raw age"] = firm_fe(
    d7, ["zerodebt", "young", "zd_young", "lnA", "small", "medplus"],
    "WITHIN-firm using the RAW self-reported establishment year (birth year NOT "
    "pinned, so 315 firm-pairs cross the age-5 line backwards). Reported for "
    "transparency: it is LARGER than col (2), not attenuated, so pinning the "
    "birth year is a conservative choice, not a result-manufacturing one")

ROWKEY = ["zerodebt", "young", "zerodebt:young", "zd_young", "a02", "a35", "a610",
          "lnA", "age", "importer", "small", "medplus"]
ROWLAB = ["ZeroDebt", "Young (age<=5)", "ZeroDebt x Young", "ZeroDebt x Young",
          "Age 0-2", "Age 3-5", "Age 6-10", "ln(Assets)", "Firm age", "Import",
          "Small (10-49)", "Medium+ (50+)"]
COLS = ["(1) pooled", "(2) firm FE", "(3) firm FE, lean", "(4) firm FE, age bins",
        "(5) firm FE, >=3 waves", "(6) firm FE + import", "(7) firm FE, raw age"]


def ptable(cols, title):
    print("\n  %s" % title)
    hdr = "  %-22s" % "" + "".join("%16s" % c.replace("firm FE", "FE") for c in cols)
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    shown = set()
    for rk, rl in zip(ROWKEY, ROWLAB):
        if not any(rk in RES[c] for c in cols):
            continue
        if rl in shown and rk == "zd_young":
            rl = ""
        shown.add(rl)
        l1 = "  %-22s" % rl
        l2 = "  %-22s" % ""
        for c in cols:
            if rk in RES[c]:
                b, se, p = RES[c][rk]
                l1 += "%16s" % ("%8.3f%s" % (b, STAR(p)))
                l2 += "%16s" % ("(%.3f)" % se)
            else:
                l1 += "%16s" % "-"
                l2 += "%16s" % ""
        print(l1)
        print(l2)
    print("  %-22s" % "N" + "".join("%16s" % format(RES[c]["_N"], ",") for c in cols))
    print("  %-22s" % "firms" + "".join("%16s" % format(RES[c]["_G"], ",") for c in cols))
    print("  %-22s" % "R2 / within-R2" + "".join("%16.3f" % RES[c]["_R2"] for c in cols))
    print("  %-22s" % "firm FE" + "".join("%16s" % ("NO" if "pooled" in c else "YES") for c in cols))
    print("  %-22s" % "year FE" + "".join("%16s" % "YES" for c in cols))
    print()
    for c in cols:
        print("  [%s] %s" % (c, RES[c]["_note"]))


ptable(COLS, "Dependent variable: share of equipment used/second-hand at purchase, "
             "in PERCENTAGE POINTS (0-100)")

b2y, s2y, p2y = RES["(2) firm FE"]["young"]
b2i, s2i, p2i = RES["(2) firm FE"]["zd_young"]
b2z, s2z, p2z = RES["(2) firm FE"]["zerodebt"]
b1y = RES["(1) pooled"]["young"][0]
b1i = RES["(1) pooled"]["zerodebt:young"][0]

bfree, sfree, pfree = lincom(RES["(2) firm FE"], {"young": 1.0, "zd_young": 1.0})
Wbin, kbin, dfbin, pbin = waldtest(RES["(4) firm FE, age bins"],
                                   ["a02", "a35", "a610"])

print("""
--- 5. HOW TO READ THIS ------------------------------------------------
  Column (1) is BETWEEN firms.  It compares different firms in the same
  sector, year and province, so it licenses only:
      "Older firms report lower used-capital shares than younger firms."

  Column (2) is WITHIN firms.  The same firm is followed across waves and
  alpha_i absorbs everything permanent about it (founder, location, sector,
  technology, its habitual way of answering the question).  A within-firm
  sentence -- "firms reduce their used-capital share as they age" -- has to be
  earned HERE, and only for the margins that actually come out significant.

  All coefficients below are in PERCENTAGE POINTS of the used-machinery share.

  Young, within firm:          %+6.2f pp (se %.2f, p = %.3f)
  Young, between firms:        %+6.2f pp
  ZeroDebt x Young, within:    %+6.2f pp (se %.2f, p = %.3f)
  ZeroDebt x Young, between:   %+6.2f pp
  ZeroDebt main effect, within:%+6.2f pp (se %.2f, p = %.3f)
  Young + ZeroDebt x Young:    %+6.2f pp (se %.2f, p = %.3f)

  VERDICT, one margin at a time:
   * Unconditional age.  Within firm, Young alone is %s (p = %.3f) and the
     three-bin age profile is jointly %s.  The unconditional sentence
     "firms reduce their used-capital share as they age" is therefore
     %s by the within-firm evidence and should NOT be written.
   * Age interacted with financing.  For DEBT-FREE firms the within-firm young
     premium is %.1f pp (p = %.3f), and %.1f pp (p = %.3f) among firms seen in
     at least three waves.  THIS is the sentence the panel supports:
     "within the same firm, the used-capital share is higher in the years when
      the firm is both young and carrying no debt, and falls as it ages out of
      that state."
   * So the referee's correction stands, and it lands on the AGE claim, not on
     the financing claim -- which is the claim the model actually needs.
""" % (b2y, s2y, p2y, b1y,
       b2i, s2i, p2i, b1i, b2z, s2z, p2z,
       bfree, sfree, pfree,
       "significant" if p2y < .10 else "INSIGNIFICANT", p2y,
       "significant" if pbin < .10 else "INSIGNIFICANT",
       "supported" if (p2y < .10 or pbin < .10) else "NOT supported",
       bfree, pfree,
       lincom(RES["(5) firm FE, >=3 waves"], {"young": 1.0, "zd_young": 1.0})[0],
       lincom(RES["(5) firm FE, >=3 waves"], {"young": 1.0, "zd_young": 1.0})[2]))

print("--- 6. TESTS THE PAPER SHOULD REPORT ---------------------------------")
for lab, res in [("(2) firm FE", RES["(2) firm FE"]),
                 ("(5) firm FE, >=3 waves", RES["(5) firm FE, >=3 waves"])]:
    b, se, p = lincom(res, {"young": 1.0, "zd_young": 1.0})
    print("  %-22s Young + ZeroDebt x Young  (the total within-firm effect of "
          "being young FOR A DEBT-FREE FIRM):" % lab)
    print("  %-22s   %+.2f pp   se %.2f   p = %.3f" % ("", b, se, p))

print("\n  (4) joint Wald test that the whole within-firm age profile is flat")
print("      H0: Age0-2 = Age3-5 = Age6-10 = 0")
print("      F(%d, %d) = %.2f, p = %.3f -> %s" % (kbin, dfbin, Wbin, pbin,
      "REJECT: the profile declines within firm" if pbin < .10 else
      "CANNOT reject: NO within-firm age profile once alpha_i is absorbed"))
bins = [RES["(4) firm FE, age bins"][k_][0] for k_ in ("a02", "a35", "a610")]
print("      point estimates, relative to age 10+: %.2f / %.2f / %.2f / 0.00 pp "
      "(monotone: %s)" % (bins[0], bins[1], bins[2],
                          "yes" if bins[0] >= bins[1] >= bins[2] >= 0 else "no"))

W2, k2, dfr2, p2 = waldtest(RES["(2) firm FE"], ["young", "zd_young"])
print("\n  (2) joint Wald test that the young margin is absent altogether")
print("      H0: Young = ZeroDebt x Young = 0")
print("      F(%d, %d) = %.2f, p = %.3f -> %s" % (k2, dfr2, W2, p2,
      "REJECT: something happens at the young margin within firm" if p2 < .10
      else "cannot reject"))
print()

print("""  Caveats to state in the paper, in this order:
   1. Identification is off the firms listed in section 3 above -- those that
      cross age five between waves and those whose debt status changes.  Firms
      with no such change contribute nothing to b2 and b3.
   2. Age itself cannot appear: Age_it = t - BirthYear_i is an exact linear
      combination of the firm and year effects.  Column (4) recovers a
      nonlinear age profile, but only relative to the linear trend the fixed
      effects already absorb -- do not read its bins as raw age differences.
   3. ln(Assets) and the size dummies are themselves outcomes of firm growth.
      Column (3) drops them; report it so the reader can see the difference.
   4. The 2015 accounts report total, not physical, assets; the year effects
      absorb the level break but not any change in its dispersion.
   5. Standard errors here are d.o.f.-corrected for the absorbed firm effects
      (the areg/reghdfe convention).  Without that correction they are about
      %.0f%% smaller and the within estimates look more precise than they are.

  SUPERSEDES the firm-FE column of motivation_analysis.py.  That column reported
  Young = +3.3 pp*** within firm, but it (a) omitted the ZeroDebt x Young
  interaction, so the debt-free young premium loaded onto Young alone, (b) used
  the raw self-reported establishment year, and (c) did not correct the d.o.f.
  for the absorbed firm effects.  With all three fixed, Young alone is +0.8 pp
  (n.s.) and the effect sits on the interaction instead.  Do not quote the
  +3.3 pp figure.
""" % (100 * (1 - 1 / np.sqrt(len(dF) / float(len(dF) - dF.fid.nunique())))))

print("=" * 78)
print("  7. AGE BETWEEN FIRMS, SIZE WITHIN FIRMS")
print("=" * 78)

print("""
  PANEL A -- AGE, BETWEEN FIRMS (the only age evidence six waves can support)
""")
raw = df.dropna(subset=["used_share", "age"]).copy()
raw["bin"] = pd.cut(raw["age"], [-0.5, 2.5, 5.5, 10.5, 200],
                    labels=["0-2", "3-5", "6-10", "10+"])
tt = raw.groupby("bin", observed=True)["used_share"].agg(["mean", "count"])
print("  raw means (all firms, zeros included):")
for b_, r_ in tt.iterrows():
    print("    age %-5s %6.1f%%   N = %s" % (b_, r_["mean"], format(int(r_["count"]), ",")))

needAge = ["used_share", "a02", "a35", "a610", "small", "medplus", "lnA",
           "wave", "sector", "prov", "fid"]
dAge = df.dropna(subset=needAge)
mAge, RES["age between"] = pooled(
    dAge, "used_share ~ a02 + a35 + a610 + small + medplus + lnA + C(wave) "
          "+ C(sector) + C(prov)",
    ["a02", "a35", "a610", "small", "medplus", "lnA"],
    "BETWEEN-firm age profile, omitted bin = 10+ years; sector+year+province FE")
print("\n  conditional on size, assets, sector, year and province "
      "(omitted bin: older than 10):")
for k_, lab_ in [("a02", "age 0-2"), ("a35", "age 3-5"), ("a610", "age 6-10")]:
    b_, s_, p_ = RES["age between"][k_]
    print("    %-9s %+6.2f pp  (se %.2f, p = %.3f)%s" % (lab_, b_, s_, p_, STAR(p_)))

print("""
  PANEL B -- SIZE, WITHIN FIRMS
  The calibration moment beta_size is the slope of the used share on LOG ASSETS
  among firms with a POSITIVE used share (betasize_regression.py, 2013).  Below
  is that same gradient with firm fixed effects, i.e. asked of the same firm as
  its assets grow.
""")

d13 = df[(df.wave == 2013) & (df.used_share > 0)].dropna(subset=["used_share", "lnA"])
m13 = smf.ols("used_share ~ lnA", d13).fit(cov_type="HC1")
print("  (a) AS CALIBRATED  2013 cross-section, used>0, no controls, no FE")
print("      lnA %+6.3f pp per log point (se %.3f, p = %.3f)%s   N = %s"
      % (m13.params["lnA"], m13.bse["lnA"], m13.pvalues["lnA"],
         STAR(m13.pvalues["lnA"]), format(int(m13.nobs), ",")))
print("      [the calibration target is -2.69 pp; reproduced here as %+.2f]"
      % m13.params["lnA"])

dIntAll = df[df.used_share > 0].dropna(subset=["used_share", "lnA", "wave", "fid"])
mIntP, RES["size between"] = pooled(
    dIntAll, "used_share ~ lnA + C(wave)", ["lnA"],
    "BETWEEN-firm, all waves, used>0, year FE")
b_, s_, p_ = RES["size between"]["lnA"]
print("\n  (b) BETWEEN FIRMS   all waves, used>0, year FE")
print("      lnA %+6.3f pp per log point (se %.3f, p = %.3f)%s   N = %s, firms %s"
      % (b_, s_, p_, STAR(p_), format(RES["size between"]["_N"], ","),
         format(RES["size between"]["_G"], ",")))

dInt = dIntAll[dIntAll.groupby("fid")["fid"].transform("size") >= 2].copy()
mIntW, RES["size within"] = firm_fe(
    dInt, ["lnA"],
    "WITHIN-firm, used>0, firm FE + year FE -- the calibration moment asked of "
    "the same firm as its assets grow")
b_, s_, p_ = RES["size within"]["lnA"]
print("\n  (c) WITHIN FIRMS    used>0, firm FE + year FE")
print("      lnA %+6.3f pp per log point (se %.3f, p = %.3f)%s   N = %s, firms %s"
      % (b_, s_, p_, STAR(p_), format(RES["size within"]["_N"], ","),
         format(RES["size within"]["_G"], ",")))

dFull = df.dropna(subset=["used_share", "lnA", "wave", "fid"]).copy()
dFull = dFull[dFull.groupby("fid")["fid"].transform("size") >= 2]
mFullW, RES["size within, zeros"] = firm_fe(
    dFull, ["lnA"],
    "WITHIN-firm, ALL firms (zeros kept), firm FE + year FE -- avoids selecting "
    "on a function of the dependent variable")
b_, s_, p_ = RES["size within, zeros"]["lnA"]
print("\n  (d) WITHIN FIRMS    all firms incl. zeros, firm FE + year FE")
print("      lnA %+6.3f pp per log point (se %.3f, p = %.3f)%s   N = %s, firms %s"
      % (b_, s_, p_, STAR(p_), format(RES["size within, zeros"]["_N"], ","),
         format(RES["size within, zeros"]["_G"], ",")))

dEmp = df.dropna(subset=["used_share", "emp", "wave", "fid"]).copy()
dEmp["lnE"] = np.log(dEmp["emp"])
dEmp = dEmp[(dEmp.used_share > 0)]
dEmp = dEmp[dEmp.groupby("fid")["fid"].transform("size") >= 2]
mEmpW, RES["size within, emp"] = firm_fe(
    dEmp, ["lnE"],
    "WITHIN-firm, used>0, size measured by log EMPLOYMENT rather than assets")
b_, s_, p_ = RES["size within, emp"]["lnE"]
print("\n  (e) WITHIN FIRMS    used>0, log EMPLOYMENT instead of assets")
print("      lnE %+6.3f pp per log point (se %.3f, p = %.3f)%s   N = %s, firms %s"
      % (b_, s_, p_, STAR(p_), format(RES["size within, emp"]["_N"], ","),
         format(RES["size within, emp"]["_G"], ",")))

bw = RES["size within"]["lnA"]
bb = RES["size between"]["lnA"]
print("""
  VERDICT ON THE SIZE GRADIENT
    between firms: %+.3f pp per log point of assets (p = %.3f)
    within  firms: %+.3f pp per log point of assets (p = %.3f)
    -> the gradient %s within firm.
       %s
""" % (bb[0], bb[2], bw[0], bw[2],
       "SURVIVES" if (bw[2] < .10 and np.sign(bw[0]) == np.sign(bb[0]))
       else "does NOT survive",
       ("The calibration moment is a genuine within-firm relationship: the same "
        "firm sheds used equipment as it grows."
        if (bw[2] < .10 and np.sign(bw[0]) == np.sign(bb[0]))
        else "beta_size is a CROSS-SECTIONAL fact. Small firms use more used "
             "equipment than big firms, but a given firm growing bigger does not "
             "measurably shed it over ten years. State the moment as a "
             "cross-sectional target, which is how the model uses it.")))

print("""  Caveat on (c): conditioning on used>0 selects on a function of the
  dependent variable, and firms move in and out of that sample between waves.
  Column (d) keeps the zeros and is the cleaner within-firm object; the
  calibration moment conditions on used>0 only because the model has no
  extensive margin.
""")

print("=" * 78)
print("  8. HEADLINE: THE COLLATERAL CHANNEL, WITHIN FIRM")
print("=" * 78)

dH = df.dropna(subset=["used_share", "zerodebt", "lnA", "small", "medplus",
                       "wave", "fid"]).copy()
dH = dH[dH.groupby("fid")["fid"].transform("size") >= 2]
lnA_bar = dH["lnA"].mean()
dH["lnAc"] = dH["lnA"] - lnA_bar
dH["zd_lnAc"] = dH["zerodebt"] * dH["lnAc"]

mH, RES["headline"] = firm_fe(
    dH, ["zerodebt", "lnAc", "zd_lnAc", "small", "medplus"],
    "WITHIN-firm collateral channel: firm FE + year FE, ln(Assets) centred at "
    "the sample mean so ZeroDebt is read at average assets")

print("\n  UsedShare_it = alpha_i + lambda_t + b1 ZeroDebt_it + b2 lnA_it")
print("                 + b3 (ZeroDebt x lnA)_it + Gamma' Size_it + eps_it\n")
for k_, lab_ in [("zerodebt", "ZeroDebt (at mean assets)"),
                 ("lnAc", "ln(Assets), centred"),
                 ("zd_lnAc", "ZeroDebt x ln(Assets)"),
                 ("small", "Small (10-49)"),
                 ("medplus", "Medium+ (50+)")]:
    b_, s_, p_ = RES["headline"][k_]
    print("    %-26s %+7.3f pp  (se %.3f, p = %.3f)%s" % (lab_, b_, s_, p_, STAR(p_)))
print("    %-26s %s obs, %s firms, within-R2 %.3f"
      % ("", format(RES["headline"]["_N"], ","),
         format(RES["headline"]["_G"], ","), RES["headline"]["_R2"]))

q = dH["lnA"].quantile([.10, .25, .50, .75, .90])
print("\n  Implied ZeroDebt premium at points of the asset distribution")
print("  (b1 + b3 * (lnA - mean lnA); the collateral channel in one line):")
for lab_, qq in zip(["p10", "p25", "p50", "p75", "p90"], q):
    b_, s_, p_ = lincom(RES["headline"], {"zerodebt": 1.0, "zd_lnAc": qq - lnA_bar})
    print("    %-4s (lnA = %5.2f)  %+6.2f pp  (se %.2f, p = %.3f)%s"
          % (lab_, qq, b_, s_, p_, STAR(p_)))

dHp = df.dropna(subset=["used_share", "zerodebt", "lnA", "small", "medplus",
                        "wave", "sector", "prov", "fid"]).copy()
dHp["lnAc"] = dHp["lnA"] - lnA_bar
mHp, RES["headline pooled"] = pooled(
    dHp, "used_share ~ zerodebt + lnAc + zerodebt:lnAc + small + medplus "
         "+ C(wave) + C(sector) + C(prov)",
    ["zerodebt", "lnAc", "zerodebt:lnAc", "small", "medplus"],
    "BETWEEN-firm counterpart: sector+year+province FE, no firm FE")
bp_, sp_, pp_ = RES["headline pooled"]["zerodebt:lnAc"]
print("\n  between-firm counterpart of the interaction: %+.3f pp (se %.3f, "
      "p = %.3f)%s" % (bp_, sp_, pp_, STAR(pp_)))

HROWS = [(("zerodebt",), r"Zero debt"),
         (("lnAc",), r"$\ln$(Assets), centred"),
         (("zerodebt:lnAc", "zd_lnAc"), r"Zero debt $\times\log(\mathrm{Assets})$"),
         (("small",), r"Small (10--49)"), (("medplus",), r"Medium+ (50+)")]
HCOLS = [("headline pooled", "(1) Between firms"), ("headline", "(2) Within firms")]
hl = [r"\begin{table}[H]", r"\centering",
      r"\caption{Credit Exclusion, Collateral Capacity and Used Capital}",
      r"\label{tab:collateral}", r"\begin{tabular}{lcc}", r"\toprule",
      " & " + " & ".join(h for _, h in HCOLS) + r" \\", r"\midrule"]
for keys, rl in HROWS:
    c1, c2 = [rl], [""]
    for ck, _ in HCOLS:
        hit = next((k for k in keys if k in RES[ck]), None)
        if hit is None:
            c1.append("---")
            c2.append("")
            continue
        b_, s_, p_ = RES[ck][hit]
        c1.append("$%.3f^{%s}$" % (b_, STAR(p_)) if STAR(p_) else "$%.3f$" % b_)
        c2.append("$(%.3f)$" % s_)
    hl.append(" & ".join(c1) + r" \\")
    hl.append(" & ".join(c2) + r" \\")
hl.append(r"\addlinespace")
hl.append(r"\multicolumn{3}{l}{\emph{Implied zero-debt premium, by asset percentile}} \\")
for lab_, qq in zip(["p10", "p25", "p50", "p75", "p90"], q):
    b_, s_, p_ = lincom(RES["headline"], {"zerodebt": 1.0, "zd_lnAc": qq - lnA_bar})
    bp2, sp2, pp2 = (RES["headline pooled"]["zerodebt"][0]
                     + RES["headline pooled"]["zerodebt:lnAc"][0] * (qq - lnA_bar),
                     float("nan"), float("nan"))
    hl.append("\\quad %s & $%.2f$ & $%.2f^{%s}$ \\\\" % (lab_, bp2, b_, STAR(p_))
              if STAR(p_) else "\\quad %s & $%.2f$ & $%.2f$ \\\\" % (lab_, bp2, b_))
hl += [r"\addlinespace",
       "Observations & " + " & ".join("$%s$" % texnum(RES[c]["_N"])
                                      for c, _ in HCOLS) + r" \\",
       "Firms & " + " & ".join("$%s$" % texnum(RES[c]["_G"])
                               for c, _ in HCOLS) + r" \\",
       "$R^{2}$ (within for FE) & " + " & ".join("$%.3f$" % RES[c]["_R2"]
                                                 for c, _ in HCOLS) + r" \\",
       r"Firm FE & No & Yes \\", r"Year FE & Yes & Yes \\",
       r"Sector, province FE & Yes & Absorbed \\",
       r"\bottomrule", r"\end{tabular}", r"\end{table}"]
OUT_TEX2 = os.path.join(HERE, "collateral_channel_table.tex")
with open(OUT_TEX2, "w", encoding="utf-8") as fh:
    fh.write("\n".join(hl) + "\n")

MROWS = [(("zerodebt:lnAc", "zd_lnAc"), r"Zero debt $\times$ log assets"),
         (("zerodebt:young", "zd_young"), r"Zero debt $\times$ young"),
         (("zerodebt",), r"Zero debt"),
         (("lnA", "lnAc"), r"Log assets"),
         (("young",), r"Young (age $\le 5$)"),
         (("small",), r"Small (10--49)"),
         (("medplus",), r"Medium+ (50+)")]
MCOLS = [("headline", "(1) Asset", "interaction"),
         ("(2) firm FE", "(2) Young-firm", "interaction")]
MTOT = {"headline": {"lnAc": 1.0, "zd_lnAc": 1.0},
        "(2) firm FE": {"young": 1.0, "zd_young": 1.0}}

MPLAIN = ["Zero debt x log assets", "Zero debt x young", "Zero debt",
          "Log assets", "Young (age<=5)", "Small (10-49)", "Medium+ (50+)"]

print("=" * 78)
print("  THE TWO MODELS THE PAPER KEEPS  (both WITHIN firm)")
print("=" * 78)
print("""
  Dependent variable: used-capital share, percentage points.
  These are TWO SEPARATE regressions, not one.  A dash means the variable is not
  in that regression at all:
    (1) interacts Zero debt with log assets  -- the collateral margin
    (2) interacts Zero debt with young       -- the age margin
  Both include firm effects, year effects and the size dummies; sector and
  province are absorbed by the firm effect.  Standard errors in parentheses,
  clustered by firm and d.o.f.-corrected for the absorbed firm effects.
""")
hdr = "  %-26s" % "" + "".join("%18s" % h1 for _, h1, _ in MCOLS)
print(hdr)
print("  %-26s" % "" + "".join("%18s" % h2 for _, _, h2 in MCOLS))
print("  " + "-" * (len(hdr) - 2))
for (keys, _), rl in zip(MROWS, MPLAIN):
    l1 = "  %-26s" % rl
    l2 = "  %-26s" % ""
    for ck, _, _ in MCOLS:
        hit = next((k for k in keys if k in RES[ck]), None)
        if hit is None:
            l1 += "%18s" % "---"
            l2 += "%18s" % ""
            continue
        b_, s_, p_ = RES[ck][hit]
        l1 += "%18s" % ("%8.3f%s" % (b_, STAR(p_)))
        l2 += "%18s" % ("(%.3f)" % s_)
    print(l1)
    print(l2)
print("  " + "-" * (len(hdr) - 2))
for lab_, key_ in [("Observations", "_N"), ("Firms", "_G")]:
    print("  %-26s" % lab_ + "".join("%18s" % format(RES[c][key_], ",")
                                     for c, _, _ in MCOLS))
print("  %-26s" % "within R2" + "".join("%18.3f" % RES[c]["_R2"]
                                        for c, _, _ in MCOLS))
for lab_ in ["Firm fixed effects", "Year fixed effects"]:
    print("  %-26s" % lab_ + "".join("%18s" % "Yes" for _ in MCOLS))
print("  %-26s" % "Sector, province" + "".join("%18s" % "Absorbed" for _ in MCOLS))

ml = [r"\begin{table}[H]", r"\centering", r"\small",
      r"\setlength{\tabcolsep}{8pt}",
      r"\caption{Financing Frictions and Used Capital}",
      r"\label{tab:emp_panel}", r"\begin{tabular}{lcc}", r"\toprule",
      r"\multicolumn{3}{l}{\emph{Dependent variable:} used-capital share "
      r"(percentage points)} \\",
      r"\addlinespace[2pt]",
      " & " + " & ".join(h1 for _, h1, _ in MCOLS) + r" \\",
      " & " + " & ".join(h2 for _, _, h2 in MCOLS) + r" \\", r"\midrule"]
for keys, rl in MROWS:
    c1, c2 = [rl], [""]
    for ck, _, _ in MCOLS:
        hit = next((k for k in keys if k in RES[ck]), None)
        if hit is None:
            c1.append("---")
            c2.append("")
            continue
        b_, s_, p_ = RES[ck][hit]
        c1.append("$%.3f^{%s}$" % (b_, STAR(p_)) if STAR(p_) else "$%.3f$" % b_)
        c2.append("$(%.3f)$" % s_)
    ml.append(" & ".join(c1) + r" \\")
    ml.append(" & ".join(c2) + r" \\")

print("\n  beta_2 + beta_3, the total within-firm gradient for a no-debt firm")
print("  (computed for reference; NOT reported as a row of the table):")
for ck, h1, _ in MCOLS:
    b_, s_, p_ = lincom(RES[ck], MTOT[ck])
    print("    %-16s %+7.3f pp  (se %.3f, p = %.3f)%s"
          % (h1, b_, s_, p_, STAR(p_)))
ml += [r"\midrule",
       "Observations & " + " & ".join("$%s$" % texnum(RES[c]["_N"])
                                      for c, _, _ in MCOLS) + r" \\",
       "Firms & " + " & ".join("$%s$" % texnum(RES[c]["_G"])
                               for c, _, _ in MCOLS) + r" \\",
       "$R^{2}$ (within) & " + " & ".join("$%.3f$" % RES[c]["_R2"]
                                          for c, _, _ in MCOLS) + r" \\",
       r"Firm effects & Yes & Yes \\",
       r"Year effects & Yes & Yes \\",
       r"Sector, province & Absorbed & Absorbed \\",
       r"\bottomrule", r"\end{tabular}", r"\end{table}"]
with open(os.path.join(HERE, "main_table.tex"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(ml) + "\n")

bH, sH, pH = RES["headline"]["zd_lnAc"]
print("""
  READING
    The interaction is %+.3f pp per log point (p = %.3f)%s.  Within the same
    firm, being cut off from credit raises the used-capital share most when the
    firm has little to pledge, and that premium %s as it accumulates assets.
    This is the collateral constraint of the model, measured directly, and it
    does not rely on firm age -- which is what six waves over ten years cannot
    identify anyway.
""" % (bH, pH, STAR(pH),
       "shrinks" if bH < 0 else "GROWS (wrong sign for the mechanism)"))

TEXCOLS = ["(1) pooled", "(2) firm FE", "(3) firm FE, lean", "(4) firm FE, age bins"]
TEXHEAD = ["(1) Pooled", "(2) Firm FE", "(3) Firm FE, lean", "(4) Firm FE, age bins"]
TEXROWS = [(("zerodebt",), r"Zero debt"), (("young",), r"Young (age $\le 5$)"),
           (("zerodebt:young", "zd_young"), r"Zero debt $\times$ Young"),
           (("a02",), r"Age 0--2"), (("a35",), r"Age 3--5"),
           (("a610",), r"Age 6--10"),
           (("lnA",), r"$\ln$(Assets)"), (("age",), r"Firm age"),
           (("small",), r"Small (10--49)"), (("medplus",), r"Medium+ (50+)")]

lines = [r"\begin{table}[H]", r"\centering",
         r"\caption{Used Capital, Financing and Firm Age: Between- and Within-Firm Estimates}",
         r"\label{tab:emp_firmfe}",
         r"\begin{tabular}{l" + "c" * len(TEXCOLS) + "}", r"\toprule",
         " & " + " & ".join(TEXHEAD) + r" \\", r"\midrule"]
for keys, rl in TEXROWS:
    if not any(k in RES[c] for k in keys for c in TEXCOLS):
        continue
    c1, c2 = [rl], [""]
    for c in TEXCOLS:
        hit = next((k for k in keys if k in RES[c]), None)
        if hit is None:
            c1.append("---")
            c2.append("")
            continue
        b, se, p = RES[c][hit]
        c1.append("$%.3f^{%s}$" % (b, STAR(p)) if STAR(p) else "$%.3f$" % b)
        c2.append("$(%.3f)$" % se)
    lines.append(" & ".join(c1) + r" \\")
    lines.append(" & ".join(c2) + r" \\")
lines += [r"\addlinespace",
          "Observations & " + " & ".join("$%s$" % texnum(RES[c]["_N"]) for c in TEXCOLS) + r" \\",
          "Firms & " + " & ".join("$%s$" % texnum(RES[c]["_G"]) for c in TEXCOLS) + r" \\",
          "$R^2$ (within for FE) & " + " & ".join("$%.3f$" % RES[c]["_R2"] for c in TEXCOLS) + r" \\",
          "Firm FE & " + " & ".join("No" if "pooled" in c else "Yes" for c in TEXCOLS) + r" \\",
          "Year FE & " + " & ".join("Yes" for _ in TEXCOLS) + r" \\",
          "Sector, province FE & " + " & ".join(
              "Yes" if "pooled" in c else "Absorbed" for c in TEXCOLS) + r" \\",
          r"\bottomrule", r"\end{tabular}", r"\end{table}"]
with open(OUT_TEX, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")

print("=" * 78)
print("  9. FILES WRITTEN BY THIS RUN")
print("=" * 78)
GENERATED = [(os.path.join(HERE, "main_table.tex"),
              "Table 1 of the paper (tab:emp_panel) -- within-firm estimates"),
             (OUT_TEX, "between/within comparison, all seven specifications"),
             (OUT_TEX2, "collateral channel with the asset-percentile premia")]
for path_, what_ in GENERATED:
    print("\n  %s\n  (%s)\n  %s" % (os.path.basename(path_), what_, "-" * 74))
    with open(path_, encoding="utf-8") as fh:
        for ln in fh.read().splitlines():
            print("  | %s" % ln)
print("\n  %s   (this file: the full log above)" % os.path.basename(OUT_TXT))
print("\nDONE.")

sys.stdout = _tee.stdout
_tee.close()
print("results saved to %s" % OUT_TXT)
