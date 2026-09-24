import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
OUT = os.path.join(HERE, "credit_constraints_results.txt")

lines = []


def say(s=""):
    print(s)
    lines.append(s)


R = {
    2011: dict(used="q30ab_11", new="q30aa_11", slf="q30ac_11", mach="q28_11",
               scale=100, applied="q80_11", problem="q81_11", whynot="q87_11",
               emp="EAq1x10", assets="EAq1q10", debt="EAq1v10"),
    2013: dict(used="q30ab_13", new="q30aa_13", slf="q30ac_13", mach="q28_13",
               scale=100, applied="q80_13", problem="q81_13", whynot="q87_13",
               emp="EAq1x12", assets="EAq1q12", debt="EAq1v12"),
    2015: dict(used="q31ab_15", new="q31aa_15", slf="q31ac_15", mach="q29_15",
               scale=100, applied="q70_15", problem="q71_15", whynot="q77_15",
               emp="EAq1k14", assets="EAq1h14", debt="EAq1i14"),
}

frames = []
whynot_tab = {}
for yr, m in sorted(R.items()):
    path = os.path.join(DATA, "enterprise %d.dta" % yr)
    cols = ["id"] + [m[k] for k in
                     ("used", "new", "slf", "mach", "applied", "problem",
                      "whynot", "emp", "assets", "debt")]
    d = pd.read_stata(path, convert_categoricals=False, columns=cols)
    d = d[~d["id"].isna()].drop_duplicates(subset="id", keep="first")
    num = lambda c: pd.to_numeric(d[c], errors="coerce")

    sc = float(m["scale"])
    u, n, s = num(m["used"]), num(m["new"]), num(m["slf"])
    tot = pd.concat([n, u, s], axis=1).sum(axis=1, min_count=3)
    in_block = (num(m["mach"]).between(2, 5) & n.notna() & u.notna() & s.notna()
                & tot.between(0.99 * sc, 1.01 * sc))

    w = pd.DataFrame(dict(fid=d["id"].astype(float).values, wave=yr))
    w["used_pp"] = (100.0 * (u / sc).clip(0, 1)).values
    w["in_block"] = in_block.values
    ap = num(m["applied"])
    w["applied"] = np.where(ap.isin([1]), 1.0, np.where(ap.isin([0, 2]), 0.0, np.nan))
    pr = num(m["problem"])
    w["problem"] = np.where(pr.isin([1]), 1.0, np.where(pr.isin([0, 2]), 0.0, np.nan))
    w["whynot"] = num(m["whynot"]).values
    emp = num(m["emp"]).where(lambda x: x > 0)
    w["size"] = pd.cut(emp, [0, 9.5, 49.5, np.inf],
                       labels=["micro", "small", "medium"]).values
    a = num(m["assets"]).where(lambda x: x > 0)
    w["lnA"] = np.log(a).values
    dbt = num(m["debt"]).where(lambda x: x >= 0)
    w["zerodebt"] = (dbt == 0).astype(float).where(dbt.notna()).values
    frames.append(w)

    lab = pd.read_stata(path.replace("enterprise_data", "enterprise_data"),
                        convert_categoricals=False, columns=[m["whynot"]])
    whynot_tab[yr] = num(m["whynot"]).value_counts(dropna=True).sort_index()

P = pd.concat(frames, ignore_index=True)

say("=" * 78)
say("CREDIT ACCESS AND USED CAPITAL  --  Vietnam SME survey, 2011/2013/2015")
say("=" * 78)

WHY = {1: "Inadequate collateral", 2: "Don't want to incur debt",
       3: "Process too difficult", 4: "Didn't need one",
       5: "Interest rate too high", 6: "Already heavily indebted",
       7: "Other"}
say("")
say("--- 1. STATED REASON FOR NOT APPLYING FOR A FORMAL LOAN ----------------")
say("  (q87_11, q87_13, q77_15; firms that did not apply since the last round)")
say("")
say("  %-28s %9s %9s %9s %9s" % ("reason", "2011", "2013", "2015", "pooled"))
say("  " + "-" * 68)
tot_by_yr = {y: whynot_tab[y][whynot_tab[y].index.isin(WHY)].sum() for y in R}
pooled_tot = sum(tot_by_yr.values())
for k, name in WHY.items():
    cells = []
    psum = 0
    for y in sorted(R):
        c = float(whynot_tab[y].get(k, 0))
        psum += c
        cells.append("%.1f%%" % (100 * c / tot_by_yr[y]))
    say("  %-28s %9s %9s %9s %9s"
        % (name, cells[0], cells[1], cells[2], "%.1f%%" % (100 * psum / pooled_tot)))
say("  " + "-" * 68)
say("  %-28s %9s %9s %9s %9s"
    % ("N answering", "{:,}".format(int(tot_by_yr[2011])),
       "{:,}".format(int(tot_by_yr[2013])), "{:,}".format(int(tot_by_yr[2015])),
       "{:,}".format(int(pooled_tot))))
coll = sum(float(whynot_tab[y].get(1, 0)) for y in R)
say("")
say("  READING.  Inadequate collateral is %.1f%% of stated reasons pooled"
    % (100 * coll / pooled_tot))
say("  (%d firms of %d).  The modal answer is 'didn't need one'.  The draft's"
    % (int(coll), int(pooled_tot)))
say("  sentence -- 'a large share of those that never applied ... lacked")
say("  collateral' -- is NOT supported and must be rewritten.  What the data")
say("  do support: most non-applicants are self-financing by choice or price")
say("  (no need, unwilling to take on debt, rate too high), and a further")
say("  %.1f%% are turned away by the process or by existing debt."
    % (100 * sum(float(whynot_tab[y].get(k, 0)) for y in R for k in (3, 6)) / pooled_tot))

say("")
say("--- 2. USED SHARE AMONG APPLICANTS, BY REPORTED DIFFICULTY --------------")
say("  (q81/q71 'problems getting the loan', asked of applicants only)")
say("")
A = P[(P["applied"] == 1) & P["problem"].notna() & P["used_pp"].notna()].copy()


def gap(df, label, controls):
    f = "used_pp ~ problem" + controls
    r = smf.ols(f, data=df).fit(cov_type="HC1")
    b, se, p = r.params["problem"], r.bse["problem"], r.pvalues["problem"]
    stars = "***" if p < .01 else "**" if p < .05 else "*" if p < .10 else ""
    say("  %-46s %+6.2f pp (se %.2f)%-3s  N=%s"
        % (label, b, se, stars, "{:,}".format(int(r.nobs))))
    return b, se, int(r.nobs)


for yr in sorted(R):
    s = A[A["wave"] == yr]
    gap(s, "%d, raw difference" % yr, "")
for yr in sorted(R):
    s = A[A["wave"] == yr].dropna(subset=["lnA", "size"])
    gap(s, "%d, + log assets and size class" % yr, " + lnA + C(size)")
gap(A, "pooled 2011-2015, raw difference", " + C(wave)")
gap(A.dropna(subset=["lnA", "size"]),
    "pooled 2011-2015, + log assets and size class", " + lnA + C(size) + C(wave)")
say("")
say("  Same, on the equipment block's own universe (machinery users with a")
say("  complete new/used/self triple):")
AB = A[A["in_block"]]
for yr in sorted(R):
    gap(AB[AB["wave"] == yr], "%d, raw difference" % yr, "")
gap(AB, "pooled 2011-2015, raw difference", " + C(wave)")
say("")
say("  Share of applicants reporting difficulty: " +
    ", ".join("%d %.1f%% (N=%s)"
              % (y, 100 * A.loc[A.wave == y, "problem"].mean(),
                 "{:,}".format(int((A.wave == y).sum()))) for y in sorted(R)))

say("")
say("--- 3. THE USED SHARE ON BOTH UNIVERSES --------------------------------")
say("  %-34s %10s %10s" % ("", "raw item", "in block"))
for yr in sorted(R):
    s = P[P["wave"] == yr]
    say("  %-34s %10s %10s"
        % ("%d  N" % yr, "{:,}".format(int(s["used_pp"].notna().sum())),
           "{:,}".format(int((s["in_block"] & s["used_pp"].notna()).sum()))))
    say("  %-34s %10.1f %10.1f"
        % ("     mean used share, pp", s["used_pp"].mean(),
           s.loc[s["in_block"], "used_pp"].mean()))
say("=" * 78)

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("\n[saved] %s" % OUT)
