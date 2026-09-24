import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "enterprise_data"))
OUT = os.path.join(HERE, "knku_aggregate_results.txt")

ROUNDS = {
    2005: ("q30aa_05", "q30ab_05", "q30ac_05", None, None),
    2007: ("q30aa_07", "q30ab_07", "q30ac_07", "EAq1t06", "EAq1q06"),
    2009: ("q30aa_09", "q30ab_09", "q30ac_09", "EAq1t08", "EAq1q08"),
    2011: ("q30aa_11", "q30ab_11", "q30ac_11", "EAq1t10", "EAq1q10"),
    2013: ("q30aa_13", "q30ab_13", "q30ac_13", "EAq1t12", "EAq1q12"),
    2015: ("q30aa_15", "q30ab_15", "q30ac_15", "EAq1t14", "EAq1q14"),
}
BENCHMARK = 2013

lines = []

def say(s=""):
    print(s)
    lines.append(s)

def load(year, cols):
    path = os.path.join(DATA, f"enterprise {year}.dta")
    if not os.path.exists(path):
        return None
    head = pd.read_stata(path, convert_categoricals=False, iterator=True)
    have = set(head.variable_labels().keys())
    cols = [c for c in cols if c and c in have]
    if not cols:
        return None
    return pd.read_stata(path, convert_categoricals=False, columns=cols)

def knku(df, sn, su, w=None):
    m = df[sn].notna() & df[su].notna()
    if w is not None:
        m &= df[w].notna() & (df[w].astype(float) > 0)
    if m.sum() < 50:
        return m.sum(), np.nan, np.nan
    new, used = df.loc[m, sn].astype(float), df.loc[m, su].astype(float)
    if w is None:
        tn, tu = new.mean(), used.mean()
    else:
        wt = df.loc[m, w].astype(float)
        tn, tu = (new * wt).sum(), (used * wt).sum()
    return int(m.sum()), tn / tu, tu / (tn + tu)

say("=" * 78)
say("NEW-to-USED CAPITAL RATIO (Kn/Ku)  --  Vietnam SME survey")
say("Target for: zeta = (delta+kappa)/KnKu_ratio, an AGGREGATE technology.")
say("Estimator must therefore be aggregate Kn / aggregate Ku, i.e. firm shares")
say("weighted by each firm's equipment stock -- NOT the mean of firm shares.")
say("=" * 78)
say()

say("--- 1. Aggregate Kn/Ku by survey round -------------------------------------")
say(f"{'round':>7} {'N':>7} {'mean_new':>9} {'mean_used':>10} "
    f"{'AGGREGATE':>10} {'agg used sh':>12}")
rows = {}
for yr, (sn, su, ss, wequip, wass) in sorted(ROUNDS.items()):
    df = load(yr, [sn, su, ss, wequip, wass])
    if df is None or sn not in df or su not in df:
        say(f"{yr:>7}   (variables not found)")
        continue
    n_eq, r_eq, _ = knku(df, sn, su, None)
    w = wequip if (wequip and wequip in df) else (wass if (wass and wass in df) else None)
    if w is None:
        say(f"{yr:>7} {n_eq:>7} {df[sn].mean():>9.2f} {df[su].mean():>10.2f} "
            f"{'n/a':>10} {'(no weight var)':>12}")
        continue
    n_ag, r_ag, us_ag = knku(df, sn, su, w)
    mark = "  <== benchmark" if yr == BENCHMARK else ""
    say(f"{yr:>7} {n_ag:>7} {df[sn].mean():>9.2f} {df[su].mean():>10.2f} "
        f"{r_ag:>10.4f} {us_ag:>12.4f}{mark}")
    rows[yr] = dict(n=n_ag, eq=r_eq, agg=r_ag, used=us_ag, weight=w)

say()
say("--- 2. Benchmark year %d: robustness to the weighting variable -------------" % BENCHMARK)
sn, su, ss, wequip, wass = ROUNDS[BENCHMARK]
alt = [(wequip, "Machinery/Equipment, market price   <== PREFERRED"),
       (ROUNDS[BENCHMARK][3].replace("12", "11") if wequip else None,
        "Machinery/Equipment, one year earlier"),
       ("q73ac_13", "Equipment/machinery value (self-reported, mVND)"),
       (wass, "TOTAL physical assets (cruder proxy)")]
df = load(BENCHMARK, [sn, su, ss] + [a for a, _ in alt if a])
say(f"{'weight':>12}  {'N':>6}  {'Kn/Ku':>8}  {'used share':>11}   description")
for a, desc in alt:
    if not a or a not in df:
        continue
    n, r, us = knku(df, sn, su, a)
    say(f"{a:>12}  {n:>6}  {r:>8.4f}  {us:>11.4f}   {desc}")

say()
say("--- 3. HEADLINE Kn/Ku MOMENT -----------------------------------------------")
if BENCHMARK in rows:
    b = rows[BENCHMARK]
    say(f"    benchmark round {BENCHMARK}, weighted by {b['weight']}")
    say(f"      AGGREGATE Kn/Ku                    : {b['agg']:.4f}   <== CALIBRATION TARGET")
    say(f"      aggregate used share of equipment  : {b['used']:.4f}")
    say()
    for kap, dl in [(0.1434, 0.06)]:
        say(f"      implied zeta = (delta+kappa)/KnKu = ({dl}+{kap})/{b['agg']:.4f} "
            f"= {(dl + kap) / b['agg']:.4f}")
    say()
    say(f"    ==> KnKu_ratio = {b['agg']:.4f}")

say("=" * 78)
say("Method note for the paper:")
say("  Firm-level equipment shares (q30aa/q30ab, self-constructed excluded) are")
say("  aggregated using each firm's machinery and equipment stock at market")
say("  prices, so the statistic is the economy-wide new-to-used capital ratio")
say("  rather than the mean across firms. The distinction matters because larger")
say("  firms hold most of the equipment and are more likely to buy used, so the")
say("  unweighted mean of firm shares overstates Kn/Ku.")
say("=" * 78)

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n[saved] {OUT}")
