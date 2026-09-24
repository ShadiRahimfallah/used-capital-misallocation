import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from targets import TARGETS

import pandas as pd
import pyreadstat

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
PANEL = os.path.join(DATA, "vsme_panel.dta")
OUT_TXT = os.path.join(HERE, "dtoy_results.txt")

DTOY_TARGET = TARGETS["DtoY"]

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, s):
        self.f.write(s)
        self.stdout.write(s)

    def flush(self):
        self.f.flush()
        self.stdout.flush()

sys.stdout = Tee(OUT_TXT)

print("=" * 78)
print("DEBT-TO-OUTPUT RATIO (DtoY)  --  Vietnam SME survey")
print("Definition: aggregate sum(debt) / sum(value_added), matching the model's")
print("gross-debt aggregate estimator (simulate.m / moments_from_invariant.m).")
print("=" * 78)

df, _ = pyreadstat.read_dta(PANEL, usecols=["year", "round", "debt", "value_added"])
d = pd.to_numeric(df["debt"], errors="coerce")
va = pd.to_numeric(df["value_added"], errors="coerce")

print("\n--- 1. Variables --------------------------------------------------------")
print("    debt        : n=%d, exact zeros=%.1f%% (unlevered firms, kept), "
      "negatives=%d (dropped)" % (d.notna().sum(), 100 * (d == 0).mean(), (d < 0).sum()))
print("    value_added : n=%d, <=0=%d (dropped)" % (va.notna().sum(), (va <= 0).sum()))

def dtoy(mask):
    dd, vv = d[mask], va[mask]
    ok = dd.notna() & vv.notna() & (vv > 0) & (dd >= 0)
    return dd[ok].sum() / vv[ok].sum(), int(ok.sum())

print("\n--- 2. Aggregate DtoY by survey round -----------------------------------")
print("    %-8s %7s %9s" % ("round", "N", "DtoY"))
R = df["round"]
for r in sorted(R.dropna().unique()):
    val, n = dtoy(R == r)
    print("    %-8d %7d %9.4f" % (int(r), n, val))

print("\n--- 3. HEADLINE DtoY MOMENT ---------------------------------------------")
groups = [
    ("2013 round", R == 2013),
    ("2015 round", R == 2015),
    ("2013 + 2015 rounds", R.isin([2013, 2015])),
    ("all rounds pooled", R.notna()),
]
head = None
for name, mask in groups:
    val, n = dtoy(mask)
    if name == "all rounds pooled":
        head = val
    print("    %-20s N=%5d  DtoY = %.4f" % (name, n, val))

print("\n    ==> DATA DtoY (all rounds)  : %.4f" % head)
print("        calibration target      : %.4f" % DTOY_TARGET)
print("        -> the target is the rounded 2013+2015 value; the model value\n           is reported in Table 3, not here.")
print("=" * 78)
print("DONE.  Full log written to dtoy_results.txt")
