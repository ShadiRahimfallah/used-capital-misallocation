import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from targets import TARGETS

import pandas as pd
import pyreadstat

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")
OUT_TXT = os.path.join(HERE, "exit_results.txt")

ROUNDS = [2005, 2007, 2009, 2011, 2013, 2015]
YEARS_SPAN = ROUNDS[-1] - ROUNDS[0]

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
print("FIRM EXIT ('DEATH') RATE  --  Vietnam SME survey, 2005-2015")
print("Method: Berkel, Rand, Tarp & Trifkovic (2020), 'The Vietnam SME Data,")
print("2005-15', Ch. 2 in Rand & Tarp (eds.), OUP/UNU-WIDER, Table 2.3.")
print("=" * 78)

def load_ids(year):
    f = os.path.join(DATA, "enterprise %d.dta" % year)
    df, _ = pyreadstat.read_dta(f, usecols=["id"], encoding="latin1")
    s = pd.to_numeric(df["id"], errors="coerce").dropna().astype("int64")
    return set(s.unique())

print("\n--- 1. Firms interviewed per survey round (persistent id `id`) --------")
present = {}
for r in ROUNDS:
    present[r] = load_ids(r)
    print("    %d : %5d firms" % (r, len(present[r])))
print("    (cf. paper 'Full sample': 2,746 / 2,507 / 2,552 / 2,461 / 2,493 / 2,603)")

print("\n--- 2. Attrition of the 2005 cohort (paper Table 2.3) ----------------")
print("    %-24s %7s %9s %8s %7s"
      % ("panel", "remain", "attrited", "round-%", "cum-%"))

N0 = len(present[ROUNDS[0]])
prev = present[ROUNDS[0]]
label = "2005"
print("    %-24s %7d" % ("2005", N0))

per_round_attr = []
for r in ROUNDS[1:]:
    cur = prev & present[r]
    attrited = len(prev) - len(cur)
    round_rate = attrited / len(prev)
    cum_rate = 1 - len(cur) / N0
    per_round_attr.append(round_rate)
    label += "-%s" % str(r)[2:]
    print("    %-24s %7d %9d %7.1f%% %6.1f%%"
          % (label + " panel", len(cur), attrited, 100 * round_rate, 100 * cum_rate))
    prev = cur

Nf = len(prev)
S10 = Nf / N0
annual_death = 1 - S10 ** (1.0 / YEARS_SPAN)

print("\n--- 3. HEADLINE EXIT MOMENT ------------------------------------------")
print("    survivors 2005 -> 2015 : %5d / %-5d = %.1f%%   (paper 982/2,746 = 35.8%%)"
      % (Nf, N0, 100 * S10))
print("    ten-year attrition     : %.1f%%                    (paper 64%%)"
      % (100 * (1 - S10)))
print("    per-round (2-yr) attr. : %.1f%% - %.1f%%             (paper 15-21%%)"
      % (100 * min(per_round_attr), 100 * max(per_round_attr)))
print("    ==> ANNUAL DEATH RATE  : %.2f%%   ~  0.10          (paper 9-10%%)"
      % (100 * annual_death))
print("\n    Model calibration target:  firm_exit_target = %.5f" % TARGETS["exit"])
print("=" * 78)
print("DONE.  Full log written to exit_results.txt")
