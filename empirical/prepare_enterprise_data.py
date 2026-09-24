import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "SME", "SMEdata")
OUT = os.path.join(HERE, "enterprise_data")

SRC_NAME = "So lieu DN %d_Final_clean.dta"
ID_STEM = {2005: "q1_05", 2007: "q1_07", 2009: "q1_09",
           2011: "q1_11", 2013: "q1_13", 2015: "q1_15"}

os.makedirs(OUT, exist_ok=True)

for yr, stem in sorted(ID_STEM.items()):
    src = os.path.join(SRC, SRC_NAME % yr)
    dst = os.path.join(OUT, "enterprise %d.dta" % yr)
    d = pd.read_stata(src, convert_categoricals=False)
    if stem not in d.columns:
        raise SystemExit("%d: enterprise code %s not in %s" % (yr, stem, src))
    d.insert(0, "id", pd.to_numeric(d[stem], errors="coerce"))
    n_id = int(d["id"].notna().sum())
    n_uniq = int(d["id"].dropna().nunique())
    d.to_stata(dst, write_index=False, version=118)
    print("  %d  %6d rows  %6d with id  %6d unique  -> %s"
          % (yr, len(d), n_id, n_uniq, os.path.basename(dst)))

print("\n[staged] %s" % OUT)
