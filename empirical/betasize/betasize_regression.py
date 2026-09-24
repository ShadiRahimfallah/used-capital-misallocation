import pandas as pd, numpy as np, re, os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from targets import TARGETS

HERE = os.path.dirname(os.path.abspath(__file__))
ENT = os.path.join(HERE, "..", "enterprise_data")
PANEL = os.path.join(ENT, "vsme_panel.dta")
OUT_TXT = os.path.join(HERE, "betasize_results.txt")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8"); self.stdout = sys.stdout
    def write(self, s): self.f.write(s); self.stdout.write(s)
    def flush(self): self.f.flush(); self.stdout.flush()

sys.stdout = Tee(OUT_TXT)

def ols_slope(y, x):
    m = np.isfinite(y) & np.isfinite(x)
    y, x = y[m], x[m]
    n = y.size
    if n < 30: return np.nan, np.nan, n
    xd = x - x.mean(); yd = y - y.mean()
    b = (xd*yd).sum()/(xd*xd).sum()
    a = y.mean() - b*x.mean()
    e = y - a - b*x
    se = np.sqrt((e*e).sum()/(n-2)/(xd*xd).sum())
    return b, se, n

for year in (2011, 2013):
    f = os.path.join(ENT, f"enterprise {year}.dta")
    with pd.io.stata.StataReader(f) as rdr:
        lab = rdr.variable_labels()
    df = pd.read_stata(f, convert_categoricals=False)
    print("="*78)
    print(f"FILE: enterprise {year}.dta   ({df.shape[0]} firms, {df.shape[1]} vars)")

    used_cols = [c for c in df.columns if re.match(r"q30a", c, re.I)]
    print("\n  q30a* variables (used-capital question) and their labels:")
    for c in used_cols:
        print(f"    {c:12s} : {lab.get(c)}")

    size_cols = [c for c in df.columns if re.match(r"EAq1q", c, re.I)]
    print("\n  EAq1q* variables (asset/size) and their labels:")
    for c in size_cols[:6]:
        print(f"    {c:12s} : {lab.get(c)}")

    ub = [c for c in used_cols if re.match(r"q30ab", c, re.I)]
    if not ub:
        print("  !! no q30ab found"); continue
    ucol = ub[0]
    us = pd.to_numeric(df[ucol], errors="coerce").to_numpy(float)
    if np.nanmax(us) > 1.5: us = us/100.0

    scol = size_cols[0] if size_cols else None
    if scol is None:
        print("  !! no EAq1q found"); continue
    A = pd.to_numeric(df[scol], errors="coerce").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lnA = np.where(A > 0, np.log(A), np.nan)

    print(f"\n  chosen: used share = {ucol}  |  size = log({scol})")
    print(f"  used-share distribution: N reporting={np.isfinite(us).sum()}, "
          f"share(us>0)={np.nanmean(us>0):.3f}, mean(us|us>0)={np.nanmean(np.where(us>0,us,np.nan)):.3f}")

    b0, se0, n0 = ols_slope(us, lnA)
    b1, se1, n1 = ols_slope(np.where(us > 0, us, np.nan), lnA)

    print(f"\n  beta_size regressions (used share on log assets):")
    print(f"    UNCONDITIONAL : slope={b0:+.4f} (se {se0:.4f}, N={n0})  <- extensive-margin artifact")
    print(f"    INTENSIVE     : slope={b1:+.4f} (se {se1:.4f}, N={n1})  <- the MODEL moment")

    m = np.isfinite(us) & np.isfinite(A) & (A > 0)
    agg_us = (us[m]*A[m]).sum()/A[m].sum()
    print(f"\n  aggregate (asset-weighted) used share = {agg_us:.3f}"
          f"  -> implied Kn/Ku = {(1-agg_us)/max(agg_us,1e-9):.2f}   (calibration target: {TARGETS['KnKu']:.4f})")
print("="*78)

print("\n" + "=" * 78)
print("SIZE GRADIENT OF USED-CAPITAL INTENSITY  (2013 benchmark, panel)")
print("=" * 78)

pan = pd.read_stata(PANEL, columns=["round", "eq_share_used", "workers", "physical_assets"])
d13 = pan[pan["round"] == 2013].copy()
us13 = pd.to_numeric(d13["eq_share_used"], errors="coerce") / 100.0
emp13 = pd.to_numeric(d13["workers"], errors="coerce")
A13 = pd.to_numeric(d13["physical_assets"], errors="coerce")

cls = pd.cut(emp13, [0, 9.5, 49.5, 299.5, np.inf],
             labels=["micro <10", "small 10-49", "medium 50-299", "large 300+"])
print("\n  used-capital share by EMPLOYMENT size class:")
print("    %-15s %7s %12s %12s" % ("class", "N", "used|all", "used|>0"))
for c in cls.cat.categories:
    u = us13[cls == c]
    mu_int = 100 * u[u > 0].mean() if (u > 0).any() else float("nan")
    print("    %-15s %7d %11.1f%% %11.1f%%"
          % (c, int((cls == c).sum()), 100 * u.mean(), mu_int))
print("    large firms are a rounding error: %d firms have >=300 workers "
      "(%d have >=100)." % (int((emp13 >= 300).sum()), int((emp13 >= 100).sum())))

g = pd.DataFrame({"us": us13, "A": A13})
g = g[np.isfinite(g.A) & (g.A > 0) & np.isfinite(g.us)]
g["dec"] = pd.qcut(np.log(g.A), 10, labels=False) + 1
print("\n  used-capital share by ASSET decile (1=smallest ... 10=largest):")
print("    %-7s %6s %10s %10s" % ("decile", "N", "used|all", "used|>0"))
allm = g.groupby("dec")["us"].agg(["mean", "count"])
intm = g[g.us > 0].groupby("dec")["us"].mean()
for dd in range(1, 11):
    print("    %-7d %6d %9.1f%% %9.1f%%"
          % (dd, int(allm.loc[dd, "count"]), 100 * allm.loc[dd, "mean"], 100 * intm.loc[dd]))

print("""
  Reading:
    - INTENSIVE margin (used>0, the model moment): used share DECLINES from the
      smallest to the largest firms (asset decile ~75% -> ~53%; micro ~70% ->
      medium ~47%). This monotone within-SME decline IS beta_size = -0.027, and
      it is spread over ~480 firms/decile -- NOT driven by the 3 large firms.
    - UNCONDITIONAL (zeros included): used share RISES with size, because bigger
      firms are more likely to own ANY used equipment (extensive margin). The
      model has no extensive margin, so the moment conditions on used>0.
    - SME truncation at the top only FLATTENS the gradient; the economy-wide
      slope would if anything be steeper.""")
print("=" * 78)
