import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "enterprise_data")

V = {
    2005: dict(used="q34b_05", scale=100, estab="q6ab_05"),
    2007: dict(used="q30b_07", scale=1, estab="q6a_07"),
    2009: dict(used="q30ab_09", scale=1, estab="q6a_09"),
    2011: dict(used="q30ab_11", scale=100, estab="q6a_11"),
    2013: dict(used="q30ab_13", scale=100, estab="q6a_13"),
    2015: dict(used="q31ab_15", scale=100, estab="q6a_15"),
}

frames = []
for yr, m in V.items():
    f = os.path.join(DATA, "enterprise %d.dta" % yr)
    d = pd.read_stata(f, convert_categoricals=False,
                      columns=["id", m["used"], m["estab"]])
    d = d[~d["id"].isna()].drop_duplicates(subset="id", keep="first")
    w = pd.DataFrame({"fid": d["id"].astype(int), "wave": yr})
    w["used_share"] = (pd.to_numeric(d[m["used"]], errors="coerce")
                       / m["scale"]).clip(0, 1).values
    est = pd.to_numeric(d[m["estab"]], errors="coerce")
    w["birthyear"] = est.where((est >= 1900) & (est <= yr)).values
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

b = df.dropna(subset=["used_share", "agebin"])
g = b.groupby("agebin", observed=True)["used_share"]
mean = 100 * g.mean()
sd = 100 * g.std()
n = g.count()
se = sd / np.sqrt(n)

ORDER = ["0-5", "6-10", "11+"]
LAB = {"0-5": "Young\n(0-5 years)", "6-10": "6-10 years", "11+": "11+ years"}

print("Used-capital share by firm age, Vietnam SME survey 2005-2015")
print("  %-12s %8s %7s %9s" % ("age group", "mean %", "s.e.", "N"))
for k in ORDER:
    print("  %-12s %8.1f %7.2f %9s"
          % (k, mean[k], se[k], "{:,}".format(int(n[k]))))
print("  total N = {:,}".format(int(n.sum())))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#0d5a9e"
INK = "#1a1a1a"
MUT = "#5a5a5a"
GRID = "#dcdcdc"

plt.rcParams.update({
    "font.size": 8,
    "axes.edgecolor": "#404040",
    "axes.linewidth": 0.8,
})

fig, ax = plt.subplots(figsize=(4.2, 3.0), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

x = np.arange(len(ORDER))
vals = [mean[k] for k in ORDER]
errs = [1.96 * se[k] for k in ORDER]

ax.bar(x, vals, width=0.62, color=BLUE, edgecolor="#25405c", linewidth=0.6,
       zorder=3)
ax.errorbar(x, vals, yerr=errs, fmt="none", ecolor=INK, elinewidth=0.9,
            capsize=3, capthick=0.9, zorder=4)

for xi, v in zip(x, vals):
    ax.text(xi, v + 1.4, "%.1f" % v, ha="center", va="bottom",
            fontsize=8, color=INK, zorder=5)

ax.set_xticks(x)
ax.set_xticklabels([LAB[k] for k in ORDER], color=INK)
ax.set_xlabel("Firm age", color=INK)
ax.set_ylabel("Used capital (% of equipment)", color=INK)
ax.set_ylim(0, max(vals) * 1.25)
ax.yaxis.grid(True, color=GRID, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(colors=MUT, labelcolor=INK, length=0)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_used_by_age.pdf"), format="pdf",
            bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(HERE, "fig_used_by_age.png"), dpi=300,
            bbox_inches="tight", facecolor="white")
print("\n[saved] fig_used_by_age.pdf (vector) + .png (300 dpi)")
