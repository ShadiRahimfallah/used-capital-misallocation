import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from targets import TARGETS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_TXT = os.path.join(HERE, "data_moments_results.txt")

ENT = os.path.join(HERE, "enterprise_data")
HAVE_DATA = os.path.isdir(ENT) and bool(os.listdir(ENT))
if not HAVE_DATA:
    print("\n  Running r_target (public World Bank data) before exiting.\n")
    subprocess.run([sys.executable, "r_target.py"],
                   cwd=os.path.join(HERE, "r_target"))
    sys.exit(
        "\n  The Vietnam SME survey microdata is not present.\n"
        "  Expected the enterprise rounds in:\n    %s\n\n"
        "  Download them from\n"
        "    https://www.wider.unu.edu/database/viet-nam-sme-database\n"
        "  See the data-availability note in README.md.\n\n"
        "  The moments computed from them are already recorded, so nothing needs\n"
        "  to be re-run to read the results:\n"
        "    data_moments_results.txt        the consolidated scorecard\n"
        "    <moment>/<moment>_results.txt   the per-moment detail\n"
        "    motivation/used_capital_table.txt   Table 1 of the paper\n"
        "    motivation/fig_used_by_age.pdf      Figure 1 of the paper\n"
        "    motivation/firm_fe_results.txt      the within-firm regressions\n"
        "    r_target/r_target_results.txt       the real interest rate target\n" % ENT
    )

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

MODULES = [
    ("exit",      "exit",     "compute_exit_rate.py"),
    ("DtoY",      "dtoy",     "compute_dtoy.py"),
    ("KnKu",      "knku",     "compute_knku_aggregate.py"),
    ("top10",     "top10",    "compute_top10.py"),
    ("beta_size", "betasize", "betasize_regression.py"),
    ("r",         "r_target", "r_target.py"),
]

EXHIBITS = [
    ("Table 1",             "motivation", "used_capital_table.py"),
    ("Figure 1",            "motivation", "fig_used_by_age.py"),
    ("within-firm regs",    "motivation", "firm_fe_regression.py"),
]

SCORE = [
    ("exit",      0.0978,  TARGETS["exit"],  "Berkel, Rand, Tarp & Trifkovic (2020) -- SME 2005 cohort, annual death rate"),
    ("DtoY",      0.6019,  TARGETS["DtoY"],  "SME accounts -- aggregate sum(debt)/sum(value added), 2013+2015"),
    ("KnKu",      2.5699,  TARGETS["KnKu"],  "SME -- aggregate Kn/Ku, 2013, weighted by machinery at market price"),
    ("top10",     0.6321,  TARGETS["top10"], "SME -- top-10% of firms' share of PAID employment, 2013 raw round"),
    ("wGini",     0.372,   TARGETS["wGini"], "Doan, Ha, Tran & Yang (2023) -- VHLSS wage Gini, 2010 (SOURCED)"),
    ("beta_size", -0.0269, TARGETS["su"],    "SME -- used-share on log assets, intensive margin, 2013"),
    ("r",         0.0310,  TARGETS["r"],     "World Bank WDI -- Viet Nam real deposit rate, 2013"),
]

def banner(txt):
    print("\n" + "#" * 78)
    print("#  " + txt)
    print("#" * 78)

def run_module(name, sub, script):
    banner("MOMENT: %s   (%s/%s)" % (name, sub, script))
    res = subprocess.run(
        [sys.executable, script],
        cwd=os.path.join(HERE, sub),
        capture_output=True, text=True, encoding="utf-8", errors="replace",
    )
    if res.stdout:
        print(res.stdout.rstrip())
    if res.returncode != 0:
        print("  [WARN] %s exited with code %d" % (script, res.returncode))
        if res.stderr:
            print("  stderr:\n" + res.stderr.rstrip())

print("=" * 78)
print("VIETNAM SME CALIBRATION -- DATA MOMENTS")
print("Runs every moment subroutine, then prints the data-vs-target scorecard.")
print("=" * 78)

for name, sub, script in MODULES:
    run_module(name, sub, script)

for name, sub, script in EXHIBITS:
    banner("EXHIBIT: %s   (%s/%s)" % (name, sub, script))
    res = subprocess.run([sys.executable, script], cwd=os.path.join(HERE, sub),
                         capture_output=True, text=True, encoding="utf-8",
                         errors="replace")
    if res.stdout:
        print(res.stdout.rstrip())
    if res.returncode != 0:
        print("  [WARN] %s exited with code %d" % (script, res.returncode))

banner("MOMENT: wGini   (wgini/wgini_source.md)  --  SOURCED, not computed")
print("  Wage Gini is NOT computed from the SME survey. The model has occupational")
print("  choice, so wGini is the ECONOMY-WIDE wage-earner Gini (a household survey,")
print("  VHLSS), which we do not hold. We cite the published moment:")
print("    Doan, Ha, Tran & Yang (2023), 'Dynamics of wage inequality over the")
print("    prolonged economic transformation: The case of Vietnam,' Economic")
print("    Analysis and Policy 78, 816-834.  Wage Gini: 0.353 (1998) -> 0.372 (2010)")
print("    -> 0.285 (2020).  Target = the 2010 value 0.372 ~ 0.37.")
print("  See wgini/wgini_source.md.")

banner("CONSOLIDATED SCORECARD  (data vs calibration target)")
print("  %-10s %9s %9s %8s   %s"
      % ("moment", "data", "target", "gap", "source"))
print("  " + "-" * 104)
for name, data, target, source in SCORE:
    gap = 100.0 * (target - data) / abs(data)
    print("  %-10s %9.4f %9.4f %+7.1f%%   %s"
          % (name, data, target, gap, source))
print("""
  Notes:
    - gap = (target - data) / |data|: how far the calibration target sits from
      the data moment computed above.  Targets are the rounded data values.
    - MODEL values are deliberately NOT reproduced here.  They live in
      paper_tables.m, Table 3, which reads them from the saved equilibria.
      Quoting them in this package would let the data appendix drift out of
      step with the model tables.
    - KnKu is the aggregate new-to-used capital ratio, firm shares weighted by
      each firm's machinery and equipment at market price.  That is the
      estimator the model computes, and it is what pins zeta.
    - top5 and ent share are omitted: top5 is not a target; ent share is a model
      OUTCOME, not a data moment.""")
print("=" * 78)
print("DONE.  Combined log in data_moments_results.txt;")
print("       per-moment detail in each <folder>/<name>_results.txt.")
