import json
import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "worldbank_raw.json")
OUT_TXT = os.path.join(HERE, "r_target_results.txt")

SERIES = {
    "FR.INR.DPST":       "Deposit interest rate (%)",
    "NY.GDP.DEFL.KD.ZG": "Inflation, GDP deflator (annual %)",
    "FR.INR.RINR":       "Real interest rate, lending (%)",
    "FR.INR.LEND":       "Lending interest rate (%)",
}
YEARS = [str(y) for y in range(2005, 2017)]
BENCHMARK_YEAR = "2013"
MEAN_WINDOW = ["2012", "2013", "2014", "2015", "2016"]
DEFLATOR_MAX = 20.0


def fetch():
    out = {}
    try:
        for ind in SERIES:
            url = ("https://api.worldbank.org/v2/country/VNM/indicator/%s"
                   "?format=json&date=2005:2016&per_page=100" % ind)
            rows = json.load(urllib.request.urlopen(url, timeout=30))[1]
            out[ind] = {r["date"]: r["value"] for r in rows}
        with open(CACHE, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
        return out, "World Bank API (cached to worldbank_raw.json)"
    except Exception as exc:
        if os.path.exists(CACHE):
            with open(CACHE, encoding="utf-8") as fh:
                return json.load(fh), "worldbank_raw.json (API unreachable: %s)" % type(exc).__name__
        raise SystemExit("No API access and no cache at %s" % CACHE)


D, SOURCE = fetch()
dep = D["FR.INR.DPST"]
defl = D["NY.GDP.DEFL.KD.ZG"]
lend_r = D["FR.INR.RINR"]
lend_n = D["FR.INR.LEND"]


def real_dep(y):
    if dep.get(y) is None or defl.get(y) is None:
        return None
    return dep[y] - defl[y]


O = []
O.append("=" * 78)
O.append("  PROVENANCE OF THE REAL INTEREST RATE TARGET")
O.append("  R_TARGET = 0.031 in benchmark/parameters_benchmark.m")
O.append("=" * 78)
O.append("")
O.append("  Source: %s" % SOURCE)
O.append("  Series: " + ", ".join("%s = %s" % (k, v) for k, v in SERIES.items()))
O.append("")
O.append("  real deposit rate = FR.INR.DPST - NY.GDP.DEFL.KD.ZG")
O.append("")
O.append("  %-6s %11s %11s %13s %13s %10s" %
         ("year", "deposit", "deflator", "REAL deposit", "real lending", "usable"))
O.append("  " + "-" * 72)
for y in YEARS:
    rd = real_dep(y)
    ok = "" if defl.get(y) is None else ("yes" if defl[y] <= DEFLATOR_MAX else "NO")
    O.append("  %-6s %11s %11s %13s %13s %10s"
             % (y,
                "--" if dep.get(y) is None else "%.3f" % dep[y],
                "--" if defl.get(y) is None else "%.3f" % defl[y],
                "--" if rd is None else "%.3f" % rd,
                "--" if lend_r.get(y) is None else "%.3f" % lend_r[y],
                ok))
O.append("")
O.append("  'usable' is no where GDP-deflator inflation exceeds %.0f percent, at"
         % DEFLATOR_MAX)
O.append("  which point the real rate is dominated by the deflator and turns")
O.append("  negative; 2011 is the clear case at %.1f percent inflation." % defl["2011"])
O.append("")

rb = real_dep(BENCHMARK_YEAR)
O.append("-" * 78)
O.append("  THE TARGET")
O.append("-" * 78)
O.append("  %s: %.3f - %.3f = %.3f percent   ->   R_TARGET = %.3f"
         % (BENCHMARK_YEAR, dep[BENCHMARK_YEAR], defl[BENCHMARK_YEAR], rb,
            round(rb, 1) / 100))
O.append("")
win = [real_dep(y) for y in MEAN_WINDOW if real_dep(y) is not None]
O.append("  Check on the choice of year. The %s-%s mean real deposit rate is"
         % (MEAN_WINDOW[0], MEAN_WINDOW[-1]))
O.append("  %.2f percent, against %.2f in %s, so the target does not depend on"
         % (sum(win)/len(win), rb, BENCHMARK_YEAR))
O.append("  taking a single year. %s is used because the used-share gradient, the"
         % BENCHMARK_YEAR)
O.append("  top-decile employment share and the new-to-used ratio are all %s SME"
         % BENCHMARK_YEAR)
O.append("  survey moments.")
O.append("")

lw = [lend_r[y] for y in MEAN_WINDOW if lend_r.get(y) is not None]
O.append("-" * 78)
O.append("  WHY THE DEPOSIT RATE AND NOT THE LENDING RATE")
O.append("-" * 78)
O.append("  real deposit  %s: %.2f percent    %s-%s mean %.2f"
         % (BENCHMARK_YEAR, rb, MEAN_WINDOW[0], MEAN_WINDOW[-1], sum(win)/len(win)))
O.append("  real lending  %s: %.2f percent    %s-%s mean %.2f"
         % (BENCHMARK_YEAR, lend_r[BENCHMARK_YEAR], MEAN_WINDOW[0],
            MEAN_WINDOW[-1], sum(lw)/len(lw)))
O.append("  gap           %s: %.2f points" % (BENCHMARK_YEAR, lend_r[BENCHMARK_YEAR] - rb))
O.append("")
O.append("  Both are correct, and the gap IS Vietnamese financial")
O.append("  underdevelopment: the return to capital sits above the US benchmark")
O.append("  of roughly 4.5 percent because capital is scarce, while the household")
O.append("  savings return sits below it because of repression and a wide")
O.append("  intermediation spread.")
O.append("")
O.append("  The model has ONE interest rate. Savers earn r, borrowers pay r, and")
O.append("  the friction is a quantity constraint -- the collateral limit theta --")
O.append("  not a price wedge. beta is identified by the household Euler equation,")
O.append("  so it is pinned by what savers earn, and the deposit rate is the right")
O.append("  counterpart. Do NOT average the two: the average corresponds to")
O.append("  nothing in the model. Report robustness over r in [0.03, 0.06]")
O.append("  instead, which spans the savings return to the borrowing cost.")
O.append("")
O.append("  METHOD vs NUMBER. Calibrating beta to match the real interest rate,")
O.append("  rather than assigning beta, follows Tan (2025). The number is the")
O.append("  Vietnamese series above.")
O.append("=" * 78)

txt = "\n".join(O)
print(txt)
with open(OUT_TXT, "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
print("\n[saved] %s" % OUT_TXT)
