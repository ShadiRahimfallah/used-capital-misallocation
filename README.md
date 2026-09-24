# Replication package

**Financial Frictions, Used Capital, and Misallocation**
Shadi Rahimfallah, Department of Economics, York University
shadirf@yorku.ca

This package contains the code that produces every table and figure in the paper,
together with the stored model solutions the figure and table scripts read.

## Data

The empirical moments are computed from the **Viet Nam SME Survey** (2005, 2007,
2009, 2011, 2013, 2015), collected by CIEM, ILSSA, the Development Economics
Research Group at the University of Copenhagen, and UNU-WIDER.

The survey files are **not included here** — they are distributed by UNU-WIDER and
may not be redistributed. They are supplied free of charge at

> https://www.wider.unu.edu/database/viet-nam-sme-database

on completion of UNU-WIDER's data request form. To run the empirical scripts, place
the enterprise files (`enterprise 2005.dta` … `enterprise 2015.dta`) in
`empirical/enterprise_data/`. The survey instruments, which document every variable
used, are openly downloadable from the same page.

## What is here

```
benchmark/            model with new and used capital: solution, calibration,
                      market clearing, invariant distribution, diagnostics
model_with_no_used/   the single-vintage economy, recalibrated to the same targets
empirical/            the survey moments: used-capital shares, firm fixed-effect
                      regressions, exit rate, debt-to-output, top-decile employment
                      share, wage Gini, new-to-used ratio, size gradient
figures.m             draws Figure 1 and Figure 2 from stored results
tables.m              prints every table in the paper
make_readme.py        regenerates README.pdf
FIGURES.txt           transcript: the numbers behind every plotted curve
TABLES.txt            transcript: the numbers behind every table
```

The two Appendix C comparison economies — the tail-pinned comparison
(`cf_calib_B/`) and the tail diagnostic (`cf_etap/`) — are not included in this
package. They are available from the author on request. Without them,
`tables.m` still produces every other table; only the Appendix C rows for these
two economies are left blank.

## Where each table and figure comes from

| Paper | Script | Output |
|---|---|---|
| Table 1 (used share by firm age) | `empirical/motivation/used_capital_table.py` | `used_capital_table.txt` |
| Table 2 (financing frictions and used capital) | `empirical/motivation/firm_fe_regression.py` | `main_table.tex`, `firm_fe_results.txt` |
| Section 2 loan-applicant and non-applicant numbers | `empirical/motivation/credit_constraints.py` | `credit_constraints_results.txt` |
| Table 3 targets | `empirical/compute_data_moments.py`; interest rate in `empirical/r_target/`; wage Gini sourced in `empirical/wgini/` | `data_moments_results.txt`, `paper_targets.m` |
| Tables 3–7, model columns | `tables.m` | `TABLES.txt` |
| Figures 1–2 | `figures.m` | `fig_mechanism.pdf`, `fig_mrpk_by_wealth.pdf`, `FIGURES.txt` |

`TABLES.txt` echoes the survey tables from `used_capital_table.txt`. The table it
labels "Table 2" (used share by employment size) is a memo table that is not in
the paper; the paper's Table 2 is `empirical/motivation/main_table.tex`.

Requirements: MATLAB for the model (`.m`), Python 3 with pandas and statsmodels for
the empirical scripts (`.py`). The `.mat` files hold solved equilibria, so the
figures and tables can be reproduced in about a minute without re-solving anything.
The stored results were produced on 64-bit Windows. Re-solving the benchmark from
scratch (`benchmark/verify_benchmark.m`) takes about forty minutes.

## Order to run

1. `empirical/` — with the survey files in place, these produce the calibration
   targets and the regressions in Section 2 of the paper. See `empirical/README.md`.
2. MATLAB, from the package root: `>> figures` writes `fig_mechanism.pdf` (Figure 1,
   both panels), `fig_mrpk_by_wealth.pdf` (Figure 2), `fig_usedcapital_tilt.pdf`
   (Figure 1's right panel on its own), and the transcript `FIGURES.txt`.
3. MATLAB: `>> tables` prints every table in the paper and writes `TABLES.txt`.

Re-solving the model from scratch, rather than reading the stored `.mat` files, is
driven from `benchmark/calibrate_newton.m` and, for the economy without used
capital, `model_with_no_used/calibrate_cf.m`. `README.pdf` documents the model
code in detail.

## Code sources and acknowledgements

The solution method follows the computational approaches of

- Midrigan, V. and D. Y. Xu (2014), "Finance and Misallocation: Evidence from
  Plant-Level Data," *American Economic Review* 104(2), 422–458; and
- Tan, E., "Entrepreneurial Investment Dynamics and the Wealth Distribution,"
  *American Economic Journal: Macroeconomics*, forthcoming.

The code in this package was written for this paper, drawing on the methods
of these papers.

Function approximation uses Miranda and Fackler's CompEcon toolbox
(Miranda, M. J. and P. L. Fackler, 2002, *Applied Computational Economics and
Finance*, MIT Press). The toolbox is not redistributed here; download it
separately and place the `CompEcon` folder next to the model folders
(see `startup.m`).

## Licence

MIT.
