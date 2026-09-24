# Data moments and the empirical evidence

Everything in the paper that comes from the Vietnam SME survey rather than from
the model: **Tables 1 and 2**, the loan-applicant numbers in Section 2, and the
**targets in Table 3** (the six calibration targets).

---

## Data availability — read this first

The survey microdata is **not included in this package.** The Vietnam SME
survey (enterprise modules, waves 2005–2015) is collected by CIEM and
distributed by UNU-WIDER. Download it from

<https://www.wider.unu.edu/database/viet-nam-sme-database>

Every script below expects the rounds in `enterprise_data/`:

```
enterprise_data/
    enterprise 2005.dta   enterprise 2007.dta   enterprise 2009.dta
    enterprise 2011.dta   enterprise 2013.dta   enterprise 2015.dta
    vsme_panel.dta                       the linked panel built from the rounds
    SME Employee 20YY_clean.dta          employee modules (used by wgini checks)
```

**Nothing needs to be re-run to check the paper.** Every number the paper
reports is already written to a `*_results.txt` file next to the script that
produced it, and those files are committed here. They are the audit trail: a
reader can verify every figure in the paper without holding the data. Running
`compute_data_moments.py` without the data exits with an explanation rather than
a traceback.

---

## The one command

```bash
python compute_data_moments.py
```

Runs every moment subroutine in turn and prints a consolidated scorecard of
**data value vs calibration target**, then writes it to
`data_moments_results.txt`. That file is the single place to check that the
"data" column of paper Table 3 is what the survey actually says.

Current state — every target reproduces its data moment:

| moment | data | target | gap |
|---|---|---|---|
| exit | 0.0978 | 0.0978 | +0.0% |
| DtoY | 0.6019 | 0.6000 | −0.3% |
| KnKu | 2.5699 | 2.5699 | +0.0% |
| top10 | 0.6321 | 0.6321 | +0.0% |
| wGini | 0.3720 | 0.3720 | +0.0% |
| beta_size | −0.0269 | −0.0269 | +0.0% |

Every target is the rounded data value; the only visible gap is DtoY, where
0.6019 is rounded to 0.6000. The same six numbers are hardcoded in
`../paper_targets.m`, which is what the model tables score against, so if a
moment here changes that file must change with it.

---

## What is in each folder

| folder | produces | paper location |
|---|---|---|
| `motivation/` | **Tables 1 and 2**, and the loan-applicant numbers | Section 2 |
| `dtoy/` | debt-to-output = 0.6019 | Table 3 |
| `top10/` | top-decile employment share = 0.6321 | Table 3 |
| `exit/` | firm exit rate = 0.0978 | Table 3 |
| `betasize/` | used-share size gradient = −0.0269 | Table 3 |
| `knku/` | new-to-used capital ratio | sets ζ, not a scored target |
| `wgini/` | wage Gini = 0.372 — **sourced, not computed** | Table 3 |
| `reentry/` | whether failed owners restart | supports the entry-margin discussion |

### motivation/ — paper Tables 1 and 2

Two scripts, run independently of `compute_data_moments.py`:

```bash
cd motivation
python motivation_analysis.py     # -> motivation_results.txt, motivation.doc, fig_used_by_age*.png
python firm_fe_regression.py      # -> main_table.tex, firm_fe_table.tex,
                                  #    collateral_channel_table.tex, firm_fe_results.txt
```

`main_table.tex` **is** Table 2 of the paper, ready to `\input`. The
specification is the used-capital share in percentage points on firm and year
fixed effects, in between and within variants. Two coefficients carry the
section's argument:

- **Zero debt × log(Assets)**: −1.216\*\*\* between, −1.067\*\* within
- **Zero debt × Young**: 2.860\* between, 3.885\* within

The within-firm estimates are the ones that matter. The unconditional age effect
(**Young** alone, +0.755) is **not** significant once firm fixed effects absorb
the cross-sectional composition, so the paper does not claim an unconditional
age gradient; it claims the financing interaction, which survives.

Two traps if these are re-estimated. Birth year is not constant within firm in
this panel — only 31.9% of firms report it consistently across rounds, so firm
age must be built from a single round rather than averaged. And `statsmodels`
does not degrees-of-freedom-correct absorbed fixed effects, so standard errors
from an absorbed fit should be checked against a dummy-variable fit.

### knku/ — the new-to-used capital ratio

`compute_knku_aggregate.py` computes the **aggregate** Kn/Ku: firm-level
equipment shares weighted by each firm's machinery and equipment at market
price. The 2013 benchmark value is **2.5699**, and that is the calibration
target `KnKu_ratio` in every parameter file.

The estimator must be aggregate rather than an average across firms because the
model's ζ is an aggregate technology, fixed by the stationary supply condition
ζ = (δ+κ)·(Ku/Kn). It identifies the economy-wide ratio of new to used capital,
not the ratio faced by the typical firm.

The ratio is not scored as a calibration moment in the usual sense: it pins ζ
directly rather than being a moment the calibration searches to hit.

### wgini/ — sourced, not computed

There is no script here, and that is deliberate. The wage Gini in the paper is
the **economy-wide** Gini among wage earners, which requires a household survey
(VHLSS). This package holds the enterprise modules, whose employee data covers a
different population and would give the wrong object. The target is taken from
the published series in Doan, Ha, Tran and Yang (2023), 2010 value 0.372. See
`wgini/wgini_source.md`.

### reentry/ — not a target

A standalone data check, not part of the scorecard. It establishes that
Vietnamese owners whose firms die do **not** typically restart: about 5% of
founders are re-entrants, 54.6% came from wage work, and 98% of exits are
permanent. This is why the paper treats exit as a move into wage work rather
than as a revolving door.

---

## Requirements

Python 3.9+, with `pandas`, `numpy` and `statsmodels`. Developed on Python
3.11 / pandas 3.0 / statsmodels 0.14.

All paths are resolved relative to each script's own location, so the package
runs from any directory without editing.
