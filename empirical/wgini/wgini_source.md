# Wage Gini moment — source

**Target: wage Gini = 0.37.** This moment is **not computed from the Vietnam SME
survey** — it is taken from a published study.

**We use this paper:**
Doan, T., V. Ha, T. Tran, and J. Yang (2023), "Dynamics of wage inequality over
the prolonged economic transformation: The case of Vietnam," *Economic Analysis
and Policy* 78, 816–834. doi:10.1016/j.eap.2023.04.014. (PDF: `teaching/vietnam.pdf`.)

They compute the Gini of **hourly wages among wage earners** from the Vietnam
Household Living Standards Survey (VHLSS/GSO). Reported series (p. 828):
the wage Gini rose from **0.353 (1998)** to **0.372 (2010)**, then fell to
**0.285 (2020)**. We use the **2010 value, 0.372 ≈ 0.37** — the peak and the
wave closest to our 2013 calibration benchmark.

**Why sourced, not computed:** the model has occupational choice, so the wage
Gini is the Gini of worker wages **economy-wide** — the right counterpart is the
national wage-earner Gini from a household survey (VHLSS), which we do not hold.
The SME survey's employee module is only a small manufacturing-worker sample
(wrong population), so we cite the published moment instead.

**Model vs target:** model wGini ≈ 0.3486 (−6% vs 0.37).
