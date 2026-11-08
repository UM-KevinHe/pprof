# pprof

**pprof** builds quality measures that profile health care providers. It
fits the risk-adjustment model with one effect per provider, turns the
fit into standardized measures, tests every provider against a
reference, calibrates the tests with an empirical null, measures
reliability, selects risk factors, and draws the displays a measure
report needs. One function per analysis step, the same call for every
outcome type.

## Provider profiling

Provider profiling identifies the health care providers whose quality of
care is unusually good or unusually poor, and it is a standard
instrument of accountability in United States health policy. The Centers
for Medicare and Medicaid Services (CMS) profiles every
Medicare-certified dialysis facility, hospital, nursing home, home
health agency and transplant program in the country ([CMS quality
measures](https://www.cms.gov/medicare/quality/measures); Ash et
al. 2012). An adverse profile has consequences of four kinds.

- **Certification.** A provider can lose its Medicare approval: in 2018
  CMS terminated a hospital’s heart transplant program after its
  survival fell well below the national rate
  ([ProPublica](https://www.propublica.org/article/medicare-to-terminate-funding-for-st-lukes-heart-transplant-program-in-houston)).
- **Payment.** Dialysis facilities that miss the standards of the [ESRD
  Quality Incentive
  Program](https://www.cms.gov/medicare/quality/end-stage-renal-disease-esrd-quality-incentive-program)
  lose part of their Medicare payment; for payment year 2027 CMS expects
  to reduce the payment of 3,256 of the 7,582 facilities enrolled in
  Medicare, by about \$21.6 million in total ([Federal Register,
  November
  2025](https://www.govinfo.gov/content/pkg/FR-2025-11-24/pdf/2025-20681.pdf)).
- **Surveillance.** The bottom five percent of dialysis facilities on
  the annual ESRD Outcomes List are inspected in full by state survey
  agencies ([CMS FY2026 Mission and
  Priorities](https://www.cms.gov/files/document/fy-26-mpd.pdf)).
- **Public rating.** Hospital star ratings on [Care
  Compare](https://data.cms.gov/provider-data/topics/hospitals/overall-hospital-quality-star-rating)
  are consulted by patients and referring clinicians.

Given these stakes, the statistical procedure that produces a profile
must be right at each of its three steps. An error at any step
misclassifies a provider: it imposes these consequences on a provider
that has not earned them, or spares one that has.

1.  **Risk adjustment.** Patient outcomes depend not only on the care a
    provider delivers but on the patients themselves: age, comorbidity
    and socioeconomic status vary widely between providers, and a
    provider whose patients carry higher risk appears to perform below
    standard. Without risk adjustment, differences in patient mix are
    mistaken for differences in quality of care (Normand, Glickman and
    Gatsonis 1997). The risk-adjustment model separates the part of the
    outcome attributable to the patients from the part attributable to
    the provider; the provider effects it yields, one per provider, are
    the input to the two steps that follow.
2.  **Standardized measures.** After adjustment, providers with
    different patient populations still cannot be compared on their raw
    outcomes. A standardized measure holds every provider against a
    national norm on a single interpretable scale: the ratio, or the
    difference, between the outcome under the provider’s care and the
    outcome expected under the norm for the provider’s own patients. A
    ratio above one, or a difference above zero, is an outcome higher
    than expected (He 2019).
3.  **Hypothesis testing.** A standardized measure gives the direction
    and size of a provider’s departure from the norm, but a departure
    can arise by chance: a small provider can depart widely through
    random variation alone. A test decides whether the departure is more
    than chance. Only a significant departure is reported as a flag, and
    the flag is what triggers investigation, penalty or public report
    (Jones and Spiegelhalter 2011; Xia et al. 2022).

### Where the steps fail

Each step requires choices, and they are made separately for each
measure by the analysts and clinicians of a measure developer, who know
what a measure is meant to answer but are rarely specialists in the
methods that answer it. The defaults of general-purpose regression
software fail at exactly the places that matter in profiling, and they
fail silently: the code runs and returns coefficients of conventional
appearance.

- **The provider effects are high-dimensional and the providers
  unequal.** A single measure carries an effect for each of thousands of
  providers, each estimated from that provider’s patients alone, so the
  number of parameters grows with the number of providers while the
  information per parameter does not (Neyman and Scott 1948). Small
  providers are estimated with bias, and a provider with no events, or
  only events, produces an estimate that diverges to infinity (Albert
  and Anderson 1984). General-purpose software expands the effects into
  thousands of indicator columns and sweeps them at every iteration,
  which at the scale of a national measure is not acceptable; `pprof`
  updates each provider from its own records, because the provider block
  of the information matrix is diagonal (He et al. 2013; Wu et
  al. 2022), and offers a Firth-corrected fit for the providers whose
  outcomes are all identical (Firth 1993).
- **Dropping or shrinking the provider effects changes the answer.** A
  pooled model without provider effects is tractable but biased whenever
  case mix is related to the provider, which it almost always is
  (Neuhaus and Kalbfleisch 1998). Random effects replace the thousands
  of parameters with one variance, but shrink every estimate toward the
  centre, which attenuates the very extremes profiling seeks, and assume
  that provider effects are independent of case mix (Kalbfleisch and
  Wolfe 2013). Fixed effects remain unbiased under such correlation; the
  package fits fixed effects by default and random and correlated random
  effects for comparison.
- **Which risk factors to carry is itself a decision.** Comorbidities
  enter claims data as thousands of diagnosis codes grouped into
  hundreds of candidate conditions; a model that includes them all
  overfits, one that omits a genuine risk factor penalizes the providers
  that treat it, and CMS guidance asks for a parsimonious model of
  clinically relevant and statistically associated factors. `pprof`
  selects risk factors by lasso and group lasso with the provider
  effects kept unpenalized.
- **Indirect and direct standardization answer different questions.**
  Indirect standardization evaluates each provider on its own patients
  and is the natural way to monitor an individual provider, but the
  results are not comparable between providers; direct standardization
  evaluates every provider on one reference population and is suited to
  ranking, but extrapolates each provider’s model to patients it may
  never treat (He 2019; Varewyck et al. 2014). The same provider can be
  above the norm under one and below it under the other. Both are
  available, and a directly standardized measure is tested on its own
  scale.
- **The classical tests fail at small providers.** The Wald test divides
  the estimate by its standard error, so at a provider whose outcomes
  are all identical both diverge and the statistic cannot flag it; the
  score test avoids this but relies on a large-sample approximation in
  the size of a single provider, which a facility of a dozen patients
  does not satisfy; an exact test removes the approximation (He et
  al. 2013; Wu et al. 2022). All three are provided, with the flags
  reported separately for providers worse and better than expected,
  since only one tail usually carries consequences.
- **Patients are observed repeatedly.** Most measures record each
  patient many times, as months, discharges or periods at risk, so a
  variance that treats the observations as independent is too small; the
  within-patient correlation must be estimated from the provider’s own
  patients, which is unstable at small providers (Mancl and DeRouen
  2001). The package keeps the patient behind every row and offers the
  patient-clustered robust variance and the stabilized variance that
  pools the correlation across patients (Pan 2001).
- **The theoretical null flags large providers.** Referring a statistic
  to N(0, 1) presumes that only sampling variation remains after risk
  adjustment, but the risk factors come from claims generated for
  billing, and residual unmeasured risk stays in the provider effect.
  Large providers are estimated precisely, so the theoretical null flags
  them for differences of no clinical consequence, and the probability
  of a flag rises with size. An empirical null estimates the centre and
  spread of the null from the bulk of the observed statistics, within
  strata of provider size, so that a provider is flagged when it is
  unusual among providers of its size (Efron 2004; Kalbfleisch and Wolfe
  2013; Xia et al. 2022; Hartman and He 2024).

## What the package does

Every measure in `pprof` is built on the same linear predictor. For
patient *j* of provider *i*, observed on occasion *k* (a patient-month,
a discharge, a period at risk),

\\\eta\_{ijk} = \gamma_i + Z\_{ijk}^\top \beta ,\\

with the provider effect γ_(*i*) and the risk factors *Z*_(*ijk*). The
outcome determines the model through which η acts: a logistic model for
a binary or binomial outcome, a linear model for a continuous one, a
proportional-hazards model for a time to event. Once the outcome type,
the unit of observation, the standardization convention and the
reference distribution are fixed, everything else (estimation, variance,
testing, flagging) is common machinery. The package is therefore
organised by **analysis step**, not by quality measure, and a measure,
whatever the outcome it is built on, is assembled from the same
components:

    pp_data()  →  fit_*()  →  standardize()  →  test()  →  null_model() / calibrate()  →  plot()
      data         model       measures          flags       reference distribution       displays

Each choice above is an argument: the statistic, the variance and the
reference distribution of
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md), the
type and norm of
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md),
the grouping and estimator of
[`null_model()`](https://um-kevinhe.github.io/pprof/reference/null_models.md).
Every step is one function with S3 methods for the model classes, so the
same [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md)
call works for a logistic, a linear or a Cox fit, a measure can be
re-specified along one axis with everything else held fixed, and two
measures differ only where they are declared to differ. The
measure-specific conventions (eligibility rules, attribution, reporting)
stay in the analysis script that calls the package.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("UM-KevinHe/pprof", ref = "Yubo")
```

A C++17 compiler is required (Rtools on Windows, Xcode command-line
tools on macOS). OpenMP is used when available.

## A complete analysis in ten lines

Patients nested in providers, a binary outcome per patient-month, a few
risk factors:

``` r

library(pprof)

d  <- sim_kecc(n_providers = 200, seed = 1)           # or your own data frame
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"),     # declare the columns
              provider = "provider", patient = "patient", cutoff = 10)

fit <- fit_logistic_fe(pd)                            # logistic model, one effect per provider
summary(fit)                                          # the risk-factor coefficients

sm  <- standardize(fit, type = "indirect")            # observed / expected per provider
tst <- test(fit, statistic = "score", variance = "robust",   # one test per provider,
            null_model = "empirical")                 # calibrated by an empirical null
tst                                                   # flags: 1 above, -1 below, 0 not different

plot(tst)                                             # caterpillar plot
plot(tst, type = "funnel")                            # funnel plot with the test's limits
plot(standardize(fit), test = tst)                    # observed against expected
```

## Which functions for which outcome

| Outcome | Example measures | Model | Standardized measure | Tests of the provider effects |
|----|----|----|----|----|
| binary (0/1 per row) | readmission within 30 days of a discharge; a fistula in use in a patient-month | [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md) (SerBIN, Firth), [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md), [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md) | indirect ratio or rate; direct rate | exact, score, Wald × naive, robust, stabilized; likelihood ratio |
| binomial (events out of *n* trials per row) | hospital days out of days at risk in a patient-month | the logistic fitters with `n` | as binary | as binary |
| continuous | a laboratory value, a length of stay | [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md), [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md), [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md) | indirect and direct differences | score, Wald × naive, robust |
| time to event (counting-process rows) | death; hospital admissions as recurrent events | [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md) (stratified by provider, two-stage) | observed / expected events, indirect and direct | mid-p Poisson, exact Poisson, log-ratio Wald |
| any of the above, choosing risk factors | measure development | [`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md) (lasso, group lasso), [`cv()`](https://um-kevinhe.github.io/pprof/reference/cv.md) |  |  |

Every test result can be calibrated with an empirical null
([`null_model()`](https://um-kevinhe.github.io/pprof/reference/null_models.md),
[`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md)),
and every model’s
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md)
result can be tested at the measure level with
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md).
Reliability is
[`iur()`](https://um-kevinhe.github.io/pprof/reference/iur.md); the
flags under alternative choices are collected by
[`flag_stability()`](https://um-kevinhe.github.io/pprof/reference/flag_stability.md).

## Documentation

The site <https://um-kevinhe.github.io/pprof/> has a *Get started* page
with the core concepts and one complete analysis, one article per
analysis step (data, risk-adjustment models, variable selection,
standardized measures, tests, reference distributions, reliability,
plots), and the function reference grouped the same way.

## Background

The methods are those of the CMS dialysis quality measures developed at
the University of Michigan Kidney Epidemiology and Cost Center and the
papers behind them. This is version 2.0.0 (development), which
supersedes `pprof` 1.0.3 on CRAN, `pprofR` and `grplasso`; `NEWS.md`
records what changed and why. See `citation("pprof")` and the references
in each help page.

## Getting help

Problems and questions: open an issue at
<https://github.com/UM-KevinHe/pprof/issues> or write to
<ybshao@umich.edu> or <kevinhe@umich.edu>.

## References

Albert A, Anderson JA (1984). On the existence of maximum likelihood
estimates in logistic regression models. *Biometrika*, 71(1), 1-10.

Ash AS, Fienberg SE, Louis TA, Normand S-LT, Stukel TA, Utts J (2012).
*Statistical issues in assessing hospital performance.* Committee of
Presidents of Statistical Societies, report to CMS.

Efron B (2004). Large-scale simultaneous hypothesis testing: the choice
of a null hypothesis. *Journal of the American Statistical Association*,
99(465), 96-104.

Firth D (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80(1), 27-38.

Hartman N, He K (2024). Individualized empirical null estimation for
exact tests of healthcare quality. *Statistics in Medicine*, 43(12),
2403-2420.

He K (2019). Indirect and direct standardization for evaluating
transplant centers. *Journal of Hospital Administration*, 8(1), 9-14.

He K, Kalbfleisch JD, Li Y, Li Y (2013). Evaluating hospital readmission
rates in dialysis facilities; adjusting for hospital effects. *Lifetime
Data Analysis*, 19, 490-512.

He K, Kalbfleisch JD, Yang Y, Li Y (2019). Inter-unit reliability for
nonlinear models. *Statistics in Medicine*, 38(5), 844-854.

He K, Schaubel DE (2014). Methods for comparing center-specific survival
outcomes using direct standardization. *Statistics in Medicine*, 33(12),
2048-2061.

Jones HE, Spiegelhalter DJ (2011). The identification of “unusual”
health-care providers from a hierarchical model. *The American
Statistician*, 65(3), 154-163.

Kalbfleisch JD, Wolfe RA (2013). On monitoring outcomes of medical
providers. *Statistics in Biosciences*, 5, 286-302.

Mancl LA, DeRouen TA (2001). A covariance estimator for GEE with
improved small-sample properties. *Biometrics*, 57(1), 126-134.

Neuhaus JM, Kalbfleisch JD (1998). Between- and within-cluster covariate
effects in the analysis of clustered data. *Biometrics*, 54(2), 638-645.

Neyman J, Scott EL (1948). Consistent estimates based on partially
consistent observations. *Econometrica*, 16(1), 1-32.

Normand S-LT, Glickman ME, Gatsonis CA (1997). Statistical methods for
profiling providers of medical care: issues and applications. *Journal
of the American Statistical Association*, 92(439), 803-814.

Pan W (2001). On the robust variance estimator in generalised estimating
equations. *Biometrika*, 88(3), 901-906.

Spiegelhalter DJ (2005). Funnel plots for comparing institutional
performance. *Statistics in Medicine*, 24(8), 1185-1202.

Varewyck M, Goetghebeur E, Eriksson M, Vansteelandt S (2014). On
shrinkage and model extrapolation in the evaluation of clinical center
performance. *Biostatistics*, 15(4), 651-664.

Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
inference for profiling health care providers. *Statistics in Medicine*,
41(15), 2840-2853.

Xia L, He K, Li Y, Kalbfleisch JD (2022). Accounting for total variation
and robustness in profiling health care providers. *Biostatistics*,
23(1), 257-273.
