# pprofR 0.2.0 (2026-09-24)

* New provider tests in `test_fe()`: `"stabilized_wald"`, `"robust_score"`
  and `"stabilized_score"` (patient-clustered; need `PatientID`). They follow
  Feiyang's `robust_tests.R`, extended to binomial outcomes, and reproduce it
  (binary) and the DASH `aoh_tests_en.R` (binomial) to 1e-15.
* Direct standardization (`SR_output(stdz = "direct")`, `direct_exp()`,
  `confint()`) merges rows sharing a linear predictor, avoids `exp()` in its
  inner loop and sums in long double: about 15 times faster with continuous
  covariates (199 s to 13 s on 2.6 million rows and 3000 providers), far
  more with categorical ones, and more accurate. Results change only by
  rounding (at most 3e-12 relative in the checks). For a rate within about
  1e-7 of 0 or 1, the former sums could be off by about 1e-7 relative in
  1 - rate; the long double sums reach double precision.
* Removed the standard score test (`test_fe(score_modified = FALSE)`), which
  no measure used. `test = "score"` (beta held fixed) is unchanged; passing
  `score_modified = FALSE` is now an error.
* The patient-level sums of the robust Wald test run in C++ in the same order
  and precision as before: results are unchanged bit for bit.

# pprofR 0.1.0 (2026-09-23)

* First release: shared engine for the AOH, PPPW, SFR and SRR measures.
