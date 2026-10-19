# Small simulated data sets shared by the tests.
sim_fe <- function(m = 25, n_range = c(15, 40), binomial = FALSE, seed = 1) {
  set.seed(seed)
  n_pat <- sample(n_range[1]:n_range[2], m, replace = TRUE)
  prov <- sprintf("P%03d", rep(seq_len(m), n_pat))
  pat <- sprintf("K%05d", seq_along(prov))
  rec <- sample(1:3, length(prov), replace = TRUE)
  prov <- rep(prov, rec)
  pat <- rep(pat, rec)
  n <- length(prov)
  Z <- data.frame(x1 = round(rnorm(n), 3), x2 = rbinom(n, 1, 0.4), x3 = round(runif(n, -1, 1), 3))
  gam <- stats::setNames(rnorm(m, -0.3, 0.5), sprintf("P%03d", seq_len(m)))
  eta <- gam[prov] + 0.5 * Z$x1 - 0.4 * Z$x2 + 0.3 * Z$x3
  N <- if (binomial) sample(5:20, n, replace = TRUE) else rep(1, n)
  Y <- stats::rbinom(n, N, stats::plogis(eta))
  list(N = N, Y = data.frame(Y = Y), Z = Z, ProvID = data.frame(prov = prov),
       PatientID = data.frame(pat = pat))
}

prep_fe <- function(d, binomial = FALSE) {
  suppressWarnings(suppressMessages(
    if (binomial) {
      fe_data_prep(N = d$N, Y = d$Y, Z = d$Z, ProvID = d$ProvID, PatientID = d$PatientID,
                   check = FALSE, message = FALSE)
    } else {
      fe_data_prep(Y = d$Y, Z = d$Z, ProvID = d$ProvID, PatientID = d$PatientID,
                   check = FALSE, message = FALSE)
    }))
}

# Reference fit: glm with one dummy per provider.
glm_fe <- function(d) {
  df <- data.frame(Y = d$Y$Y, N = d$N, prov = factor(d$ProvID$prov), d$Z)
  stats::glm(cbind(Y, N - Y) ~ 0 + prov + x1 + x2 + x3, family = stats::binomial(), data = df,
             control = stats::glm.control(epsilon = 1e-14, maxit = 100))
}
