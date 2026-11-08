#' Prepare data for a provider model
#'
#' Checks and arranges the data every model and test in the package reads:
#' the outcome, the covariates, the provider identifier and, when the data
#' have repeated observations, the patient identifier and the number of
#' trials. Rows are sorted by provider and patient, provider sizes are
#' counted, and providers below the size `cutoff` or with no or all events
#' are flagged.
#'
#' @param data a data frame.
#' @param y name of the outcome column: 0/1 for a binary outcome, the number
#'   of events out of `n` trials for a binomial outcome, or any numeric value
#'   for a continuous outcome.
#' @param x names of the covariate columns. Covariates must be numeric; code
#'   factors as 0/1 columns first.
#' @param provider name of the provider identifier column.
#' @param patient name of the patient identifier column, or `NULL` when every
#'   row is an independent unit. Required by the patient-clustered tests.
#' @param n name of the column with the number of trials per row (binomial
#'   outcome), or `NULL`.
#' @param outcome `"auto"` (default), `"binary"`, `"binomial"` or
#'   `"continuous"`. With `"auto"` the outcome is binomial when `n` is given,
#'   binary when every value of `y` is 0 or 1, and continuous otherwise.
#' @param cutoff providers with fewer than `cutoff` patients (rows when
#'   `patient` is `NULL`) are kept in the data but marked `included = FALSE`
#'   and left out of model fitting. Default 10.
#' @param check run [pp_check()] first (missing values and zero-variance
#'   covariates are errors; near-zero variance, high correlation and large
#'   VIF are warnings). Default `TRUE`.
#' @param message print a short summary. Default `TRUE`.
#'
#' @details Rows are ordered by `provider`, then `patient`, with
#' [order()]'s radix method, which sorts character identifiers in the C
#' locale and is stable: a patient's rows keep their input order, so sort
#' the data by time within patient first when the stabilized tests will be
#' used. Provider identifiers are converted to character.
#'
#' Only the named columns are kept.
#'
#' @return An object of class `pp_data`: a list with
#' \item{data}{the sorted data frame with the named columns plus the logical
#'   columns `included`, `no_events` and `all_events` (the last two are
#'   `FALSE` for a continuous outcome).}
#' \item{vars}{the column names `y`, `x`, `provider`, `patient`, `n` and the
#'   `outcome` type.}
#' \item{sizes}{one row per provider in the sorted order: `provider`,
#'   `rows`, `patients` (`NA` without a patient column), `included`,
#'   `no_events`, `all_events`.}
#' \item{cutoff}{the cutoff used.}
#'
#' @seealso [pp_check()], [sim_kecc()]
#' @examples
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", cutoff = 10, message = FALSE)
#' pd
#' head(pd$sizes)
#' @export
pp_data <- function(data, y, x, provider, patient = NULL, n = NULL,
                    outcome = c("auto", "binary", "binomial", "continuous"), cutoff = 10,
                    check = TRUE, message = TRUE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  outcome <- match.arg(outcome)
  one_name(y, "y")
  one_name(provider, "provider")
  if (!is.character(x) || length(x) < 1) {
    stop("`x` must be one or more column names.", call. = FALSE)
  }
  if (!is.null(patient)) one_name(patient, "patient")
  if (!is.null(n)) one_name(n, "n")
  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff < 0) {
    stop("`cutoff` must be a single non-negative number.", call. = FALSE)
  }
  cols <- c(y, x, provider, patient, n)
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(cols)) stop("The same column is named twice.", call. = FALSE)

  d <- data[, cols, drop = FALSE]
  if (!is.numeric(d[[y]])) stop("`y` must be numeric.", call. = FALSE)
  not_numeric <- x[!vapply(d[x], is.numeric, logical(1))]
  if (length(not_numeric)) {
    stop("Covariate(s) not numeric: ", paste(not_numeric, collapse = ", "),
         ". Code factors as 0/1 columns first.", call. = FALSE)
  }
  if (outcome == "auto") {
    outcome <- if (!is.null(n)) {
      "binomial"
    } else if (all(d[[y]] %in% c(0, 1))) {
      "binary"
    } else {
      "continuous"
    }
  }
  if (outcome == "binomial" && is.null(n)) {
    stop("A binomial outcome needs the trials column `n`.", call. = FALSE)
  }
  if (outcome != "binomial" && !is.null(n)) {
    stop("`n` is only used with a binomial outcome.", call. = FALSE)
  }
  trials <- if (is.null(n)) rep(1, nrow(d)) else d[[n]]
  if (outcome == "binomial") {
    if (!is.numeric(trials) || any(trials < 1 | trials != round(trials), na.rm = TRUE)) {
      stop("`n` must hold positive whole numbers of trials.", call. = FALSE)
    }
  }
  if (outcome != "continuous") {
    if (any(d[[y]] < 0 | d[[y]] > trials | d[[y]] != round(d[[y]]), na.rm = TRUE)) {
      stop("`y` must hold whole numbers of events between 0 and the number of trials.",
           call. = FALSE)
    }
  }
  if (check) pp_check(d, y = y, x = x, provider = provider, n = n, message = message)

  d[[provider]] <- as.character(d[[provider]])
  if (!is.null(patient)) d[[patient]] <- as.character(d[[patient]])
  ord <- if (is.null(patient)) {
    order(d[[provider]], method = "radix")
  } else {
    order(d[[provider]], d[[patient]], method = "radix")
  }
  d <- d[ord, , drop = FALSE]
  rownames(d) <- NULL
  trials <- trials[ord]

  prov <- d[[provider]]
  levels <- unique(prov)
  prov_f <- factor(prov, levels = levels)
  rows <- as.integer(table(prov_f))
  patients <- if (is.null(patient)) {
    rep(NA_integer_, length(levels))
  } else {
    as.integer(tapply(d[[patient]], prov_f, function(p) length(unique(p))))
  }
  size <- if (is.null(patient)) rows else patients
  included <- size >= cutoff
  if (outcome == "continuous") {
    no_events <- all_events <- rep(FALSE, length(levels))
  } else {
    events <- as.numeric(rowsum(d[[y]], prov_f, reorder = FALSE)[, 1])
    total <- as.numeric(rowsum(trials, prov_f, reorder = FALSE)[, 1])
    no_events <- included & events == 0
    all_events <- included & events == total
  }

  d$included <- rep(included, rows)
  d$no_events <- rep(no_events, rows)
  d$all_events <- rep(all_events, rows)
  sizes <- data.frame(provider = levels, rows = rows, patients = patients,
                      included = included, no_events = no_events, all_events = all_events,
                      stringsAsFactors = FALSE)

  out <- structure(list(data = d,
                        vars = list(y = y, x = x, provider = provider, patient = patient, n = n,
                                    outcome = outcome),
                        sizes = sizes, cutoff = cutoff),
                   class = "pp_data")
  if (message) print(out)
  invisible(out)
}

one_name <- function(value, what) {
  if (!is.character(value) || length(value) != 1 || is.na(value) || !nzchar(value)) {
    stop("`", what, "` must be a single column name.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.pp_data <- function(x, ...) {
  s <- x$sizes
  unit <- if (is.null(x$vars$patient)) "rows" else "patients"
  cat("Prepared data: ", nrow(x$data), " rows, ", nrow(s), " providers",
      if (!is.null(x$vars$patient)) paste0(", ", sum(s$patients), " patients"), "\n", sep = "")
  cat("  outcome ", x$vars$y,
      switch(x$vars$outcome, binary = " (binary)", continuous = " (continuous)",
             binomial = paste0(" out of ", x$vars$n)),
      "; ", length(x$vars$x), " covariate", if (length(x$vars$x) != 1) "s", "\n", sep = "")
  cat("  ", sum(s$included), " providers with at least ", x$cutoff, " ", unit, " included; ",
      sum(!s$included), " excluded\n", sep = "")
  inc <- x$data$included
  if (x$vars$outcome == "continuous") {
    cat("  mean outcome among included rows: ",
        format(mean(x$data[[x$vars$y]][inc]), digits = 4), "\n", sep = "")
  } else {
    cat("  included providers with no events: ", sum(s$no_events),
        "; with all events: ", sum(s$all_events), "\n", sep = "")
    trials <- trials_of(x)
    cat("  event rate among included rows: ",
        format(100 * sum(x$data[[x$vars$y]][inc]) / sum(trials[inc]), digits = 3), "%\n",
        sep = "")
  }
  invisible(x)
}

# Number of trials per row of a pp_data object (1 for a binary outcome).
trials_of <- function(pd) {
  if (is.null(pd$vars$n)) rep(1, nrow(pd$data)) else as.numeric(pd$data[[pd$vars$n]])
}

# The outcome type of a pp_data object or a fit (older objects without the
# field are binary or binomial).
outcome_of <- function(x) {
  v <- x$vars
  if (!is.null(v$outcome)) return(v$outcome)
  if (is.null(v$n)) "binary" else "binomial"
}

# A pp_data object restricted to the rows `keep` (logical or integer), as the
# cross-validation needs: the provider sizes are recounted but the `included`
# flag of each provider is kept.
pp_subset <- function(pd, keep) {
  d <- pd$data[keep, , drop = FALSE]
  rownames(d) <- NULL
  v <- pd$vars
  prov <- d[[v$provider]]
  levels <- unique(prov)
  prov_f <- factor(prov, levels = levels)
  rows <- as.integer(table(prov_f))
  patients <- if (is.null(v$patient)) {
    rep(NA_integer_, length(levels))
  } else {
    as.integer(tapply(d[[v$patient]], prov_f, function(p) length(unique(p))))
  }
  included <- pd$sizes$included[match(levels, pd$sizes$provider)]
  if (outcome_of(pd) == "continuous") {
    no_events <- all_events <- rep(FALSE, length(levels))
  } else {
    trials <- if (is.null(v$n)) rep(1, nrow(d)) else as.numeric(d[[v$n]])
    events <- as.numeric(rowsum(d[[v$y]], prov_f, reorder = FALSE)[, 1])
    total <- as.numeric(rowsum(trials, prov_f, reorder = FALSE)[, 1])
    no_events <- included & events == 0
    all_events <- included & events == total
  }
  d$included <- rep(included, rows)
  d$no_events <- rep(no_events, rows)
  d$all_events <- rep(all_events, rows)
  structure(list(data = d, vars = v,
                 sizes = data.frame(provider = levels, rows = rows, patients = patients,
                                    included = included, no_events = no_events,
                                    all_events = all_events, stringsAsFactors = FALSE),
                 cutoff = pd$cutoff),
            class = "pp_data")
}
