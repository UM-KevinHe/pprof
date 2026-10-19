#' Data preparation for the fixed provider effects model
#'
#' Combines outcome, number of trials, provider identifier, optional patient
#' identifier and covariates into one data frame sorted by provider (and
#' patient), flags small providers and providers with no or all events.
#'
#' @param N number of trials for each row (binomial outcome), a vector or
#'   one-column data frame. If missing, every row is a single Bernoulli trial
#'   (`N = 1`), the binary case.
#' @param Y outcome: number of events for each row (0/1 when `N = 1`), a
#'   one-column data frame or a vector.
#' @param Z a numeric matrix or data frame of covariates.
#' @param ProvID provider identifier, a one-column data frame or a vector.
#' @param PatientID optional patient identifier. Required for the
#'   patient-clustered robust Wald tests.
#' @param cutoff providers with fewer than `cutoff` patients (or rows when
#'   `PatientID` is missing) are labelled `included = 0` and are left out of
#'   model fitting. Default 10.
#' @param check a logical; run [data_check()]-style checks (missingness,
#'   variation, correlation, VIF). Default `TRUE`.
#' @param message a logical; print progress information. Default `TRUE`.
#'
#' @details The returned data frame keeps all providers; small ones are marked
#' `included = 0`, and [logis_fe()] only uses rows with `included == 1`.
#' Columns `no.events` and `all.events` flag included providers whose
#' outcomes are all 0 or all equal to the number of trials.
#'
#' @return An object of class `data_prep`: a list with
#' \item{data}{the sorted data frame.}
#' \item{char_list}{column names: `N.char`, `Y.char`, `prov.char`, `Z.char`
#'   and, when supplied, `PatientID.char`.}
#'
#' @export
fe_data_prep <- function(N, Y, Z, ProvID, PatientID, cutoff = 10, check = TRUE, message = TRUE) {
  if (missing(N)) {
    N <- rep(1, NROW(Y))
    if (message == TRUE) message("No number of trials provided. Defaulting to 1.")
  }

  # check dimensions of the input data
  if (NROW(N) != NROW(Y) | NROW(Y) != NROW(ProvID) | NROW(ProvID) != NROW(Z)) {
    stop("Dimensions of the input data do not match!!", call. = FALSE)
  }

  N <- data.frame(N)
  if (missing(PatientID)) {
    data <- as.data.frame(cbind(N, Y, ProvID, Z))
  } else {
    if (NROW(PatientID) != NROW(Y)) stop("Dimensions of the input data do not match!!", call. = FALSE)
    data <- as.data.frame(cbind(N, Y, ProvID, PatientID, Z))
    PatientID.char <- colnames(data)[4]
  }
  N.char <- colnames(data)[1]
  Y.char <- colnames(data)[2]
  prov.char <- colnames(data)[3]
  Z.char <- colnames(Z)
  if (is.null(Z.char)) stop("Covariates 'Z' must have column names!", call. = FALSE)

  if (check) {
    check_vars(data, N.char, Y.char, Z.char, prov.char, message = message)
  }

  prov.size.r <- as.integer(table(data[, prov.char])) # provider sizes (based on # records)
  if (missing(PatientID)) {
    data <- data[order(factor(data[, prov.char])), ] # sort data by provider ID
    prov.size.p <- prov.size.r
  } else {
    data <- data[order(factor(data[, prov.char]), factor(data[, PatientID.char])), ] # sort by provider & patient
    prov.size.p <- as.integer(table(unique(data[, c(prov.char, PatientID.char)])[, 1])) # provider sizes (based on # patients)
  }
  prov.size.long <- rep(prov.size.p, prov.size.r) # provider sizes assigned to records

  data$included <- 1 * (prov.size.long >= cutoff) # indicator of providers kept for model fitting
  if (message == TRUE) warning(sum(prov.size.p < cutoff), " out of ", length(prov.size.p),
                               " providers considered small and filtered out!", immediate. = TRUE, call. = FALSE)

  inc <- data$included == 1
  prov.list <- unique(data[inc, prov.char]) # a reduced list of provider IDs
  prov.inc <- factor(data[inc, prov.char])
  prov.no.events <- # providers with no events
    prov.list[sapply(split(data[inc, Y.char], prov.inc), sum) == 0]
  data$no.events <- 0
  data$no.events[data[, prov.char] %in% c(prov.no.events)] <- 1
  if (message == TRUE) message(paste(length(prov.no.events), "out of", length(prov.list),
                                     "remaining providers with no events."))
  prov.all.events <- # providers with all events (Y = N for every row)
    prov.list[sapply(split(data[inc, N.char] - data[inc, Y.char], prov.inc), sum) == 0]
  data$all.events <- 0
  data$all.events[data[, prov.char] %in% c(prov.all.events)] <- 1
  if (message == TRUE) message(paste(length(prov.all.events), "out of", length(prov.list),
                                     "remaining providers with all events."))
  if (message == TRUE) message(paste0("After screening, ",
                                      round(sum(data[inc, Y.char]) / sum(data[inc, N.char]) * 100, 2),
                                      "% of all trials exhibit occurrences of events"))

  char_list <- list(N.char = N.char,
                    Y.char = Y.char,
                    prov.char = prov.char,
                    Z.char = Z.char)
  if (!missing(PatientID)) {
    char_list$PatientID.char <- PatientID.char
  }

  structure(list(data = data, char_list = char_list), class = "data_prep")
}
