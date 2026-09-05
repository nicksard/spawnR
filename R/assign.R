#' Assign parentage with a stated confidence criterion
#'
#' Runs a complete analysis: LOD scores for every offspring against every
#' candidate, followed by one of two significance criteria.
#'
#' `criterion = "pairwise"` uses the trio-specific rule of Amiri Roudbar et al.
#' (2025). For each offspring, [lod_null()] runs a backward simulation
#' conditioned on that offspring's own genotypes and missing-data pattern, and
#' the critical value is the lower-tail quantile of the resulting true-parent LOD
#' distribution. The top candidate is assigned when its LOD reaches that
#' threshold. Because a separate simulation is run per offspring, this is the
#' more expensive option, and it is the one that adapts to unusual genotypes and
#' to uneven missing data.
#'
#' @section Choosing a criterion:
#' `criterion = "posterior"` is the default and is usually the right answer. It
#' reports, for each offspring, the probability that the named candidate is the
#' parent, given the LOD scores, the size of the pool searched and
#' `prop_sampled`. Because it is a statement about one assignment rather than a
#' decision rule, it does not have to trade recall against precision the way the
#' two threshold rules below do: in simulation across pool sizes from 50 to 5000
#' and sampling fractions from 0.5 to 1, precision at a posterior of 0.95 stayed
#' between 98.5% and 100% in every cell, and the stated probabilities tracked
#' observed correctness throughout. See [parentage_posterior()] and
#' [parentage_fdr()].
#'

#' The two criteria control different things, and the difference matters as the
#' candidate pool grows. `"pairwise"` calibrates its threshold on the
#' distribution of the LOD score for a *true* parent, so it controls the
#' probability of rejecting a true parent - a candidate clears the bar when its
#' score is one a true parent would plausibly produce. It says nothing about how
#' many unrelated candidates also clear it, and with a large pool some will, by
#' chance alone. In simulations with fifteen microsatellite loci and the true
#' sire always present, the proportion of assignments that were wrong ran at 0%
#' with 100 or 500 candidates, 3.1% with 2000 and 5.6% with 5000, all at a
#' nominal 99% level.
#'
#' `"delta"` defines confidence the other way round, as the proportion of
#' assignments above a threshold that are correct, and its simulation includes
#' the whole candidate pool and the chance that the true parent is missing from
#' it. That is the quantity most parentage studies actually want to report, and
#' it is the one that stays interpretable as the pool grows.
#'
#' Use `"pairwise"` when the question is whether a specific named candidate is
#' the parent, when missing data are uneven across individuals, or when the pool
#' is small. Use `"delta"` when screening a large pool, and in general when
#' reporting an assignment rate. Running both and comparing is cheap.
#'
#' `criterion = "delta"` uses the population-level rule of Marshall et al.
#' (1998). One call to [sim_delta()] produces a single critical Delta for the
#' whole dataset, and the top candidate is assigned when the gap between the best
#' and second-best LOD reaches it. This is what CERVUS reports, and it is much
#' cheaper, but it applies an average threshold to every case.
#'
#' @param g A `"genotypes"` object.
#' @param offspring Identifiers of the offspring to assign.
#' @param candidates Identifiers of the candidate parents.
#' @param freqs An `"allele_freqs"` object.
#' @param error Genotyping error rate, one value or one per locus.
#' @param known For `type = "one_known"`, the known second parent of each
#'   offspring: one identifier, or one per offspring, `NA` allowed.
#' @param type `"both_unknown"` or `"one_known"`.
#' @param criterion `"pairwise"` or `"delta"`.
#' @param confidence Confidence levels to report. For `"pairwise"` these become
#'   significance levels `1 - confidence`; for `"delta"` they are passed to
#'   [delta_critical()]. The defaults are the two levels each source uses.
#' @param nsim Simulation replicates, per offspring for `"pairwise"` and in
#'   total for `"delta"`.
#' @param prop_sampled Prior probability that the true parent is among the
#'   candidates offered. Required for `criterion = "posterior"`, where it is a
#'   real modelling assumption that should be varied rather than guessed once;
#'   defaults to 1 for `criterion = "delta"`. A value of exactly 1 forbids the
#'   model from concluding that no candidate is the parent, so it will always
#'   return a best guess.
#' @param prop_typed For `criterion = "delta"`, the assumed proportion of
#'   genotypes scored. Defaults to the observed proportion.
#' @param max_mismatch Optional cap on mismatching loci for a candidate to
#'   remain eligible.
#'
#' @return An object of class `"parentage_assignment"`: a `data.frame` with one
#'   row per offspring giving the top candidate, its LOD, `delta`, the critical
#'   values at each confidence level, and a `confidence` column holding the
#'   highest level attained (`NA` where none is). The full LOD table and the
#'   criterion details are attached as attributes `lod_table` and `criteria`.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x}
#'
#' Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
#' (2025) Pairwise paternity assignment with forward-backward simulations:
#' refining CERVUS using trio-based likelihood and locus-specific error rates.
#' \emph{Ecology and Evolution}, \strong{15}(10), e72230.
#' \doi{10.1002/ece3.72230}
#'
#' @examples
#' sim <- simulate_population(n = 40, n_loci = 10, n_alleles = 8, seed = 9)
#' p <- allele_freqs(sim$genotypes)
#' a <- assign_parentage(sim$genotypes, sim$pedigree$offspring[1:5],
#'                       sim$sires, p, error = 0.01,
#'                       criterion = "pairwise", nsim = 400)
#' a
#' @export
assign_parentage <- function(g, offspring, candidates, freqs, error,
                             known = NULL,
                             type = c("both_unknown", "one_known"),
                             criterion = c("posterior", "pairwise", "delta"),
                             confidence = NULL, nsim = 10000L,
                             prop_sampled = NULL, prop_typed = NULL,
                             max_mismatch = NULL) {
  stopifnot(inherits(g, "genotypes"))
  type <- if (!is.null(known)) "one_known" else match.arg(type)
  criterion <- match.arg(criterion)
  if (is.null(confidence)) {
    confidence <- if (criterion == "delta") c(0.80, 0.95) else c(0.95, 0.99)
  }
  confidence <- sort(confidence)
  if (criterion == "posterior" && is.null(prop_sampled)) {
    stop("`prop_sampled` is required for criterion = 'posterior': it is the ",
         "prior probability that the true parent is among the candidates ",
         "offered. Give your best estimate, and vary it to check sensitivity.",
         call. = FALSE)
  }
  if (is.null(prop_sampled)) prop_sampled <- 1

  tab <- parentage_lod(g, offspring, candidates, freqs, error, known = known,
                       type = type, max_mismatch = max_mismatch)
  best <- delta_stat(tab)
  off <- .as_ids(g, offspring, "offspring")
  kn <- if (type == "one_known") {
    k <- as.character(known); if (length(k) == 1L) rep(k, length(off)) else k
  } else rep(NA_character_, length(off))

  if (criterion == "posterior") {
    po <- parentage_posterior(.as_lod_matrix(tab), prop_sampled = prop_sampled,
                              n_candidates = length(candidates))
    m <- match(best$offspring, po$offspring)
    best$candidate <- po$candidate[m]
    best$lod <- po$lod[m]
    best$posterior <- po$posterior[m]
    best$post_unsampled <- po$post_unsampled[m]
    pass <- vapply(confidence, function(cf)
      !is.na(best$posterior) & best$posterior >= cf, logical(nrow(best)))
    dim(pass) <- c(nrow(best), length(confidence))
    criteria <- data.frame(
      confidence = confidence,
      n_at_or_above = colSums(pass),
      expected_false = vapply(seq_along(confidence), function(j)
        sum(1 - best$posterior[pass[, j]]), numeric(1)))
    attr(best, "fdr_5pct") <- parentage_fdr(best$posterior, target = 0.05)
  } else if (criterion == "pairwise") {
    crit <- matrix(NA_real_, nrow(best), length(confidence),
                   dimnames = list(NULL, paste0("crit_", confidence)))
    for (i in seq_len(nrow(best))) {
      j <- match(best$offspring[i], off)
      nd <- lod_null(g, offspring = best$offspring[i], freqs = freqs, error = error,
                     known = if (is.na(kn[j])) NULL else kn[j],
                     method = "backward", nsim = nsim, prop_typed = prop_typed,
                     type = type)
      crit[i, ] <- lod_critical(nd, alpha = 1 - confidence)$critical_lod
    }
    best <- cbind(best, as.data.frame(crit))
    pass <- vapply(seq_along(confidence), function(j) best$lod >= crit[, j],
                   logical(nrow(best)))
    dim(pass) <- c(nrow(best), length(confidence))
    criteria <- data.frame(confidence = confidence,
                           critical = "per offspring, backward simulation",
                           stringsAsFactors = FALSE)
  } else {
    if (is.null(prop_typed)) prop_typed <- mean(!is.na(g$a1))
    ds <- sim_delta(freqs, error = error, n_candidates = length(candidates),
                    prop_sampled = prop_sampled, prop_typed = prop_typed,
                    nsim = nsim, type = type)
    dc <- delta_critical(ds, levels = confidence)
    for (j in seq_along(confidence)) {
      best[[paste0("crit_", confidence[j])]] <- dc$critical_delta[j]
    }
    pass <- vapply(seq_along(confidence), function(j) {
      cd <- dc$critical_delta[j]
      if (is.na(cd)) rep(FALSE, nrow(best)) else !is.na(best$delta) & best$delta >= cd
    }, logical(nrow(best)))
    dim(pass) <- c(nrow(best), length(confidence))
    criteria <- dc
    attr(best, "delta_sim") <- ds
  }

  lvl <- rep(NA_real_, nrow(best))
  for (j in seq_along(confidence)) lvl[pass[, j]] <- confidence[j]
  best$confidence <- lvl
  best$assigned <- !is.na(lvl)

  attr(best, "lod_table") <- tab
  attr(best, "criteria") <- criteria
  attr(best, "criterion") <- criterion
  attr(best, "type") <- type
  class(best) <- c("parentage_assignment", "data.frame")
  best
}

#' @param x A `"parentage_assignment"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname assign_parentage
#' @export
print.parentage_assignment <- function(x, ...) {
  cat("<parentage_assignment>", nrow(x), "offspring, criterion =",
      attr(x, "criterion"), ", type =", attr(x, "type"), "\n")
  n <- sum(x$assigned)
  cat("  assigned: ", n, " of ", nrow(x), " (",
      format(100 * n / max(nrow(x), 1L), digits = 3), "%)\n", sep = "")
  print.data.frame(as.data.frame(x), digits = 4, row.names = FALSE, ...)
  invisible(x)
}

#' @param object A `"parentage_assignment"` object.
#' @rdname assign_parentage
#' @return `summary()` returns the criterion table invisibly after printing a
#'   short report.
#' @export
summary.parentage_assignment <- function(object, ...) {
  cat("Parentage assignment\n")
  cat("  criterion   : ", attr(object, "criterion"), "\n", sep = "")
  cat("  configuration: ", attr(object, "type"), "\n", sep = "")
  cat("  offspring   : ", nrow(object), "\n", sep = "")
  cat("  assigned    : ", sum(object$assigned), "\n", sep = "")
  tb <- table(factor(object$confidence))
  if (length(tb)) {
    cat("  by confidence level:\n")
    for (nm in names(tb)) cat("    ", nm, ": ", tb[[nm]], "\n", sep = "")
  }
  cat("\nCriteria:\n")
  print(attr(object, "criteria"))
  invisible(attr(object, "criteria"))
}
