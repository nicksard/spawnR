#' Population-level simulation of the Delta criterion
#'
#' Reproduces the simulation of Marshall et al. (1998), which underlies the
#' confidence levels reported by CERVUS. A true parent and an offspring are
#' generated under Hardy-Weinberg proportions, a pool of unrelated candidates is
#' added, the true parent is included in that pool only with probability
#' `prop_sampled`, genotyping error and missing data are applied, and the Delta
#' statistic of the resulting assignment is recorded together with whether the
#' assignment was correct.
#'
#' The output feeds [delta_critical()], which converts the simulated
#' Delta-and-correctness pairs into critical values at the conventional relaxed
#' (80%) and strict (95%) confidence levels.
#'
#' This is the population-averaged alternative to the trio-specific criterion of
#' [lod_null()]. Both are provided so that the two can be compared on the same
#' data, which is the comparison Amiri Roudbar et al. (2025) make.
#'
#' @param freqs An `"allele_freqs"` object supplying the reference frequencies.
#' @param error Genotyping error rate assumed when computing LOD scores, and
#'   used to generate errors, one value or one per locus.
#' @param n_candidates Number of candidate parents offered for each offspring.
#' @param prop_sampled Probability that the true parent is among the candidates.
#' @param prop_typed Probability that any given genotype is scored.
#' @param nsim Number of simulated assignments.
#' @param type `"both_unknown"` or `"one_known"`; for the latter a known second
#'   parent is generated and used in the likelihood.
#' @param assumed_error Error rate used in the likelihood if it should differ
#'   from the rate used to generate errors. Defaults to `error`.
#'
#' @return An object of class `"delta_sim"`: a `data.frame` with columns `delta`,
#'   `lod` (of the top candidate), `correct` (whether the top candidate is the
#'   true parent) and `sampled` (whether the true parent was in the pool).
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
#' @examples
#' sim <- simulate_population(n = 30, n_loci = 8, n_alleles = 6, seed = 5)
#' p <- allele_freqs(sim$genotypes)
#' ds <- sim_delta(p, error = 0.01, n_candidates = 20, prop_sampled = 0.9,
#'                 nsim = 300)
#' delta_critical(ds)
#' @export
sim_delta <- function(freqs, error, n_candidates = 100L, prop_sampled = 1,
                      prop_typed = 1, nsim = 10000L,
                      type = c("both_unknown", "one_known"),
                      assumed_error = NULL) {
  type <- match.arg(type)
  nl <- length(freqs)
  e_true <- .expand_error(error, nl)
  e_used <- if (is.null(assumed_error)) e_true else .expand_error(assumed_error, nl, "assumed_error")
  nsim <- as.integer(nsim); nc <- as.integer(n_candidates)
  if (nsim < 1L || nc < 1L) stop("`nsim` and `n_candidates` must be positive.", call. = FALSE)
  if (prop_sampled < 0 || prop_sampled > 1) stop("`prop_sampled` must lie in [0, 1].", call. = FALSE)

  sampled <- stats::runif(nsim) < prop_sampled
  lodsum <- matrix(0, nsim, nc)
  n_comp <- matrix(0L, nsim, nc)

  for (l in seq_len(nl)) {
    p <- freqs[[l]]
    k <- length(p)
    dr <- function(m) sample.int(k, m, TRUE, prob = p)

    f1 <- dr(nsim); f2 <- dr(nsim)                       # true parent
    if (type == "one_known") { k1 <- dr(nsim); k2 <- dr(nsim) } else { k1 <- k2 <- NULL }
    o_from_f <- ifelse(stats::runif(nsim) < 0.5, f1, f2)
    o_from_o <- if (type == "one_known") ifelse(stats::runif(nsim) < 0.5, k1, k2) else dr(nsim)

    # candidate pool: column 1 is the true parent when he was sampled
    C1 <- matrix(dr(nsim * nc), nsim, nc)
    C2 <- matrix(dr(nsim * nc), nsim, nc)
    C1[sampled, 1L] <- f1[sampled]
    C2[sampled, 1L] <- f2[sampled]

    og <- .apply_missing(.apply_error(o_from_f, o_from_o, e_true[l], p), prop_typed)
    cg <- .apply_missing(.apply_error(as.vector(C1), as.vector(C2), e_true[l], p), prop_typed)
    C1 <- matrix(cg$a1, nsim, nc); C2 <- matrix(cg$a2, nsim, nc)
    if (type == "one_known") {
      kg <- .apply_missing(.apply_error(k1, k2, e_true[l], p), prop_typed)
      k1 <- kg$a1; k2 <- kg$a2
    }

    # The contribution depends on the candidate only through his genotype, so
    # evaluate the likelihood once per distinct genotype (G of them) over the
    # nsim simulated offspring, then reduce the nsim-by-nc block to a lookup.
    tb <- .locus_tab(k)
    Vl <- matrix(0, nsim, tb$G + 1L)
    Vn <- matrix(0L, nsim, tb$G + 1L)
    for (a in seq_len(tb$G)) {
      v <- lod_locus(og$a1, og$a2, tb$gi[a], tb$gj[a], p = p, error = e_used[l],
                     k1 = k1, k2 = k2, type = type)
      fin <- is.finite(v)
      v[!fin] <- 0
      Vl[, a] <- v
      Vn[, a] <- as.integer(fin)
    }
    ccode <- match(pmin(C1, C2) + (pmax(C1, C2) - 1L) * k, tb$gi + (tb$gj - 1L) * k)
    ccode[is.na(ccode)] <- tb$G + 1L
    idx <- cbind(rep(seq_len(nsim), nc), as.vector(ccode))
    lodsum <- lodsum + matrix(Vl[idx], nsim, nc)
    n_comp <- n_comp + matrix(Vn[idx], nsim, nc)
  }
  lodsum[n_comp == 0L] <- -Inf

  top  <- max.col(lodsum, ties.method = "first")
  best <- lodsum[cbind(seq_len(nsim), top)]
  tmp <- lodsum
  tmp[cbind(seq_len(nsim), top)] <- -Inf
  second <- if (nc >= 2L) apply(tmp, 1L, max) else rep(-Inf, nsim)

  npos <- rowSums(lodsum > 0)
  delta <- ifelse(npos >= 2L, best - second, ifelse(npos == 1L, best, NA_real_))
  delta[!is.finite(best)] <- NA_real_

  out <- data.frame(delta = delta, lod = ifelse(is.finite(best), best, NA_real_),
                    correct = sampled & top == 1L & npos >= 1L,
                    sampled = sampled)
  attr(out, "settings") <- list(n_candidates = nc, prop_sampled = prop_sampled,
                                prop_typed = prop_typed, nsim = nsim, type = type,
                                error = e_true)
  class(out) <- c("delta_sim", "data.frame")
  out
}

#' Critical Delta values at stated confidence
#'
#' Applies the definition of confidence used by Marshall et al. (1998): the
#' confidence attached to a threshold is the proportion of simulated assignments
#' with Delta at or above that threshold in which the top candidate really is the
#' true parent. The critical value at a given confidence level is therefore the
#' most permissive threshold that still attains it. The conventional relaxed and
#' strict levels of 80% and 95% are the defaults.
#'
#' @param x A `"delta_sim"` object from [sim_delta()].
#' @param levels Confidence levels.
#' @param min_n Minimum number of simulated assignments that must lie above a
#'   threshold before it is eligible, which stops the extreme tail of the
#'   simulation from producing an unstable criterion.
#' @param ... Ignored.
#'
#' @return A `data.frame` with columns `confidence`, `critical_delta`,
#'   `n_above` and `prop_correct`. `critical_delta` is `NA` when no threshold
#'   attains the level, which happens when too few of the true parents were
#'   sampled or the markers are too weak.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' @examples
#' sim <- simulate_population(n = 30, n_loci = 8, n_alleles = 6, seed = 5)
#' p <- allele_freqs(sim$genotypes)
#' ds <- sim_delta(p, error = 0.01, n_candidates = 20, nsim = 300)
#' delta_critical(ds, levels = c(0.80, 0.95))
#' @export
delta_critical <- function(x, levels = c(0.80, 0.95), min_n = 20L, ...) {
  stopifnot(inherits(x, "delta_sim"))
  d <- x[!is.na(x$delta), , drop = FALSE]
  if (nrow(d) == 0L) {
    return(data.frame(confidence = levels, critical_delta = NA_real_,
                      n_above = 0L, prop_correct = NA_real_))
  }
  o <- order(d$delta, decreasing = TRUE)
  dd <- d$delta[o]
  cum <- cumsum(d$correct[o]) / seq_along(o)
  res <- lapply(levels, function(lv) {
    ok <- which(cum >= lv & seq_along(cum) >= min_n)
    if (length(ok) == 0L) {
      data.frame(confidence = lv, critical_delta = NA_real_, n_above = 0L,
                 prop_correct = NA_real_)
    } else {
      i <- max(ok)
      data.frame(confidence = lv, critical_delta = dd[i], n_above = i,
                 prop_correct = cum[i])
    }
  })
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

#' @param x A `"delta_sim"` object.
#' @param ... Ignored.
#' @rdname sim_delta
#' @export
print.delta_sim <- function(x, ...) {
  s <- attr(x, "settings")
  cat("<delta_sim>", s$nsim, "assignments,", s$n_candidates, "candidates, type =", s$type, "\n")
  cat("  true parent sampled: ", format(100 * mean(x$sampled), digits = 3), "%",
      ";  top candidate correct: ", format(100 * mean(x$correct), digits = 3), "%\n", sep = "")
  cat("  Delta: median ", format(stats::median(x$delta, na.rm = TRUE), digits = 4),
      ", max ", format(max(x$delta, na.rm = TRUE), digits = 4), "\n", sep = "")
  invisible(x)
}
