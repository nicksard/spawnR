#' Assign mothers, fathers and parent trios with Delta confidence
#'
#' One call for the common design in which candidate mothers and candidate
#' fathers are sampled separately and every offspring needs its most likely
#' mother, its most likely father, and its most likely parent pair, each with a
#' statement of statistical confidence and a count of the loci at which each
#' parent is excluded outright.
#'
#' Three analyses are run on the same data:
#'
#' * **Mother alone** - every candidate mother scored against the offspring
#'   with the father unknown (`type = "both_unknown"` in [parentage_lod()]).
#' * **Father alone** - the same for every candidate father.
#' * **Trio** - every mother-by-father combination scored jointly
#'   (`type = "pair"`), so the best pair is chosen on the trio likelihood
#'   rather than by pasting together the two best single parents. The two
#'   usually agree; where they do not, the trio is the better-supported answer.
#'
#' For each, Delta is the gap between the best and second-best LOD among
#' candidates with a positive score (see [delta_stat()]), and significance is
#' the population-level criterion of Marshall et al. (1998): [sim_delta()] is
#' run once per analysis with the pool sizes and sampling proportions given
#' here, and [delta_critical()] converts the simulation into a critical Delta
#' at each `confidence` level. An assignment is significant at a level when its
#' Delta reaches the critical value, and `*_confidence` holds the highest level
#' attained. For the trio, a simulated assignment only counts as correct when
#' both parents are right, so the trio confidence is a statement about the pair.
#'
#' @section Mismatches:
#' The `*_mismatches` columns count loci at which the named parent is
#' incompatible with the offspring under Mendelian inheritance, ignoring
#' genotyping error: exclusion alone, separately for each parent. For the trio,
#' `trio_mother_mismatches` and `trio_father_mismatches` give each parent's own
#' count and `trio_mismatches` counts loci at which the three genotypes cannot
#' be reconciled jointly, which can exceed either parent's count (each parent
#' may carry one of the offspring's alleles, but not the same one). These
#' counts are reported for screening and do not enter the likelihood, which
#' already allows for error; use `max_mismatch` to drop candidates outright.
#'
#' @param g A `"genotypes"` object holding the offspring and all candidates.
#' @param offspring Identifiers of the offspring.
#' @param mothers Identifiers of the candidate mothers.
#' @param fathers Identifiers of the candidate fathers.
#' @param freqs An `"allele_freqs"` object.
#' @param error Genotyping error rate, one value or one per locus.
#' @param confidence Confidence levels at which to test Delta. The defaults are
#'   the conventional relaxed and strict levels reported by CERVUS.
#' @param prop_sampled_mothers,prop_sampled_fathers Proportion of the true
#'   mothers and fathers assumed to be among the candidates. These enter the
#'   significance simulation and matter: confidence is lower when the true
#'   parent may be missing.
#' @param prop_typed Proportion of genotypes scored, for the simulation.
#'   Defaults to the observed proportion among the offspring, mothers and
#'   fathers analysed.
#' @param nsim Replicates for each single-parent simulation.
#' @param nsim_trio Replicates for the trio simulation, whose cost grows with
#'   the number of mothers times the number of fathers.
#' @param max_mismatch If not `NULL`, candidates (and, for the trio, pairs)
#'   mismatching at more than this many loci are excluded before the best is
#'   chosen and Delta computed. The same screen is applied inside the
#'   significance simulations.
#' @param min_loci If not `NULL`, a candidate (or pair) compared with the
#'   offspring at fewer than this many loci is excluded, so an offspring typed
#'   at too few loci is left unassigned. Also applied inside the significance
#'   simulations. CERVUS applies a rule of this kind; a common choice is half
#'   the loci.
#' @param output `"best"` for one row per offspring, or `"all"` for every
#'   mother-by-father combination of every offspring, ranked.
#'
#' @return With `output = "best"`, an object of class `"trio_assignment"`: a
#'   `data.frame` with one row per offspring and, for each of the prefixes
#'   `mother_`, `father_` and `trio_`, the best candidate(s), `lod`, `delta`,
#'   `n_compared`, mismatch counts and `confidence`, plus a logical
#'   `*_significant` at the lowest level tested. The critical values are
#'   attached as attribute `criteria` and the simulations as `delta_sims`.
#'
#'   With `output = "all"`, a `data.frame` with one row per offspring, mother
#'   and father: `trio_lod` and its `rank` within the offspring, the
#'   single-parent `mother_lod` and `father_lod`, the three mismatch counts and
#'   `n_compared`. The top-ranked row of each offspring also carries `delta`,
#'   `confidence` and `significant`. Rows are sorted by offspring and rank, and
#'   with many candidates this table is large; `max_mismatch` keeps it small.
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
#' sim <- simulate_population(n = 40, n_loci = 10, n_alleles = 8, seed = 3)
#' p <- allele_freqs(sim$genotypes)
#' a <- assign_trios(sim$genotypes, sim$offspring[1:6], mothers = sim$dams,
#'                   fathers = sim$sires, freqs = p, error = 0.01,
#'                   nsim = 500, nsim_trio = 200)
#' a
#' summary(a)
#'
#' # every candidate pair for the first offspring
#' all <- assign_trios(sim$genotypes, sim$offspring[1], sim$dams, sim$sires,
#'                     p, error = 0.01, nsim = 500, nsim_trio = 200,
#'                     output = "all")
#' head(all)
#' @export
assign_trios <- function(g, offspring, mothers, fathers, freqs, error,
                         confidence = c(0.80, 0.95),
                         prop_sampled_mothers = 1, prop_sampled_fathers = 1,
                         prop_typed = NULL, nsim = 10000L, nsim_trio = 2000L,
                         max_mismatch = NULL, min_loci = NULL,
                         output = c("best", "all")) {
  stopifnot(inherits(g, "genotypes"))
  output <- match.arg(output)
  confidence <- sort(confidence)
  off <- .as_ids(g, offspring, "offspring")
  mom <- unique(.as_ids(g, mothers, "mothers"))
  dad <- unique(.as_ids(g, fathers, "fathers"))
  if (!length(off) || !length(mom) || !length(dad)) {
    stop("`offspring`, `mothers` and `fathers` must each be non-empty.", call. = FALSE)
  }
  nm <- length(mom); nf <- length(dad)

  # LOD, mismatch and comparability matrices, offspring in rows. Pair columns
  # run over fathers fastest: column i + (j - 1) * nf is father i, mother j.
  mo <- parentage_lod(g, off, mom, freqs, error, type = "both_unknown", output = "matrix")
  fa <- parentage_lod(g, off, dad, freqs, error, type = "both_unknown", output = "matrix")
  pr <- parentage_lod(g, off, dad, freqs, error, type = "pair", mothers = mom,
                     output = "matrix")
  pm <- rep(mom, each = nf); pf <- rep(dad, times = nm)

  mask <- function(L, mis, ncmp, self) {
    L[self] <- NA_real_
    if (!is.null(max_mismatch)) L[mis > max_mismatch] <- NA_real_
    if (!is.null(min_loci)) L[ncmp < min_loci] <- NA_real_
    L
  }
  Lm <- mask(mo$lod, mo$mismatches, mo$n_compared, outer(off, mom, `==`))
  Lf <- mask(fa$lod, fa$mismatches, fa$n_compared, outer(off, dad, `==`))
  Lp <- mask(pr$lod, pr$mismatches, pr$n_compared,
             outer(off, pm, `==`) | outer(off, pf, `==`) |
               matrix(pm == pf, length(off), nf * nm, byrow = TRUE))

  # Significance: one population-level simulation per analysis.
  if (is.null(prop_typed)) prop_typed <- .prop_typed(g, c(off, mom, dad))
  sims <- list(
    mother = sim_delta(freqs, error, n_candidates = nm, prop_sampled = prop_sampled_mothers,
                       prop_typed = prop_typed, nsim = nsim,
                       max_mismatch = max_mismatch, min_loci = min_loci),
    father = sim_delta(freqs, error, n_candidates = nf, prop_sampled = prop_sampled_fathers,
                       prop_typed = prop_typed, nsim = nsim,
                       max_mismatch = max_mismatch, min_loci = min_loci),
    trio = sim_delta(freqs, error, n_candidates = nf, n_candidates2 = nm,
                     prop_sampled = prop_sampled_fathers,
                     prop_sampled2 = prop_sampled_mothers,
                     prop_typed = prop_typed, nsim = nsim_trio, type = "pair",
                     max_mismatch = max_mismatch, min_loci = min_loci))
  crit <- lapply(sims, delta_critical, levels = confidence)
  criteria <- do.call(rbind, lapply(names(crit), function(nm_)
    cbind(analysis = nm_, crit[[nm_]], stringsAsFactors = FALSE)))
  rownames(criteria) <- NULL

  bm <- .best_rows(Lm); bf <- .best_rows(Lf); bp <- .best_rows(Lp)
  cm <- .conf_level(bm$delta, crit$mother)
  cf <- .conf_level(bf$delta, crit$father)
  cp <- .conf_level(bp$delta, crit$trio)

  if (output == "all") {
    return(.all_trios(off, pm, pf, Lp, pr, mo, fa, mom, dad, nf, bp, cp, confidence,
                      criteria, sims))
  }

  r <- seq_along(off)
  at <- function(X, j) X[cbind(r, ifelse(is.na(j), 1L, j))] * ifelse(is.na(j), NA, 1)
  tm <- match(pm[bp$top], mom); tf <- match(pf[bp$top], dad)
  out <- data.frame(
    offspring = off,
    mother = mom[bm$top], mother_lod = bm$lod, mother_delta = bm$delta,
    mother_n_compared = as.integer(at(mo$n_compared, bm$top)),
    mother_mismatches = as.integer(at(mo$mismatches, bm$top)),
    mother_confidence = cm, mother_significant = .sig(cm),
    father = dad[bf$top], father_lod = bf$lod, father_delta = bf$delta,
    father_n_compared = as.integer(at(fa$n_compared, bf$top)),
    father_mismatches = as.integer(at(fa$mismatches, bf$top)),
    father_confidence = cf, father_significant = .sig(cf),
    trio_mother = pm[bp$top], trio_father = pf[bp$top],
    trio_lod = bp$lod, trio_delta = bp$delta,
    trio_n_compared = as.integer(at(pr$n_compared, bp$top)),
    trio_mother_mismatches = as.integer(at(mo$mismatches, tm)),
    trio_father_mismatches = as.integer(at(fa$mismatches, tf)),
    trio_mismatches = as.integer(at(pr$mismatches, bp$top)),
    trio_confidence = cp, trio_significant = .sig(cp),
    stringsAsFactors = FALSE)
  attr(out, "criteria") <- criteria
  attr(out, "delta_sims") <- sims
  attr(out, "confidence") <- confidence
  class(out) <- c("trio_assignment", "data.frame")
  out
}

# Best column of each row of a LOD matrix (NA = not a contender), with the
# runner-up and Delta defined exactly as in delta_stat().
.best_rows <- function(L) {
  Z <- L; Z[is.na(Z)] <- -Inf
  n <- nrow(Z)
  top <- max.col(Z, ties.method = "first")
  best <- Z[cbind(seq_len(n), top)]
  npos <- rowSums(Z > 0)
  second <- rep(-Inf, n)
  if (ncol(Z) >= 2L) {
    Z[cbind(seq_len(n), top)] <- -Inf
    second <- Z[cbind(seq_len(n), max.col(Z, ties.method = "first"))]
  }
  delta <- ifelse(npos >= 2L, best - second, ifelse(npos == 1L, best, NA_real_))
  none <- !is.finite(best)
  top[none] <- NA_integer_
  list(top = top, lod = ifelse(none, NA_real_, best),
       lod2 = ifelse(npos >= 2L, second, NA_real_), delta = delta, npos = npos)
}

# Highest confidence level whose critical Delta is reached; NA if none.
.conf_level <- function(delta, crit) {
  lvl <- rep(NA_real_, length(delta))
  for (j in seq_len(nrow(crit))) {
    cd <- crit$critical_delta[j]
    if (!is.na(cd)) lvl[!is.na(delta) & delta >= cd] <- crit$confidence[j]
  }
  lvl
}

.sig <- function(conf) !is.na(conf)

.all_trios <- function(off, pm, pf, Lp, pr, mo, fa, mom, dad, nf, bp, cp,
                       confidence, criteria, sims) {
  no <- length(off); np <- length(pm)
  im <- match(pm, mom); ifa <- match(pf, dad)
  d <- data.frame(offspring = rep(off, each = np),
                  mother = rep(pm, times = no), father = rep(pf, times = no),
                  stringsAsFactors = FALSE)
  d$trio_lod <- as.vector(t(Lp))
  d$mother_lod <- as.vector(t(mo$lod[, im, drop = FALSE]))
  d$father_lod <- as.vector(t(fa$lod[, ifa, drop = FALSE]))
  d$n_compared <- as.integer(t(pr$n_compared))
  d$mother_mismatches <- as.integer(t(mo$mismatches[, im, drop = FALSE]))
  d$father_mismatches <- as.integer(t(fa$mismatches[, ifa, drop = FALSE]))
  d$trio_mismatches <- as.integer(t(pr$mismatches))
  d <- d[!is.na(d$trio_lod), , drop = FALSE]
  d <- d[order(match(d$offspring, off), -d$trio_lod), , drop = FALSE]
  d$rank <- stats::ave(d$trio_lod, d$offspring, FUN = seq_along)
  first <- d$rank == 1
  k <- match(d$offspring, off)
  d$delta <- ifelse(first, bp$delta[k], NA_real_)
  d$confidence <- ifelse(first, cp[k], NA_real_)
  d$significant <- ifelse(first, .sig(cp[k]), NA)
  d <- d[c("offspring", "rank", "mother", "father", "trio_lod", "delta",
           "confidence", "significant", "mother_lod", "father_lod", "n_compared",
           "mother_mismatches", "father_mismatches", "trio_mismatches")]
  d$rank <- as.integer(d$rank)
  rownames(d) <- NULL
  attr(d, "criteria") <- criteria
  attr(d, "delta_sims") <- sims
  d
}

#' @param x A `"trio_assignment"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname assign_trios
#' @export
print.trio_assignment <- function(x, ...) {
  conf <- attr(x, "confidence")
  cat("<trio_assignment>", nrow(x), "offspring\n")
  for (w in c("mother", "father", "trio")) {
    n <- sum(x[[paste0(w, "_significant")]])
    cat(sprintf("  %-7s significant at >= %g%%: %d of %d\n", w,
                100 * min(conf), n, nrow(x)))
  }
  show <- data.frame(
    offspring = x$offspring,
    mother = x$mother, m_lod = x$mother_lod, m_delta = x$mother_delta,
    m_mis = x$mother_mismatches, m_conf = .conf_label(x$mother_confidence),
    father = x$father, f_lod = x$father_lod, f_delta = x$father_delta,
    f_mis = x$father_mismatches, f_conf = .conf_label(x$father_confidence),
    trio_mother = x$trio_mother, trio_father = x$trio_father,
    t_lod = x$trio_lod, t_delta = x$trio_delta,
    t_mis = paste(x$trio_mother_mismatches, x$trio_father_mismatches,
                  x$trio_mismatches, sep = "/"),
    t_conf = .conf_label(x$trio_confidence), stringsAsFactors = FALSE)
  print.data.frame(show, digits = 4, row.names = FALSE, ...)
  cat("  t_mis = mother/father/trio mismatching loci\n")
  invisible(x)
}

.conf_label <- function(conf) ifelse(is.na(conf), "-", paste0(100 * conf, "%"))

#' @param object A `"trio_assignment"` object.
#' @rdname assign_trios
#' @export
summary.trio_assignment <- function(object, ...) {
  conf <- attr(object, "confidence")
  cat("Mother, father and trio assignment\n")
  cat("  offspring: ", nrow(object), "\n", sep = "")
  tab <- vapply(c("mother", "father", "trio"), function(w) {
    cc <- object[[paste0(w, "_confidence")]]
    vapply(conf, function(lv) sum(!is.na(cc) & cc >= lv), numeric(1))
  }, numeric(length(conf)))
  tab <- matrix(tab, nrow = length(conf),
                dimnames = list(paste0(">= ", 100 * conf, "%"),
                                c("mother", "father", "trio")))
  cat("  assigned, by confidence level:\n")
  print(tab)
  agree <- object$mother == object$trio_mother & object$father == object$trio_father
  cat("  trio agrees with the two best single parents: ",
      sum(agree, na.rm = TRUE), " of ", nrow(object), "\n", sep = "")
  cat("\nCritical Delta:\n")
  print(attr(object, "criteria"), row.names = FALSE)
  invisible(attr(object, "criteria"))
}
