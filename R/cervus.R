#' Read a CERVUS parentage results file
#'
#' Reads the comma-separated parentage output CERVUS writes and returns it in
#' the column names this package uses, so that [compare_cervus()] can hold the
#' two implementations against each other on the same data.
#'
#' @section Column detection:
#' CERVUS's exact header text varies between analysis types and versions, and no
#' authoritative list is published, so columns are matched by pattern rather than
#' assumed. Matching ignores case, spaces and punctuation. The mapping that was
#' used is returned as the attribute `"mapping"` and printed, so you can check it
#' against your file rather than trust it. Anything not recognised is kept
#' unchanged.
#'
#' @section Pair against trio:
#' CERVUS reports a *pair* LOD for an offspring against one candidate parent and
#' a *trio* LOD when the other parent is included. These correspond to
#' `type = "both_unknown"` and `type = "one_known"` in [parentage_lod()]. Compare
#' like with like: a trio LOD will not match a pair LOD and the difference is not
#' a bug in either program.
#'
#' @param file Path to the CERVUS output file.
#' @param which Which score to extract when the file holds both, `"pair"` or
#'   `"trio"`. Defaults to whichever is present, preferring `"pair"`.
#' @param ... Passed to [utils::read.csv()].
#'
#' @return A `data.frame` with whichever of `offspring`, `candidate`,
#'   `n_compared`, `mismatches`, `lod`, `delta` and `confidence` could be found,
#'   plus any unrecognised columns. The detected mapping is attached as the
#'   attribute `"mapping"`.
#'
#' @seealso [compare_cervus()]
#'
#' @references
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x}
#'
#' @export
read_cervus_results <- function(file, which = c("pair", "trio"), ...) {
  which <- match.arg(which)
  d <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
  key <- tolower(gsub("[^a-z0-9]", "", tolower(names(d))))

  pick <- function(...) {
    pats <- c(...)
    for (p in pats) {
      i <- grep(p, key)
      if (length(i)) return(i[1L])
    }
    NA_integer_
  }
  pre <- if (which == "trio") "trio" else "pair"
  alt <- if (which == "trio") "pair" else "trio"

  map <- c(
    offspring  = pick("^offspringid$", "^offspring$", "^offspringid"),
    candidate  = pick("^candidate.*id$", "^candidate", "^father", "^mother", "^parentid$"),
    n_compared = pick(paste0("^", pre, "locicompared$"), "locicompared",
                      paste0("^", alt, "locicompared$")),
    mismatches = pick(paste0("^", pre, "locimismatch"), "locimismatch",
                      paste0("^", alt, "locimismatch")),
    lod        = pick(paste0("^", pre, "lodscore$"), paste0("^", pre, "lod$"),
                      "lodscore", "^lod$"),
    delta      = pick(paste0("^", pre, "delta$"), "^delta$"),
    confidence = pick(paste0("^", pre, "confidence$"), "confidence")
  )
  map <- map[!is.na(map)]
  if (!all(c("offspring", "candidate", "lod") %in% names(map))) {
    stop("Could not find offspring, candidate and LOD columns in '", file,
         "'. Columns present: ", paste(names(d), collapse = ", "),
         ". Rename them or pass the file through read.csv() and build the ",
         "columns by hand.", call. = FALSE)
  }
  out <- d[, map, drop = FALSE]
  names(out) <- names(map)
  out$offspring <- trimws(as.character(out$offspring))
  out$candidate <- trimws(as.character(out$candidate))
  for (nm in intersect(c("n_compared", "mismatches", "lod", "delta"), names(out))) {
    out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  }
  out <- out[nzchar(out$offspring) & nzchar(out$candidate) & !is.na(out$lod), , drop = FALSE]
  keep <- setdiff(seq_along(d), map)
  if (length(keep)) out <- cbind(out, d[rownames(out), keep, drop = FALSE])
  rownames(out) <- NULL
  attr(out, "mapping") <- stats::setNames(names(d)[map], names(map))
  attr(out, "which") <- which
  message("Mapped CERVUS columns (", which, "): ",
          paste(sprintf("%s <- '%s'", names(attr(out, "mapping")),
                        attr(out, "mapping")), collapse = ", "))
  out
}

#' Compare spawnR LOD scores against CERVUS
#'
#' Joins this package's scores to CERVUS's on offspring and candidate, and
#' reports how closely they agree. This is the external check that matters: the
#' error-model coefficients in this package were re-derived rather than
#' transcribed, because the ones printed in the literature do not form a
#' probability distribution, and the corrigendum to Kalinowski et al. (2007)
#' states that those typesetting errors never affected the CERVUS software. So
#' CERVUS implements the intended equations, and agreement with it is evidence
#' the re-derivation recovered them.
#'
#' @section Getting the inputs to match:
#' Four things must line up before a difference means anything.
#' \enumerate{
#'   \item **Allele frequencies.** Both programs must use the same ones.
#'     [allele_freqs()] floors frequencies at `1/(2n)` by default, which CERVUS
#'     does not do; pass `min_freq = 0` to turn that off, and check CERVUS's own
#'     minimum allele frequency setting.
#'   \item **Error rate.** Use the single rate CERVUS was given, not a per-locus
#'     vector, since CERVUS assumes one rate across loci.
#'   \item **Configuration.** A CERVUS pair LOD is `type = "both_unknown"`; a
#'     trio LOD with the other parent supplied is `type = "one_known"`.
#'   \item **Alternative hypothesis.** Leave `relatedness` and `ibd` unset.
#'     CERVUS's alternative is an unrelated candidate, so a relatedness-aware
#'     score is deliberately different and will not match.
#' }
#'
#' @param x A `"parentage_lod"` long table from [parentage_lod()], or any
#'   `data.frame` with `offspring`, `candidate` and `lod`.
#' @param cervus A table from [read_cervus_results()].
#' @param tol Absolute LOD difference above which a row is flagged.
#'
#' @return An object of class `"cervus_comparison"`: a list with `summary` (a
#'   one-row `data.frame` of counts, correlation and difference statistics),
#'   `joined` (the merged table with a `diff` column) and `worst` (the rows
#'   exceeding `tol`, worst first).
#'
#' @references
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x} Corrigendum:
#' \emph{Molecular Ecology}, \strong{19}, 1512. \doi{10.1111/j.1365-294x.2010.04544.x}
#'
#' @examples
#' # Standing in for a CERVUS export, to show the shape of the comparison.
#' sim <- simulate_population(n = 60, n_loci = 10, n_alleles = 6,
#'                            n_offspring = 10, error = 0.01, seed = 3)
#' p <- allele_freqs(sim$genotypes, min_freq = 0)
#' mine <- parentage_lod(sim$genotypes, sim$pedigree$offspring, sim$sires,
#'                       p, error = 0.01)
#' fake <- data.frame(offspring = mine$offspring, candidate = mine$candidate,
#'                    lod = mine$lod, mismatches = mine$mismatches)
#' compare_cervus(mine, fake)
#' @export
compare_cervus <- function(x, cervus, tol = 0.01) {
  need <- c("offspring", "candidate", "lod")
  if (!all(need %in% names(x))) stop("`x` needs columns: ", paste(need, collapse = ", "), call. = FALSE)
  if (!all(need %in% names(cervus))) stop("`cervus` needs columns: ", paste(need, collapse = ", "), call. = FALSE)
  a <- data.frame(offspring = as.character(x$offspring),
                  candidate = as.character(x$candidate),
                  lod_spawnR = as.numeric(x$lod),
                  mm_spawnR = if (!is.null(x$mismatches)) x$mismatches else NA_integer_,
                  nc_spawnR = if (!is.null(x$n_compared)) x$n_compared else NA_integer_,
                  stringsAsFactors = FALSE)
  b <- data.frame(offspring = as.character(cervus$offspring),
                  candidate = as.character(cervus$candidate),
                  lod_cervus = as.numeric(cervus$lod),
                  mm_cervus = if (!is.null(cervus$mismatches)) cervus$mismatches else NA_integer_,
                  nc_cervus = if (!is.null(cervus$n_compared)) cervus$n_compared else NA_integer_,
                  stringsAsFactors = FALSE)
  j <- merge(a, b, by = c("offspring", "candidate"))
  j$diff <- j$lod_spawnR - j$lod_cervus
  ok <- is.finite(j$diff)

  mm_agree <- if (all(is.na(j$mm_spawnR)) || all(is.na(j$mm_cervus))) NA_real_
              else mean(j$mm_spawnR == j$mm_cervus, na.rm = TRUE)
  nc_agree <- if (all(is.na(j$nc_spawnR)) || all(is.na(j$nc_cervus))) NA_real_
              else mean(j$nc_spawnR == j$nc_cervus, na.rm = TRUE)

  s <- data.frame(
    n_spawnR = nrow(a), n_cervus = nrow(b), n_matched = nrow(j),
    correlation = if (sum(ok) > 2) stats::cor(j$lod_spawnR[ok], j$lod_cervus[ok]) else NA_real_,
    mean_abs_diff = if (any(ok)) mean(abs(j$diff[ok])) else NA_real_,
    max_abs_diff = if (any(ok)) max(abs(j$diff[ok])) else NA_real_,
    n_over_tol = sum(abs(j$diff) > tol, na.rm = TRUE),
    mismatch_agreement = mm_agree, n_compared_agreement = nc_agree)
  w <- j[order(-abs(j$diff)), , drop = FALSE]
  w <- w[abs(w$diff) > tol & !is.na(w$diff), , drop = FALSE]
  rownames(w) <- NULL
  structure(list(summary = s, joined = j, worst = w, tol = tol),
            class = "cervus_comparison")
}

#' @param x A `"cervus_comparison"` object.
#' @param ... Ignored.
#' @rdname compare_cervus
#' @export
print.cervus_comparison <- function(x, ...) {
  s <- x$summary
  cat("<cervus_comparison>\n")
  cat("  rows: spawnR ", s$n_spawnR, ", CERVUS ", s$n_cervus,
      ", matched ", s$n_matched, "\n", sep = "")
  if (s$n_matched == 0L) {
    cat("  nothing matched - check that identifiers are spelled the same in both.\n")
    return(invisible(x))
  }
  cat("  LOD correlation : ", format(s$correlation, digits = 8), "\n", sep = "")
  cat("  mean |diff|     : ", format(s$mean_abs_diff, digits = 4), "\n", sep = "")
  cat("  max  |diff|     : ", format(s$max_abs_diff, digits = 4), "\n", sep = "")
  cat("  rows over tol (", x$tol, "): ", s$n_over_tol, "\n", sep = "")
  if (!is.na(s$mismatch_agreement)) {
    cat("  mismatch counts agree: ",
        format(100 * s$mismatch_agreement, digits = 4), "%\n", sep = "")
  }
  if (!is.na(s$n_compared_agreement)) {
    cat("  loci compared agree  : ",
        format(100 * s$n_compared_agreement, digits = 4), "%\n", sep = "")
  }
  if (nrow(x$worst)) {
    cat("\n  worst disagreements:\n")
    print(utils::head(x$worst[, c("offspring", "candidate", "lod_spawnR",
                                  "lod_cervus", "diff")], 6L),
          digits = 5, row.names = FALSE)
  }
  invisible(x)
}
