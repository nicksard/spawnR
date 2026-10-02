#' Genotype containers for parentage analysis
#'
#' `genotypes()` converts a wide-format table of codominant marker data into
#' the compact integer representation used throughout \pkg{spawnR}. The
#' input has one row per individual, an identifier column, and then two adjacent
#' columns per locus holding the two observed alleles.
#'
#' Alleles are recoded to consecutive integers within each locus, and the
#' original labels are retained so that results can be reported on the original
#' scale. Any allele matching `missing` (and any `NA`) is treated as untyped; if
#' either allele of a locus is untyped the whole locus is set to untyped for that
#' individual, since a half-scored codominant genotype carries no usable
#' information for a likelihood-ratio calculation.
#'
#' @param x A `data.frame` or `matrix` with one row per individual: an
#'   identifier column followed by `2 * n_loci` allele columns in locus order.
#' @param id_col Name or index of the identifier column. If `NULL`, row names
#'   are used, or sequential identifiers are generated.
#' @param loci Optional character vector of locus names, length `n_loci`. If
#'   `NULL`, names are taken from the first column of each pair with a trailing
#'   separator plus suffix removed, falling back to `"L1"`, `"L2"`, ...
#' @param missing Character vector of allele codes to treat as untyped, in
#'   addition to `NA`.
#'
#' @return An object of class `"genotypes"`: a list with elements `ids`
#'   (character), `loci` (character), `alleles` (a list of character vectors of
#'   allele labels, one per locus), and `a1`, `a2` (integer matrices of
#'   dimension `n_ind` by `n_loci` holding recoded alleles, `NA` where untyped).
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' @examples
#' raw <- data.frame(
#'   id  = c("a", "b", "c"),
#'   L1a = c(101, 103, 101), L1b = c(103, 103, 105),
#'   L2a = c(201, 201, NA),  L2b = c(205, 203, NA)
#' )
#' g <- genotypes(raw, id_col = "id")
#' g
#' n_ind(g)
#' n_loci(g)
#' @export
genotypes <- function(x, id_col = 1, loci = NULL, missing = c("", "0", "000", "NA", "*", "-9")) {
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.data.frame(x)) stop("`x` must be a data.frame or a matrix.", call. = FALSE)
  if (nrow(x) < 1L) stop("`x` must have at least one row.", call. = FALSE)

  if (is.null(id_col)) {
    ids <- if (!is.null(rownames(x))) rownames(x) else paste0("ind", seq_len(nrow(x)))
    allele_cols <- seq_len(ncol(x))
  } else {
    j <- if (is.character(id_col)) match(id_col, names(x)) else as.integer(id_col)
    if (is.na(j) || j < 1L || j > ncol(x)) stop("`id_col` does not identify a column of `x`.", call. = FALSE)
    ids <- as.character(x[[j]])
    allele_cols <- setdiff(seq_len(ncol(x)), j)
  }
  if (anyDuplicated(ids)) stop("Individual identifiers must be unique.", call. = FALSE)
  if (length(allele_cols) %% 2L != 0L) {
    stop("Expected an even number of allele columns (two per locus); got ",
         length(allele_cols), ".", call. = FALSE)
  }
  nl <- length(allele_cols) %/% 2L
  if (nl < 1L) stop("At least one locus is required.", call. = FALSE)

  odd  <- allele_cols[seq(1L, by = 2L, length.out = nl)]
  even <- allele_cols[seq(2L, by = 2L, length.out = nl)]

  if (is.null(loci)) {
    nm <- names(x)[odd]
    nm <- sub("[._ -]?[abAB12]$", "", nm)
    bad <- is.na(nm) | !nzchar(nm) | duplicated(nm)
    nm[bad] <- paste0("L", which(bad))
    loci <- nm
  }
  loci <- as.character(loci)
  if (length(loci) != nl) stop("`loci` must have length ", nl, ".", call. = FALSE)

  n <- nrow(x)
  a1 <- matrix(NA_integer_, n, nl, dimnames = list(ids, loci))
  a2 <- a1
  allele_labels <- vector("list", nl)
  names(allele_labels) <- loci

  for (l in seq_len(nl)) {
    v1 <- .clean_allele(x[[odd[l]]],  missing)
    v2 <- .clean_allele(x[[even[l]]], missing)
    untyped <- is.na(v1) | is.na(v2)
    v1[untyped] <- NA_character_
    v2[untyped] <- NA_character_
    lab <- sort(unique(c(v1, v2)), na.last = NA)
    allele_labels[[l]] <- lab
    a1[, l] <- match(v1, lab)
    a2[, l] <- match(v2, lab)
  }

  structure(list(ids = ids, loci = loci, alleles = allele_labels, a1 = a1, a2 = a2),
            class = "genotypes")
}

.clean_allele <- function(v, missing) {
  v <- as.character(v)
  v <- trimws(v)
  v[v %in% missing] <- NA_character_
  v
}

#' @param x A `"genotypes"` object.
#' @rdname genotypes
#' @return `n_ind()` and `n_loci()` return a single integer.
#' @export
n_ind <- function(x) {
  stopifnot(inherits(x, "genotypes"))
  length(x$ids)
}

#' @rdname genotypes
#' @export
n_loci <- function(x) {
  stopifnot(inherits(x, "genotypes"))
  length(x$loci)
}

#' @param i Individuals to keep: a character vector of identifiers, an integer
#'   vector of positions, or a logical vector.
#' @param ... Ignored.
#' @rdname genotypes
#' @return `subset_ind()` returns a `"genotypes"` object.
#' @export
subset_ind <- function(x, i) {
  stopifnot(inherits(x, "genotypes"))
  if (is.character(i)) {
    k <- match(i, x$ids)
    if (anyNA(k)) stop("Unknown identifier(s): ", paste(i[is.na(k)], collapse = ", "), call. = FALSE)
    i <- k
  }
  x$ids <- x$ids[i]
  x$a1  <- x$a1[i, , drop = FALSE]
  x$a2  <- x$a2[i, , drop = FALSE]
  x
}

#' @rdname genotypes
#' @export
print.genotypes <- function(x, ...) {
  cat("<genotypes>", n_ind(x), "individuals,", n_loci(x), "loci\n")
  na <- vapply(x$alleles, length, integer(1))
  cat("  alleles per locus: ", paste0(range(na), collapse = "-"),
      " (mean ", format(mean(na), digits = 3), ")\n", sep = "")
  miss <- mean(is.na(x$a1))
  cat("  missing genotypes: ", format(100 * miss, digits = 3), "%\n", sep = "")
  cat("  first ids: ", paste(utils::head(x$ids, 5L), collapse = ", "),
      if (n_ind(x) > 5L) ", ..." else "", "\n", sep = "")
  invisible(x)
}

.row_of <- function(g, id, what = "individual") {
  if (inherits(id, "genotypes")) {
    if (n_ind(id) != 1L) stop("Expected exactly one ", what, ".", call. = FALSE)
    return(list(a1 = id$a1[1L, ], a2 = id$a2[1L, ]))
  }
  i <- if (is.character(id)) match(id, g$ids) else as.integer(id)
  if (length(i) != 1L || is.na(i)) stop("Cannot locate ", what, ": ", paste(id, collapse = ", "), call. = FALSE)
  list(a1 = g$a1[i, ], a2 = g$a2[i, ])
}
