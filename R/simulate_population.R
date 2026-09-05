#' Simulate a small population with a known pedigree
#'
#' Generates a population of unrelated adults under Hardy-Weinberg proportions,
#' produces offspring by random pairing, and then applies genotyping error and
#' missing data. Intended for examples, regression tests and calibration of the
#' simulation-based confidence criteria; it is not a demographic simulator.
#'
#' Genotyping error is applied with the random-genotype-replacement model used
#' throughout the package: with probability `error` an individual's genotype at a
#' locus is discarded and replaced by a genotype drawn under Hardy-Weinberg
#' proportions. Amiri Roudbar et al. (2025) draw the per-locus rate from a range
#' rather than holding it constant, which `error` reproduces when given a vector.
#'
#' @param n Number of adults, split as evenly as possible into sires and dams.
#' @param n_loci Number of loci.
#' @param n_alleles Number of alleles per locus.
#' @param n_offspring Number of offspring. Defaults to `n / 2`.
#' @param error Genotyping error rate: one value, or one per locus.
#' @param prop_missing Probability that a genotype is untyped.
#' @param freq_decay Controls the allele frequency spectrum. Allele `i` has
#'   expected frequency proportional to `freq_decay^(i - 1)` before
#'   normalisation, so `1` gives even frequencies and smaller values give a few
#'   common and many rare alleles.
#' @param seed Optional integer seed.
#'
#' @return A list with elements `genotypes` (a `"genotypes"` object holding
#'   adults and offspring), `pedigree` (a `data.frame` with columns `offspring`,
#'   `sire`, `dam`), `sires`, `dams`, `offspring` (identifier vectors),
#'   `true_freqs` (the allele frequencies used to generate the adults) and
#'   `error` (the per-locus rates used).
#'
#' @references
#' Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
#' (2025) Pairwise paternity assignment with forward-backward simulations:
#' refining CERVUS using trio-based likelihood and locus-specific error rates.
#' \emph{Ecology and Evolution}, \strong{15}(10), e72230.
#' \doi{10.1002/ece3.72230}
#'
#' @examples
#' sim <- simulate_population(n = 20, n_loci = 5, n_alleles = 6, seed = 7)
#' sim$genotypes
#' head(sim$pedigree)
#' @export
simulate_population <- function(n = 100, n_loci = 12, n_alleles = 8,
                                n_offspring = NULL, error = 0.01,
                                prop_missing = 0, freq_decay = 0.75,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(n); n_loci <- as.integer(n_loci); n_alleles <- as.integer(n_alleles)
  if (n < 2L || n_loci < 1L || n_alleles < 2L) {
    stop("Need n >= 2, n_loci >= 1 and n_alleles >= 2.", call. = FALSE)
  }
  if (is.null(n_offspring)) n_offspring <- n %/% 2L
  n_offspring <- as.integer(n_offspring)
  e <- .expand_error(error, n_loci)

  base <- freq_decay^(seq_len(n_alleles) - 1)
  freqs <- lapply(seq_len(n_loci), function(l) {
    p <- base * stats::runif(n_alleles, 0.5, 1.5)
    p <- p / sum(p)
    names(p) <- paste0("A", seq_len(n_alleles))
    p
  })

  n_sire <- n %/% 2L
  n_dam <- n - n_sire
  sires <- sprintf("S%03d", seq_len(n_sire))
  dams  <- sprintf("D%03d", seq_len(n_dam))
  offs  <- sprintf("O%03d", seq_len(n_offspring))

  A1 <- matrix(NA_integer_, n + n_offspring, n_loci)
  A2 <- A1
  ad <- seq_len(n)
  for (l in seq_len(n_loci)) {
    A1[ad, l] <- sample.int(n_alleles, n, TRUE, prob = freqs[[l]])
    A2[ad, l] <- sample.int(n_alleles, n, TRUE, prob = freqs[[l]])
  }

  sire_of <- sample.int(n_sire, n_offspring, TRUE)
  dam_of  <- n_sire + sample.int(n_dam, n_offspring, TRUE)
  for (l in seq_len(n_loci)) {
    from_s <- ifelse(stats::runif(n_offspring) < 0.5, A1[sire_of, l], A2[sire_of, l])
    from_d <- ifelse(stats::runif(n_offspring) < 0.5, A1[dam_of, l],  A2[dam_of, l])
    A1[n + seq_len(n_offspring), l] <- from_s
    A2[n + seq_len(n_offspring), l] <- from_d
  }

  N <- n + n_offspring
  for (l in seq_len(n_loci)) {
    hit <- stats::runif(N) < e[l]
    if (any(hit)) {
      m <- sum(hit)
      A1[hit, l] <- sample.int(n_alleles, m, TRUE, prob = freqs[[l]])
      A2[hit, l] <- sample.int(n_alleles, m, TRUE, prob = freqs[[l]])
    }
    if (prop_missing > 0) {
      gone <- stats::runif(N) < prop_missing
      A1[gone, l] <- NA_integer_
      A2[gone, l] <- NA_integer_
    }
  }

  ids <- c(sires, dams, offs)
  lab <- paste0("A", seq_len(n_alleles))
  wide <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (l in seq_len(n_loci)) {
    wide[[paste0("Loc", l, "_a")]] <- lab[A1[, l]]
    wide[[paste0("Loc", l, "_b")]] <- lab[A2[, l]]
  }
  g <- genotypes(wide, id_col = "id")

  ped <- data.frame(offspring = offs, sire = sires[sire_of],
                    dam = dams[dam_of - n_sire], stringsAsFactors = FALSE)
  names(freqs) <- g$loci
  list(genotypes = g, pedigree = ped, sires = sires, dams = dams,
       offspring = offs, true_freqs = freqs, error = e)
}
