# Direct enumeration over all genotypes: the definitions the closed forms must
# reproduce.
enumerate_ne <- function(p) {
  k <- length(p)
  G <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  g1 <- G[, 1]; g2 <- G[, 2]
  P <- geno_freq(g1, g2, p); n <- length(P)
  ne1 <- ne2 <- nepp <- nesi <- 0
  for (o in 1:n) {
    for (c in 1:n) ne1 <- ne1 + P[o] * P[c] * (trans_prob(g1[o], g2[o], g1[c], g2[c], p) > 0)
    for (m in 1:n) {
      tom <- trans_prob(g1[o], g2[o], g1[m], g2[m], p)
      for (c in 1:n) {
        ok <- trans_prob_pair(g1[o], g2[o], g1[m], g2[m], g1[c], g2[c], p) > 0
        nepp <- nepp + P[o] * P[m] * P[c] * ok             # random trio
        ne2 <- ne2 + P[m] * tom * P[c] * ok                # true mother, random father
      }
    }
    nesi <- nesi + P[o] * (0.25 * P[o] + 0.5 * trans_prob(g1[o], g2[o], g1[o], g2[o], p) + 0.25)
  }
  c(ne_1p = ne1, ne_2p = ne2, ne_pp = nepp, ne_i = sum(P^2), ne_si = nesi)
}

one_locus <- function(p) {
  # a genotypes object whose single locus has exactly these frequencies
  wide <- data.frame(id = "x", L_a = "A1", L_b = "A1")
  g <- genotypes(wide, id_col = "id")
  g$alleles[[1]] <- paste0("A", seq_along(p))
  fr <- structure(list(L = setNames(p, g$alleles[[1]])), class = c("allele_freqs", "list"))
  marker_summary(g, freqs = fr)
}

test_that("closed forms equal direct enumeration", {
  for (p in list(c(0.5, 0.5), c(0.7, 0.2, 0.05, 0.05), rep(0.1, 10),
                 c(0.38, 0.23, 0.16, 0.12, 0.08, 0.03))) {
    ms <- one_locus(p)
    e <- enumerate_ne(p)
    for (v in names(e)) expect_equal(ms[[v]], e[[v]], tolerance = 1e-12, label = v)
  }
})

test_that("Jamieson & Taylor (1997) worked values are reproduced", {
  # Table 2, HTG4: exclusion 0.54 (eq 1a, = 1 - NE-2P), 0.36 (eq 2a, = 1 - NE-1P),
  # 0.73 (eq 3a, = 1 - NE-PP)
  ms <- one_locus(c(0.38, 0.23, 0.16, 0.12, 0.08, 0.03))
  expect_equal(round(1 - ms$ne_2p, 2), 0.54)
  expect_equal(round(1 - ms$ne_1p, 2), 0.36)
  expect_equal(round(1 - ms$ne_pp, 2), 0.73)
  # Table 1 and equations 1b, 2b, 3b: equal frequencies, n alleles
  for (n in c(2, 5, 10, 40)) {
    ms <- one_locus(rep(1 / n, n))
    expect_equal(1 - ms$ne_2p, (n - 1) * (n^3 - n^2 - 2 * n + 3) / n^4)
    expect_equal(1 - ms$ne_1p, (n - 1) * (n^2 - 3 * n + 3) / n^3)
    expect_equal(1 - ms$ne_pp, (n - 1) * (n^4 + n^3 - 7 * n^2 + 5 * n + 3) / n^5)
  }
  # Table 1 as printed, to its two decimals (0.995 for 40 alleles has three)
  tab1 <- rbind(c(2, 0.19, 0.13, 0.28), c(10, 0.79, 0.66, 0.93), c(40, 0.95, 0.90, 0.995))
  for (r in seq_len(nrow(tab1))) {
    ms <- one_locus(rep(1 / tab1[r, 1], tab1[r, 1]))
    got <- 1 - c(ms$ne_2p, ms$ne_1p, ms$ne_pp)
    expect_true(all(abs(got - tab1[r, 2:4]) <= 0.005 + 1e-12))
  }
})

test_that("heterozygosity, PIC and F_IS are computed from the reference", {
  raw <- data.frame(id = paste0("i", 1:4),
                    L1a = c("A", "A", "B", "A"), L1b = c("B", "A", "B", "A"),
                    L2a = c("C", "C", "D", NA),  L2b = c("D", "E", "E", NA))
  g <- genotypes(raw, id_col = "id")
  ms <- marker_summary(g)
  expect_equal(ms$n_typed, c(4L, 3L))
  expect_equal(ms$h_obs, c(1 / 4, 1))
  p1 <- c(5, 3) / 8
  expect_equal(ms$h_exp[1], 8 / 7 * (1 - sum(p1^2)))
  expect_equal(ms$pic[1], 1 - sum(p1^2) - 2 * p1[1]^2 * p1[2]^2)
  expect_equal(ms$f_is[1], 1 - ms$h_obs[1] / ms$h_exp[1])
  expect_lt(ms$f_is[2], 0)                       # all heterozygous: excess
  # no frequency floor is applied: an allele seen in nobody counts for nothing
  sub <- marker_summary(g, reference = c("i2", "i4"))
  expect_equal(sub$n_alleles[1], 1L)
  expect_equal(sub$ne_1p[1], 1)
})

test_that("combined values multiply over loci", {
  sim <- simulate_population(n = 40, n_loci = 6, n_alleles = 5, seed = 3)
  ms <- marker_summary(sim$genotypes, reference = c(sim$dams, sim$sires))
  cmb <- attr(ms, "combined")
  expect_equal(cmb[["ne_pp"]], prod(ms$ne_pp))
  expect_equal(cmb[["ne_1p"]], prod(ms$ne_1p))
  expect_equal(cmb[["mean_h_exp"]], mean(ms$h_exp))
  expect_true(all(ms$ne_pp <= ms$ne_2p & ms$ne_2p <= ms$ne_1p))
  expect_output(print(ms), "NE-PP")
  fr <- suppressMessages(allele_freqs(sim$genotypes))
  expect_equal(nrow(marker_summary(sim$genotypes, freqs = fr)), 6L)
})
