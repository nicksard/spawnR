test_that("a small reference set triggers the small-sample message", {
  sim <- simulate_population(n = 40, n_loci = 5, n_alleles = 4, seed = 8)
  six <- c(sim$sires[1:3], sim$dams[1:3])
  expect_message(allele_freqs(sim$genotypes, reference = six), "fewer than 25")
  expect_silent(allele_freqs(sim$genotypes, reference = six, min_freq = 0.01))
  expect_silent(allele_freqs(sim$genotypes))   # 60 typed individuals
})

test_that("the default floor is 1/(2n) capped at 0.01", {
  sim <- simulate_population(n = 40, n_loci = 5, n_alleles = 4, seed = 8)
  six <- c(sim$sires[1:3], sim$dams[1:3])
  p1 <- suppressMessages(allele_freqs(sim$genotypes, reference = six))
  p2 <- allele_freqs(sim$genotypes, reference = six, min_freq = 0.01)
  expect_equal(unclass(p1), unclass(p2))
  # with many typed individuals 1/(2n) is below the cap and is used as before
  big <- simulate_population(n = 400, n_loci = 3, n_alleles = 4, seed = 9)
  ntyped <- colSums(!is.na(big$genotypes$a1))
  expect_equal(unclass(allele_freqs(big$genotypes)),
               unclass(allele_freqs(big$genotypes, min_freq = 1 / (2 * ntyped[1]))))
})
