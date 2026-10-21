# A delta_sim built by hand, so the expected critical values can be worked out.
fake_sim <- function(delta, correct) {
  structure(data.frame(delta = delta, lod = delta, correct = correct,
                       sampled = TRUE),
            class = c("delta_sim", "data.frame"))
}

test_that("an all-correct simulation gives critical Delta 0, unresolved", {
  x <- fake_sim(seq(5, 50, length.out = 100), rep(TRUE, 100))
  dc <- delta_critical(x, levels = c(0.80, 0.95, 0.99))
  expect_equal(dc$critical_delta, c(0, 0, 0))
  expect_identical(dc$resolved, c(FALSE, FALSE, FALSE))
  expect_equal(dc$n_above, rep(100L, 3))
  expect_equal(dc$prop_correct, rep(1, 3))
})

test_that("wrong assignments at low Delta put the threshold above them", {
  # 80 correct at Delta 11..90, then 20 wrong at Delta 1.5..20.5 interleaved
  # with the bottom of the range: proportion correct overall is 0.8. Simulated
  # Deltas are continuous, so the values are kept free of ties.
  delta <- c(11:90, 1:20 + 0.5)
  correct <- c(rep(TRUE, 80), rep(FALSE, 20))
  x <- fake_sim(delta, correct)
  dc <- delta_critical(x, levels = c(0.95, 0.99), min_n = 20L)
  expect_true(all(dc$resolved))
  expect_true(all(dc$critical_delta > 0))
  # every wrong assignment with Delta above the threshold must still leave the
  # proportion correct at or above the level
  for (j in seq_len(nrow(dc))) {
    above <- delta >= dc$critical_delta[j]
    expect_gte(mean(correct[above]), dc$confidence[j])
    expect_equal(sum(above), dc$n_above[j])
  }
  # and the threshold is the most permissive one: one step lower fails
  o <- sort(delta, decreasing = TRUE)
  k <- dc$n_above[1]
  expect_lt(mean(correct[delta >= o[k + 1]]), dc$confidence[1])
})

test_that("a level met by the whole simulation, errors included, gives 0", {
  x <- fake_sim(c(11:100, 1:10 + 0.5), c(rep(TRUE, 90), rep(FALSE, 10)))
  dc <- delta_critical(x, levels = c(0.80, 0.95))
  expect_equal(dc$critical_delta[1], 0)
  expect_false(dc$resolved[1])
  expect_equal(dc$prop_correct[1], 0.9)
  expect_true(dc$resolved[2])
  expect_gt(dc$critical_delta[2], 0)
})

test_that("a level never reached gives NA", {
  x <- fake_sim(1:100, rep(c(TRUE, FALSE), 50))
  dc <- delta_critical(x, levels = 0.95)
  expect_true(is.na(dc$critical_delta))
  expect_true(is.na(dc$resolved))
  expect_equal(dc$n_above, 0L)
  dc0 <- delta_critical(fake_sim(rep(NA_real_, 5), rep(FALSE, 5)))
  expect_true(all(is.na(dc0$critical_delta)) && all(is.na(dc0$resolved)))
})

test_that("min_n still applies before the whole-simulation case", {
  x <- fake_sim(1:10, rep(TRUE, 10))
  expect_true(is.na(delta_critical(x, levels = 0.8, min_n = 20L)$critical_delta))
})

test_that("a powerful SNP panel gives critical Delta 0 (reported reproduction)", {
  set.seed(1)
  pop <- simulate_population(n = 60, n_loci = 150, n_alleles = 2, seed = 1)
  p <- allele_freqs(pop$genotypes)
  s <- sim_delta(p, error = 0.01, n_candidates = 4, nsim = 2000)
  expect_equal(mean(s$correct), 1)
  dc <- delta_critical(s)
  expect_equal(dc$critical_delta, c(0, 0))
  expect_identical(dc$resolved, c(FALSE, FALSE))
})

test_that("a threshold of 0 passes every positive Delta and no NA Delta", {
  crit <- data.frame(confidence = c(0.8, 0.95), critical_delta = c(0, 0))
  lv <- spawnR:::.conf_level(c(0, 0.001, 12, NA), crit)
  expect_equal(lv, c(0.95, 0.95, 0.95, NA))
})

test_that("assign_trios calls every true trio significant on a strong panel", {
  sim <- simulate_population(n = 12, n_loci = 150, n_alleles = 2,
                             n_offspring = 10, seed = 2)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  set.seed(3)
  a <- assign_trios(sim$genotypes, sim$offspring, sim$dams, sim$sires, p,
                    error = 0.01, nsim = 500, nsim_trio = 300)
  ped <- sim$pedigree[match(a$offspring, sim$pedigree$offspring), ]
  right <- a$trio_mother == ped$dam & a$trio_father == ped$sire
  expect_true(all(right))
  expect_true(all(a$trio_significant[right]))
  crit <- attr(a, "criteria")
  expect_true("resolved" %in% names(crit))
})
