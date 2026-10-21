test_that("exclusions at zero error count as -700 in simulations, as in scoring", {
  expect_equal(spawnR:::.lod_for_sum(c(1.5, NA, -Inf, Inf)), c(1.5, 0, -700, 700))
  sim <- simulate_population(n = 30, n_loci = 20, n_alleles = 4, seed = 4)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  og <- spawnR:::.row_of(sim$genotypes, sim$offspring[1])
  set.seed(1)
  v <- spawnR:::.sim_unrelated(og, p, rep(0, 20), 500)
  # an unrelated individual is almost always excluded somewhere at 20 loci;
  # before the fix those loci silently contributed 0
  expect_gt(mean(v < -600), 0.5)
})

test_that("max_mismatch inside sim_delta removes true parents carrying errors", {
  sim <- simulate_population(n = 40, n_loci = 30, n_alleles = 6, seed = 5)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  set.seed(2)
  open <- sim_delta(p, error = 0.05, n_candidates = 10, nsim = 1500)
  set.seed(2)
  strict <- sim_delta(p, error = 0.05, n_candidates = 10, nsim = 1500,
                      max_mismatch = 0)
  # with 5% error over 30 loci most true parents mismatch somewhere
  expect_lt(mean(strict$correct), mean(open$correct) - 0.2)
  set.seed(3)
  pr <- sim_delta(p, error = 0.05, n_candidates = 4, n_candidates2 = 3,
                  nsim = 400, type = "pair", max_mismatch = 0)
  set.seed(3)
  pr0 <- sim_delta(p, error = 0.05, n_candidates = 4, n_candidates2 = 3,
                   nsim = 400, type = "pair")
  expect_lt(mean(pr$correct), mean(pr0$correct) - 0.2)
})

test_that("min_loci inside sim_delta drops comparisons with too few loci", {
  sim <- simulate_population(n = 40, n_loci = 10, n_alleles = 6, seed = 6)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  s <- sim_delta(p, error = 0.01, n_candidates = 5, nsim = 200, min_loci = 11)
  expect_true(all(is.na(s$delta)) && !any(s$correct))
  s2 <- sim_delta(p, error = 0.01, n_candidates = 5, nsim = 200,
                  prop_typed = 0.6, min_loci = 4)
  expect_true(any(!is.na(s2$delta)))
  sp <- sim_delta(p, error = 0.01, n_candidates = 3, n_candidates2 = 2,
                  nsim = 100, type = "pair", min_loci = 11)
  expect_true(all(is.na(sp$delta)))
})

test_that("min_loci screens real comparisons", {
  sim <- simulate_population(n = 30, n_loci = 10, n_alleles = 6,
                             n_offspring = 5, seed = 7)
  g <- sim$genotypes
  o1 <- match(sim$offspring[1], g$ids)
  g$a1[o1, 1:7] <- NA; g$a2[o1, 1:7] <- NA      # offspring 1 typed at 3 loci
  p <- suppressMessages(allele_freqs(g))
  tab <- parentage_lod(g, sim$offspring, sim$sires, p, 0.01, min_loci = 5)
  expect_false(sim$offspring[1] %in% tab$offspring)
  expect_true(all(tab$n_compared >= 5))
  set.seed(4)
  a <- assign_trios(g, sim$offspring, sim$dams, sim$sires, p, 0.01,
                    nsim = 200, nsim_trio = 50, min_loci = 5)
  r <- a[a$offspring == sim$offspring[1], ]
  expect_true(is.na(r$mother) && is.na(r$father) && is.na(r$trio_mother))
  expect_false(r$trio_significant)
  expect_true(all(!is.na(a$trio_mother[-1])))
})

test_that("prop_typed defaults to the individuals analysed", {
  sim <- simulate_population(n = 30, n_loci = 10, n_alleles = 6,
                             n_offspring = 5, seed = 8)
  g <- sim$genotypes
  ids <- c(sim$offspring, sim$dams, sim$sires)
  expect_equal(spawnR:::.prop_typed(g, ids), 1)
  # an untyped bystander in `g` no longer drags the simulated typing rate down
  extra <- sim$offspring[5]
  g$a1[match(extra, g$ids), ] <- NA; g$a2[match(extra, g$ids), ] <- NA
  expect_equal(spawnR:::.prop_typed(g, setdiff(ids, extra)), 1)
  expect_lt(mean(!is.na(g$a1)), 1)
})
