cross_fixture <- function() {
  sim <- simulate_population(n = 40, n_loci = 12, n_alleles = 6,
                             n_offspring = 30, seed = 11)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  # each dam crossed with two sires, plus every true cross so the list is complete
  cr <- unique(rbind(
    data.frame(mother = rep(sim$dams[1:10], each = 2),
               father = sim$sires[c(1:10, 2:10, 1)]),
    data.frame(mother = sim$pedigree$dam, father = sim$pedigree$sire)))
  list(sim = sim, p = p, cr = cr)
}

test_that("crosses are validated", {
  fx <- cross_fixture()
  g <- fx$sim$genotypes
  expect_error(assign_trios(g, fx$sim$offspring[1], freqs = fx$p, error = 0.01,
                            crosses = data.frame(mother = "nobody", father = fx$sim$sires[1])),
               "not present")
  expect_error(assign_trios(g, fx$sim$offspring[1], mothers = fx$sim$dams[1],
                            fathers = fx$sim$sires, freqs = fx$p, error = 0.01,
                            crosses = fx$cr[fx$cr$mother != fx$sim$dams[1], ][1, ]),
               "must be among")
  expect_error(assign_trios(g, fx$sim$offspring[1], freqs = fx$p, error = 0.01),
               "Give `mothers`")
  cr2 <- fx$cr; cr2$weight <- -1
  expect_error(assign_trios(g, fx$sim$offspring[1], freqs = fx$p, error = 0.01,
                            crosses = cr2), "non-negative")
  p <- fx$p
  expect_error(sim_delta(p, 0.01, n_candidates = 3, n_candidates2 = 2, nsim = 10,
                         type = "pair", crosses = cbind(4, 1)), "must index")
  expect_error(sim_delta(p, 0.01, n_candidates = 3, nsim = 10, crosses = cbind(1, 1)),
               "only to type = 'pair'")
})

test_that("only listed crosses are scored, with the pair LODs of parentage_lod()", {
  fx <- cross_fixture()
  set.seed(1)
  a <- assign_trios(fx$sim$genotypes, fx$sim$offspring, freqs = fx$p,
                    error = 0.01, crosses = fx$cr, nsim = 300, nsim_trio = 200)
  key <- paste(fx$cr$mother, fx$cr$father)
  expect_true(all(paste(a$trio_mother, a$trio_father) %in% key))
  ped <- fx$sim$pedigree[match(a$offspring, fx$sim$pedigree$offspring), ]
  expect_gte(mean(a$trio_mother == ped$dam & a$trio_father == ped$sire), 0.9)
  tab <- parentage_lod(fx$sim$genotypes, fx$sim$offspring,
                       fx$cr[c("father", "mother")], fx$p, 0.01, type = "pair")
  ds <- delta_stat(tab)
  expect_equal(ds$lod, a$trio_lod)
  expect_equal(ds$delta, a$trio_delta)
})

test_that("cross posteriors are proper and the set covers the best cross", {
  fx <- cross_fixture()
  set.seed(2)
  a <- assign_trios(fx$sim$genotypes, fx$sim$offspring, freqs = fx$p,
                    error = 0.01, crosses = fx$cr, prop_sampled_crosses = 0.9,
                    nsim = 300, nsim_trio = 200, output = "all")
  tot <- tapply(a$posterior, a$offspring, sum)
  b <- assign_trios(fx$sim$genotypes, fx$sim$offspring, freqs = fx$p,
                    error = 0.01, crosses = fx$cr, prop_sampled_crosses = 0.9,
                    nsim = 300, nsim_trio = 200)
  expect_equal(as.vector(tot[b$offspring]) + b$trio_post_none, rep(1, nrow(b)))
  expect_true(all(b$trio_post_none > 0))
  first <- vapply(strsplit(b$trio_set, "; "), `[`, "", 1)
  expect_identical(first, paste(b$trio_mother, b$trio_father, sep = " x "))
  expect_true(all(b$trio_set_size >= 1L))
  fam <- attr(b, "families")
  expect_equal(nrow(fam), nrow(fx$cr))
  expect_equal(sum(fam$expected) + sum(b$trio_post_none), nrow(b))
  expect_equal(sum(fam$assigned), sum(b$trio_significant))
})

test_that("cross weights act as a prior", {
  fx <- cross_fixture()
  L <- matrix(c(2, 2), 1)                       # two crosses, equal evidence
  even <- spawnR:::.cross_posterior(L, c(1, 1), 1)
  skew <- spawnR:::.cross_posterior(L, c(3, 1), 1)
  expect_equal(as.vector(even$P), c(0.5, 0.5))
  expect_equal(as.vector(skew$P), c(0.75, 0.25))
  expect_equal(spawnR:::.cross_posterior(matrix(NA_real_, 1, 2), c(1, 1), 1)$none, NA_real_)
})

test_that("the cross-design simulation draws the true cross from the list", {
  sim <- simulate_population(n = 40, n_loci = 12, n_alleles = 6, seed = 12)
  p <- suppressMessages(allele_freqs(sim$genotypes))
  # factorial design given as crosses behaves like the all-pairs simulation
  set.seed(3)
  a <- sim_delta(p, 0.01, n_candidates = 4, n_candidates2 = 3, nsim = 800, type = "pair")
  set.seed(4)
  b <- sim_delta(p, 0.01, n_candidates = 4, n_candidates2 = 3, nsim = 800, type = "pair",
                 crosses = expand.grid(f = 1:4, m = 1:3))
  expect_lt(abs(mean(a$correct) - mean(b$correct)), 0.05)
  expect_equal(attr(b, "settings")$n_crosses, 12L)
  # nothing is correct when offspring never come from a listed cross
  d <- sim_delta(p, 0.01, n_candidates = 4, n_candidates2 = 3, nsim = 200, type = "pair",
                 crosses = cbind(1:3, 1:3), prop_sampled = 0)
  expect_false(any(d$correct))
})
