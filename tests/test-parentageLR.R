library(parentageLR)

ok <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAILED: ", label, call. = FALSE)
  cat("ok  ", label, "\n", sep = "")
}
near <- function(a, b, tol = 1e-10) all(abs(a - b) < tol)

## ---------------------------------------------------------------- genotypes
raw <- data.frame(
  id  = c("a", "b", "c", "d"),
  L1a = c("101", "103", "101", NA),   L1b = c("103", "103", "105", "101"),
  L2a = c("201", "201", "0",   "203"), L2b = c("205", "203", "205", "203"),
  stringsAsFactors = FALSE
)
g <- genotypes(raw, id_col = "id")
ok(n_ind(g) == 4L && n_loci(g) == 2L, "genotypes() shape")
ok(identical(g$loci, c("L1", "L2")), "locus names stripped from column names")
ok(identical(g$alleles[[1]], c("101", "103", "105")), "allele labels recoded and sorted")
ok(is.na(g$a1[4, 1]) && is.na(g$a2[4, 1]), "half-missing genotype is fully untyped")
ok(is.na(g$a1[3, 2]), "'0' is treated as a missing allele code")
ok(inherits(try(genotypes(raw[, 1:2], id_col = "id"), silent = TRUE), "try-error"),
   "odd number of allele columns is rejected")
ok(n_ind(subset_ind(g, c("a", "c"))) == 2L, "subset_ind() by identifier")

p <- allele_freqs(g)
ok(all(vapply(p, function(x) abs(sum(x) - 1) < 1e-12, logical(1))), "allele_freqs() sum to one")
ok(all(unlist(p) > 0), "allele_freqs() are strictly positive")

## ------------------------------------------------- Mendelian identities
pf <- c(0.4, 0.25, 0.2, 0.1, 0.05)
G <- unique(t(apply(t(utils::combn(rep(seq_along(pf), 2), 2)), 1, sort)))
o1 <- G[, 1]; o2 <- G[, 2]
Pg <- geno_freq(o1, o2, pf)
ok(near(sum(Pg), 1), "Hardy-Weinberg genotype frequencies sum to one")
for (r in seq_len(nrow(G))) {
  ok0 <- near(sum(trans_prob(o1, o2, G[r, 1], G[r, 2], pf)), 1)
  if (!ok0) stop("T(o|f) does not sum to one for parent row ", r)
}
ok(TRUE, "sum_o T(o | f) = 1 for every parental genotype")
ok(near(sum(trans_prob_pair(o1, o2, G[3, 1], G[3, 2], G[7, 1], G[7, 2], pf)), 1),
   "sum_o T(o | m, f) = 1")
oo <- G[5, ]
ok(near(sum(Pg * trans_prob_pair(oo[1], oo[2], G[3, 1], G[3, 2], G[, 1], G[, 2], pf)),
        trans_prob(oo[1], oo[2], G[3, 1], G[3, 2], pf)),
   "identity I2: marginalising the second parent gives T(o | m)")
i3 <- sum(vapply(seq_len(nrow(G)), function(i)
  Pg[i] * sum(Pg * trans_prob_pair(oo[1], oo[2], G[i, 1], G[i, 2], G[, 1], G[, 2], pf)),
  numeric(1)))
ok(near(i3, geno_freq(oo[1], oo[2], pf)),
   "identity I3: offspring of two random parents is a random individual")

ok(near(trans_prob(1L, 2L, 1L, 2L, pf), (pf[1] + pf[2]) / 2), "T table: AiAj parent, AiAj offspring")
ok(near(trans_prob(1L, 1L, 1L, 2L, pf), pf[1] / 2),           "T table: AiAj parent, AiAi offspring")
ok(near(trans_prob(1L, 3L, 1L, 2L, pf), pf[3] / 2),           "T table: AiAj parent, AiAk offspring")
ok(near(trans_prob(2L, 3L, 1L, 1L, pf), 0),                   "T table: no shared allele gives zero")
ok(near(trans_prob_pair(1L, 2L, 1L, 1L, 2L, 2L, pf), 1),      "T2: fully determined cross")
ok(is.na(trans_prob(NA_integer_, 2L, 1L, 2L, pf)),            "missing alleles propagate to NA")

## ----------------------------- the error expansion is a probability model
lik <- parentageLR:::.locus_lik
for (e in c(0, 1e-6, 0.001, 0.01, 0.05, 0.2, 0.5, 0.9, 1)) {
  for (ty in c("one_known", "both_unknown", "pair")) {
    L <- lik(o1, o2, 2L, 4L, pf, e, 1L, 3L, ty)
    if (!near(sum(L$num), 1, 1e-9)) stop("L(H1) does not sum to one: e=", e, " type=", ty)
    if (!near(sum(L$den), 1, 1e-9)) stop("L(H2) does not sum to one: e=", e, " type=", ty)
  }
}
ok(TRUE, "L(H1) and L(H2) are proper distributions at every error rate and type")

# The coefficients printed by Amiri Roudbar et al. (2025, eq. 1) do not have
# this property; this test records the discrepancy that motivated re-deriving
# them from the error model. If a future correction is published, revisit
# inst/notes/derivation.md before changing lod_locus().
e <- 0.2
Po <- geno_freq(o1, o2, pf)
Tk <- trans_prob(o1, o2, 1L, 3L, pf)
Tc <- trans_prob(o1, o2, 2L, 4L, pf)
Tkc <- trans_prob_pair(o1, o2, 1L, 3L, 2L, 4L, pf)
as_printed <- (1 - e)^3 * Tkc + e * (1 - e)^2 * (Tk + Tc + Po) +
  e^2 * (1 - e)^3 * Po + e^3 * Po
ok(abs(sum(as_printed) - 1) > 0.05,
   "published coefficients do not sum to one (documented departure)")

## ---------------------------------------------------------- LOD behaviour
ok(all(lod_locus(o1, o2, 2L, 4L, pf, 1, 1L, 3L, "one_known") == 0),
   "an error rate of one makes every locus uninformative")
a <- lod_locus(o1, o2, 2L, 4L, pf, 0, type = "both_unknown")
b <- log(trans_prob(o1, o2, 2L, 4L, pf) / Po)
fin <- is.finite(a) & is.finite(b)
ok(near(a[fin], b[fin]), "zero error reduces to the classical Marshall ratio (both unknown)")
a <- lod_locus(o1, o2, 2L, 4L, pf, 0, 1L, 3L, "one_known")
b <- log(trans_prob_pair(o1, o2, 1L, 3L, 2L, 4L, pf) / trans_prob(o1, o2, 1L, 3L, pf))
fin <- is.finite(a) & is.finite(b)
ok(near(a[fin], b[fin]), "zero error reduces to the classical ratio (one parent known)")
ok(all(is.finite(lod_locus(o1, o2, 2L, 4L, pf, 0.01, type = "both_unknown"))),
   "a positive error rate keeps every locus finite")
ok(is.na(lod_locus(NA_integer_, NA_integer_, 2L, 4L, pf, 0.01, type = "both_unknown")),
   "a missing offspring genotype gives NA")

# Regression test: a length-one denominator must recycle against a longer
# numerator. Getting this wrong silently returned NA for all but the first
# candidate.
v <- lod_locus(1L, 2L, c(1L, 1L, 3L, 4L), c(2L, 3L, 4L, 5L), pf, 0.01,
               type = "both_unknown")
ok(length(v) == 4L && !anyNA(v), "scalar denominator recycles across many candidates")

## --------------------------------------------------------- parentage_lod
sim <- simulate_population(n = 60, n_loci = 15, n_alleles = 8,
                           n_offspring = 40, error = 0.01, seed = 4321)
gg <- sim$genotypes
pp <- allele_freqs(gg)
ped <- sim$pedigree

tab <- parentage_lod(gg, ped$offspring, sim$sires, pp, error = 0.01,
                     type = "both_unknown")
ok(nrow(tab) == nrow(ped) * length(sim$sires), "one row per offspring by candidate")
ok(all(c("offspring", "candidate", "n_compared", "mismatches", "lod") %in% names(tab)),
   "parentage_lod() column names")
d <- delta_stat(tab)
truth <- ped$sire[match(d$offspring, ped$offspring)]
ok(mean(d$candidate == truth) >= 0.95,
   "top candidate is the true sire for at least 95% of offspring")
ok(all(d$delta[!is.na(d$lod2)] == d$lod[!is.na(d$lod2)] - d$lod2[!is.na(d$lod2)]),
   "Delta equals the gap between the best and second-best LOD")

key <- paste(tab$offspring, tab$candidate)
true_key <- paste(ped$offspring, ped$sire)
ok(mean(tab$lod[key %in% true_key], na.rm = TRUE) >
   mean(tab$lod[!(key %in% true_key)], na.rm = TRUE) + 10,
   "true sires score far above non-sires")
ok(mean(tab$mismatches[key %in% true_key]) < 0.5, "true sires rarely mismatch")

tabk <- parentage_lod(gg, ped$offspring, sim$sires, pp, error = 0.01,
                      known = ped$dam, type = "one_known")
dk <- delta_stat(tabk)
ok(mean(dk$candidate == ped$sire[match(dk$offspring, ped$offspring)]) >= 0.95,
   "known-mother configuration recovers the true sire")

tp <- parentage_lod(gg, ped$offspring[1:8], sim$sires, pp, error = 0.01,
                    type = "pair", mothers = sim$dams)
dp <- delta_stat(tp)
trp <- ped[match(dp$offspring, ped$offspring), ]
ok(mean(dp$candidate == trp$sire & dp$candidate2 == trp$dam) >= 0.85,
   "joint parent-pair configuration recovers both parents")

ok(nrow(parentage_lod(gg, ped$offspring[1], sim$sires, pp, error = 0.01,
                      max_mismatch = 0)) <=
   nrow(parentage_lod(gg, ped$offspring[1], sim$sires, pp, error = 0.01)),
   "max_mismatch prunes candidates")
ok(inherits(try(parentage_lod(gg, ped$offspring[1], "nope", pp, error = 0.01),
                silent = TRUE), "try-error"),
   "unknown identifiers are rejected")
ok(inherits(try(parentage_lod(gg, ped$offspring[1], sim$sires, pp, error = 2),
                silent = TRUE), "try-error"),
   "an out-of-range error rate is rejected")

## -------------------------------------------------- null distributions
nb <- lod_null(gg, offspring = ped$offspring[1], freqs = pp, error = 0.01,
               method = "backward", nsim = 1000)
nf <- lod_null(gg, offspring = ped$offspring[1], candidate = ped$sire[1],
               freqs = pp, error = 0.01, method = "forward", nsim = 1000)
ok(length(nb) == 1000L && length(nf) == 1000L, "lod_null() returns nsim values")
ok(all(is.finite(nb)) && all(is.finite(nf)), "simulated LOD scores are finite")
cv <- lod_critical(nb)
ok(nrow(cv) == 2L && cv$critical_lod[1] > cv$critical_lod[2],
   "a stricter alpha gives a lower critical LOD")
ok(near(cv$confidence, 1 - cv$alpha), "confidence is one minus alpha")

# Coverage: a true parent should fall below the alpha quantile about alpha of
# the time. Bounds are loose enough to be stable under the fixed seed.
set.seed(808)
below <- vapply(seq_len(30), function(i) {
  nbi <- lod_null(gg, offspring = ped$offspring[i], freqs = pp, error = 0.01,
                  method = "backward", nsim = 800)
  lod_i <- tab$lod[match(paste(ped$offspring[i], ped$sire[i]), key)]
  lod_i < lod_critical(nbi, alpha = 0.05)$critical_lod
}, logical(1))
ok(mean(below) <= 0.25, "backward null is not grossly anti-conservative")

## -------------------------------------------------------- Delta criterion
ds <- sim_delta(pp, error = 0.01, n_candidates = 20, prop_sampled = 0.6,
                prop_typed = 1, nsim = 1500)
ok(nrow(ds) == 1500L, "sim_delta() returns nsim assignments")
ok(all(ds$correct[ds$correct] %in% TRUE) && !any(ds$correct & !ds$sampled),
   "an unsampled true parent is never counted correct")
dc <- delta_critical(ds, levels = c(0.80, 0.95))
ok(all(is.na(dc$critical_delta) | diff(dc$critical_delta) >= 0),
   "a stricter confidence level needs at least as large a Delta")
ok(all(is.na(dc$prop_correct) | dc$prop_correct >= dc$confidence),
   "the attained proportion correct meets the requested level")

## ------------------------------------------------------- error estimation
sime <- simulate_population(n = 200, n_loci = 10, n_alleles = 8,
                            n_offspring = 400, error = 0.03, seed = 616)
pe <- allele_freqs(sime$genotypes)
prs <- rbind(
  data.frame(parent = sime$pedigree$sire, offspring = sime$pedigree$offspring),
  data.frame(parent = sime$pedigree$dam,  offspring = sime$pedigree$offspring)
)
er <- estimate_error_rates(sime$genotypes, prs, pe)
ok(nrow(er) == 10L && all(er$error > 0), "estimate_error_rates() returns one rate per locus")
ok(abs(mean(er$error) - 0.03) < 0.02, "estimated error rate is close to the truth")
ok(all(er$detectable > 0 & er$detectable < 1), "detection probability is a probability")
ern <- estimate_error_rates(sime$genotypes, prs, pe, method = "naive")
ok(all(is.na(ern$detectable)), "the naive method reports no detection probability")
ok(all(estimate_error_rates(sime$genotypes, prs, pe, min_error = 0.02)$error >= 0.02),
   "min_error is enforced")

## -------------------------------------------------------------- assignment
as1 <- assign_parentage(gg, ped$offspring[1:12], sim$sires, pp, error = 0.01,
                        criterion = "pairwise", nsim = 600)
ok(nrow(as1) == 12L, "assign_parentage() returns one row per offspring")
tr1 <- ped$sire[match(as1$offspring, ped$offspring)]
ok(all(as1$candidate[as1$assigned] == tr1[as1$assigned]),
   "every pairwise assignment made is correct")
ok(!is.null(attr(as1, "lod_table")), "the full LOD table is retained")

as2 <- assign_parentage(gg, ped$offspring[1:12], sim$sires, pp, error = 0.01,
                        criterion = "delta", nsim = 1000, prop_sampled = 0.8)
tr2 <- ped$sire[match(as2$offspring, ped$offspring)]
ok(all(as2$candidate[as2$assigned] == tr2[as2$assigned]),
   "every Delta assignment made is correct")
ok(is.data.frame(summary(as2)), "summary() returns the criterion table")

## ------------------------------------------------------- fast engine
# parentage_lod() reduces the per-candidate likelihood to a table lookup. These
# tests hold it against a direct locus-by-locus evaluation, and hold the two
# accumulation strategies against each other.
ref_total <- function(g, p, e, oid, cid, kid = NULL, type = "both_unknown") {
  oi <- match(oid, g$ids); ci <- match(cid, g$ids)
  ki <- if (is.null(kid)) NA_integer_ else match(kid, g$ids)
  v <- vapply(seq_len(n_loci(g)), function(l)
    lod_locus(g$a1[oi, l], g$a2[oi, l], g$a1[ci, l], g$a2[ci, l], p[[l]], e[l],
              k1 = if (is.null(kid)) NULL else g$a1[ki, l],
              k2 = if (is.null(kid)) NULL else g$a2[ki, l], type = type),
    numeric(1))
  sum(v[is.finite(v)])
}
gf <- simulate_population(n = 60, n_loci = 12, n_alleles = 7, n_offspring = 30,
                          error = 0.01, prop_missing = 0.08, seed = 1234)
gg2 <- gf$genotypes; pp2 <- allele_freqs(gg2); ee <- rep(0.01, 12)
pedf <- gf$pedigree

t1 <- parentage_lod(gg2, pedf$offspring, gf$sires, pp2, 0.01)
ok(max(abs(vapply(seq_len(40), function(z)
  t1$lod[z] - ref_total(gg2, pp2, ee, t1$offspring[z], t1$candidate[z]),
  numeric(1)))) < 1e-9, "fast engine matches locus-by-locus reference (both unknown)")

t2 <- parentage_lod(gg2, pedf$offspring, gf$sires, pp2, 0.01,
                    known = pedf$dam, type = "one_known")
ok(max(abs(vapply(seq_len(40), function(z)
  t2$lod[z] - ref_total(gg2, pp2, ee, t2$offspring[z], t2$candidate[z],
                        t2$known[z], "one_known"),
  numeric(1)))) < 1e-9, "fast engine matches reference (one parent known)")

t3 <- parentage_lod(gg2, pedf$offspring[1:4], gf$sires, pp2, 0.01,
                    type = "pair", mothers = gf$dams)
ok(max(abs(vapply(seq_len(40), function(z)
  t3$lod[z] - ref_total(gg2, pp2, ee, t3$offspring[z], t3$candidate[z],
                        t3$candidate2[z], "pair"),
  numeric(1)))) < 1e-9, "fast engine matches reference (parent pair)")

gi <- parentageLR:::.geno_index(gg2)
oi <- match(pedf$offspring, gg2$ids); ci <- match(gf$sires, gg2$ids)
V <- parentageLR:::.build_V(gi, gi$code[oi, , drop = FALSE], NULL, pp2, ee, "both_unknown")
cc <- t(gi$code[ci, , drop = FALSE]); cc[is.na(cc)] <- gi$G
A <- parentageLR:::.accumulate(V$lod, cc, gi$G, "blas",   length(oi), 12L, length(ci))
B <- parentageLR:::.accumulate(V$lod, cc, gi$G, "lookup", length(oi), 12L, length(ci))
ok(max(abs(A - B)) < 1e-9, "the BLAS and lookup accumulation paths agree")

mx <- parentage_lod(gg2, pedf$offspring, gf$sires, pp2, 0.01, output = "matrix")
ok(inherits(mx, "parentage_lod_matrix") && all(dim(mx$lod) == c(30L, 30L)),
   "matrix output has the right shape")
z <- t1[t1$offspring == pedf$offspring[1] & t1$candidate == gf$sires[2], ]
ok(abs(z$lod - mx$lod[pedf$offspring[1], gf$sires[2]]) < 1e-9,
   "matrix and long output agree")

# The null samplers draw each locus contribution from its exact discrete
# distribution rather than simulating alleles. Check the first two moments
# against a direct allele-level simulation.
set.seed(4242)
og <- parentageLR:::.row_of(gg2, pedf$offspring[1])
nn <- parentageLR:::.sim_backward(og, NULL, pp2, ee, 30000, 1, "both_unknown")
alt <- replicate(30000, {
  tot <- 0
  for (l in seq_len(12)) {
    o1 <- og$a1[l]; o2 <- og$a2[l]
    if (is.na(o1)) next
    pl <- pp2[[l]]
    a <- if (stats::runif(1) < 0.5) o1 else o2
    b <- sample.int(length(pl), 1L, prob = pl)
    if (stats::runif(1) < ee[l]) {
      a <- sample.int(length(pl), 1L, prob = pl); b <- sample.int(length(pl), 1L, prob = pl)
    }
    v <- lod_locus(o1, o2, a, b, pl, ee[l], type = "both_unknown")
    if (is.finite(v)) tot <- tot + v
  }
  tot
})
ok(abs(mean(nn) - mean(alt)) < 0.15 * stats::sd(alt),
   "null sampler reproduces the allele-level mean")
ok(abs(stats::sd(nn) - stats::sd(alt)) < 0.15 * stats::sd(alt),
   "null sampler reproduces the allele-level spread")


cat("\nAll tests passed.\n")
