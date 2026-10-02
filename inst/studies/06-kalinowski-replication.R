suppressMessages(invisible(lapply(list.files("../../R", full.names=TRUE), source)))
# Kalinowski, Taper & Marshall (2007), stated simulation parameters:
#   "allele frequencies in the population were [0.25, 0.25, 0.2, 0.15, 0.05,
#    0.05, 0.02, 0.01, 0.01, 0.005, 0.005]"
#   "100 unrelated adult males and 100 unrelated adult females"
#   "a genotyping error rate of 0.01"
#   "One hundred thousand simulations were run for each set of parameters."
# Reported: "When six loci were genotyped, and both the actual and assumed
#   genotyping error rates were 0.01 ... The corresponding result for the
#   corrected equations was 73%" at a Delta 0.99 confidence level.
af <- c(0.25, 0.25, 0.2, 0.15, 0.05, 0.05, 0.02, 0.01, 0.01, 0.005, 0.005)
stopifnot(abs(sum(af) - 1) < 1e-12)
cat("allele frequencies sum to", sum(af), "with", length(af), "alleles\n\n")

run <- function(nloc, nsim = 40000, seed = 2007) {
  set.seed(seed)
  p <- replicate(nloc, {v <- af; names(v) <- paste0("A", seq_along(af)); v},
                 simplify = FALSE)
  class(p) <- c("allele_freqs", "list")
  ds <- sim_delta(p, error = 0.01, n_candidates = 100, prop_sampled = 1,
                  prop_typed = 1, nsim = nsim, type = "one_known")
  cv <- delta_critical(ds, levels = 0.99)$critical_delta
  keep <- !is.na(ds$delta) & ds$delta >= cv
  data.frame(loci = nloc, crit99 = cv,
             assigned = mean(keep),
             correct_of_all = mean(keep & ds$correct),
             precision = mean(ds$correct[keep]))
}
cat("Delta 0.99 criterion, 100 candidate males, mother known, error 0.01\n")
cat(sprintf("%5s %9s %10s %16s %11s\n","loci","crit99","assigned","correct/all off","precision"))
res <- do.call(rbind, lapply(c(4, 5, 6, 7, 8), run))
for (i in seq_len(nrow(res))) with(res[i,], cat(sprintf(
  "%5d %9.3f %9.1f%% %15.1f%% %10.1f%%\n", loci, crit99, 100*assigned,
  100*correct_of_all, 100*precision)))
cat("\nKalinowski et al. report 73% for six loci with the corrected equations.\n")
