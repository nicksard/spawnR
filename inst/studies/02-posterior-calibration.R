suppressMessages(invisible(lapply(list.files("../R", full.names=TRUE), source)))
set.seed(202)
n_off <- 500; n_pool <- 6000; nloc <- 15; n_all <- 10; err <- 0.01
s <- simulate_population(n = 2*n_pool, n_loci = nloc, n_alleles = n_all,
                         n_offspring = n_off, error = err, prop_missing = 0.05, seed = 202)
g <- s$genotypes; p <- allele_freqs(g); ped <- s$pedigree; sires <- s$sires
LOD <- parentage_lod(g, ped$offspring, sires, p, err, output = "matrix")$lod
true_col <- match(ped$sire, sires)

cands <- c(50, 200, 1000, 5000); props <- c(1, 0.8, 0.5)
rows <- list(); calib <- list(); k <- 0; kk <- 0
for (nc in cands) for (pi_ in props) {
  set.seed(700 + nc + round(100*pi_))
  inpool <- runif(n_off) < pi_
  Lsub <- matrix(NA_real_, n_off, nc)
  tcol <- rep(NA_integer_, n_off)
  for (i in seq_len(n_off)) {
    others <- sample(setdiff(seq_along(sires), true_col[i]), nc - as.integer(inpool[i]))
    cols <- if (inpool[i]) c(true_col[i], others) else others
    Lsub[i, ] <- LOD[i, cols]
    if (inpool[i]) tcol[i] <- 1L
  }
  rownames(Lsub) <- ped$offspring; colnames(Lsub) <- paste0("C", seq_len(nc))
  mm <- structure(list(lod = Lsub, type = "both_unknown"), class = "parentage_lod_matrix")

  for (assumed in c("correct", "0.9", "0.5")) {
    pa <- switch(assumed, correct = pi_, `0.9` = 0.9, `0.5` = 0.5)
    po <- parentage_posterior(mm, prop_sampled = pa)
    topc <- as.integer(sub("C", "", po$candidate))
    correct <- inpool & topc == 1L
    hi <- po$posterior >= 0.95
    fd <- parentage_fdr(po, target = 0.05)
    sel <- if (is.na(fd$cutoff)) rep(FALSE, n_off) else po$posterior >= fd$cutoff
    k <- k + 1
    rows[[k]] <- data.frame(n_cand = nc, prop_sampled = pi_, assumed = assumed,
      p95_assigned = mean(hi), p95_precision = if (sum(hi)) mean(correct[hi]) else NA,
      p95_recall = if (sum(inpool)) mean(correct & hi) / mean(inpool) else NA,
      fdr_target = 0.05, fdr_assigned = mean(sel),
      fdr_realised = if (sum(sel)) 1 - mean(correct[sel]) else NA)
    if (assumed == "correct") {
      kk <- kk + 1
      calib[[kk]] <- data.frame(n_cand = nc, prop_sampled = pi_,
                                post = po$posterior, correct = correct)
    }
  }
}
res <- do.call(rbind, rows); cal <- do.call(rbind, calib)
saveRDS(cal, "./study2_calib.rds")
write.csv(res, "./study2_cells.csv", row.names=FALSE)

cat("=== POSTERIOR, prop_sampled correctly specified ===\n")
cat(sprintf("%6s %7s | %9s %10s %9s | %10s %11s\n","nCand","piSamp","asgn>=.95","prec>=.95","rec>=.95","FDR asgn","FDR realis"))
r <- res[res$assumed=="correct",]
for (i in seq_len(nrow(r))) with(r[i,], cat(sprintf("%6d %7.2f | %8.1f%% %9.1f%% %8.1f%% | %9.1f%% %10.1f%%\n",
  n_cand, prop_sampled, 100*p95_assigned, 100*p95_precision, 100*p95_recall, 100*fdr_assigned, 100*fdr_realised)))

cat("\n=== CALIBRATION (all cells pooled, correct prior) ===\n")
br <- c(0,.5,.8,.9,.95,.99,.999,1)
cal$bin <- cut(cal$post, br, include.lowest=TRUE)
agg <- aggregate(cbind(correct, post) ~ bin, cal, function(z) c(mean(z), length(z)))
tb <- do.call(rbind, lapply(levels(cal$bin), function(b){
  d <- cal[cal$bin==b,]; if(!nrow(d)) return(NULL)
  data.frame(bin=b, n=nrow(d), mean_post=mean(d$post), obs_correct=mean(d$correct))}))
print(tb, row.names=FALSE, digits=3)

cat("\n=== SENSITIVITY: prior misspecified (precision at posterior >= 0.95) ===\n")
cat(sprintf("%6s %7s | %10s %10s %10s\n","nCand","true pi","correct","assume .9","assume .5"))
for (nc in cands) for (pi_ in props) {
  a <- res[res$n_cand==nc & res$prop_sampled==pi_,]
  cat(sprintf("%6d %7.2f | %9.1f%% %9.1f%% %9.1f%%\n", nc, pi_,
    100*a$p95_precision[a$assumed=="correct"], 100*a$p95_precision[a$assumed=="0.9"],
    100*a$p95_precision[a$assumed=="0.5"]))
}
