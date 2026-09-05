suppressMessages(invisible(lapply(list.files("../R", full.names=TRUE), source)))

run_cells <- function(n_loci, n_alleles, err, label, seed = 100,
                      n_off = 400, n_pool = 6000,
                      cands = c(50, 200, 1000, 5000),
                      props = c(1, 0.8, 0.5), nsim = 4000) {
  set.seed(seed)
  s <- simulate_population(n = 2*n_pool, n_loci = n_loci, n_alleles = n_alleles,
                           n_offspring = n_off, error = err, prop_missing = 0.05, seed = seed)
  g <- s$genotypes; p <- allele_freqs(g); ped <- s$pedigree
  sires <- s$sires
  M <- parentage_lod(g, ped$offspring, sires, p, err, output = "matrix")
  LOD <- M$lod
  true_col <- match(ped$sire, sires)

  # trio-specific critical values, one backward simulation per offspring
  crit <- vapply(seq_len(n_off), function(i)
    lod_critical(lod_null(g, offspring = ped$offspring[i], freqs = p, error = err,
                          method = "backward", nsim = nsim), alpha = 0.01)$critical_lod,
    numeric(1))

  out <- list(); k <- 0
  for (nc in cands) for (pi_ in props) {
    set.seed(seed + nc + round(100*pi_))
    inpool <- runif(n_off) < pi_
    asg <- cor_ <- rep(FALSE, n_off)
    dasg <- dcor <- rep(FALSE, n_off)
    # Delta criterion needs a matching population simulation
    ds <- sim_delta(p, error = err, n_candidates = nc, prop_sampled = pi_,
                    prop_typed = 0.95, nsim = if (nc >= 1000) 1200 else 3000)
    dcrit <- delta_critical(ds, levels = 0.95)$critical_delta
    for (i in seq_len(n_off)) {
      others <- sample(setdiff(seq_along(sires), true_col[i]), nc - as.integer(inpool[i]))
      cols <- if (inpool[i]) c(true_col[i], others) else others
      v <- LOD[i, cols]
      o <- order(-v)
      top <- cols[o[1]]; lod1 <- v[o[1]]
      lod2 <- if (length(v) >= 2) v[o[2]] else NA
      npos <- sum(v > 0)
      dlt <- if (npos >= 2) lod1 - lod2 else if (npos == 1) lod1 else NA
      asg[i]  <- lod1 >= crit[i]
      cor_[i] <- asg[i] && inpool[i] && top == true_col[i]
      dasg[i] <- !is.na(dlt) && !is.na(dcrit) && dlt >= dcrit
      dcor[i] <- dasg[i] && inpool[i] && top == true_col[i]
    }
    k <- k + 1
    out[[k]] <- data.frame(markers = label, n_cand = nc, prop_sampled = pi_,
      pw_assigned = mean(asg), pw_precision = if(sum(asg)) mean(cor_[asg]) else NA,
      pw_recall = if(sum(inpool)) mean(cor_[inpool]) else NA,
      d_assigned = mean(dasg), d_precision = if(sum(dasg)) mean(dcor[dasg]) else NA,
      d_recall = if(sum(inpool)) mean(dcor[inpool]) else NA)
  }
  do.call(rbind, out)
}

r1 <- run_cells(15, 10, 0.01, "15 msat")
write.csv(r1, "./study1_msat.csv", row.names = FALSE)
fmt <- function(d) {
  cat(sprintf("%-8s %6s %6s | %9s %9s %8s | %9s %9s %8s\n",
     "markers","nCand","piSamp","PW asgn","PW prec","PW rec","D asgn","D prec","D rec"))
  for (i in seq_len(nrow(d))) with(d[i,], cat(sprintf(
     "%-8s %6d %6.2f | %8.1f%% %8.1f%% %7.1f%% | %8.1f%% %8.1f%% %7.1f%%\n",
     markers, n_cand, prop_sampled, 100*pw_assigned, 100*pw_precision, 100*pw_recall,
     100*d_assigned, 100*d_precision, 100*d_recall)))
}
fmt(r1)
