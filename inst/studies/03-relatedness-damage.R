suppressMessages(invisible(lapply(list.files("../R", full.names=TRUE), source)))
source("./study3_fns.R")

err <- 0.01; n_off <- 800
panels <- list(c(15,10,"15 msat x10"), c(8,5,"8 msat x5"), c(5,4,"5 msat x4"),
               c(96,2,"96 SNP"), c(300,2,"300 SNP"))
cat(sprintf("%-12s %4s %4s %6s | %9s %10s %9s | %s\n",
    "panel","nFS","nHS","piSire","asgn>=.95","precision","FP rate","calib(stated/obs)"))
res <- list(); z <- 0
for (pn in panels) {
  nloc <- as.integer(pn[1]); nall <- as.integer(pn[2]); lab <- pn[3]
  for (k in list(c(0,0), c(3,0), c(10,0), c(3,10))) {
    for (pi_ in c(1, 0.5)) {
      S <- make_set(n_off, nloc, nall, k[1], k[2], 199 - k[1] - k[2], err,
                    2000 + nloc*7 + k[1] + round(10*pi_))
      L <- lod_set(S, err)
      set.seed(31 + nloc + k[1] + round(100*pi_))
      present <- runif(n_off) < pi_
      keep <- if (all(present)) seq_len(S$ncand) else NULL
      Lp <- L
      # drop the true sire (column 1) for offspring whose sire was not sampled
      Lp[!present, 1] <- NA
      rownames(Lp) <- paste0("O", seq_len(n_off)); colnames(Lp) <- paste0("C", seq_len(S$ncand))
      mm <- structure(list(lod = Lp, type = "both_unknown"), class = "parentage_lod_matrix")
      po <- parentage_posterior(mm, prop_sampled = pi_, n_candidates = S$ncand)
      topc <- as.integer(sub("C","",po$candidate))
      correct <- present & topc == 1L
      hi <- po$posterior >= 0.95
      fp <- if (sum(hi)) mean(!correct[hi]) else NA
      z <- z + 1
      res[[z]] <- data.frame(panel=lab, n_fs=k[1], n_hs=k[2], pi=pi_,
        assigned=mean(hi), precision=if(sum(hi)) mean(correct[hi]) else NA,
        fp=fp, stated=if(sum(hi)) mean(po$posterior[hi]) else NA)
      cat(sprintf("%-12s %4d %4d %6.1f | %8.1f%% %9.1f%% %8.1f%% | %.3f / %.3f\n",
          lab,k[1],k[2],pi_,100*mean(hi),100*res[[z]]$precision,100*fp,
          res[[z]]$stated, res[[z]]$precision))
    }
  }
}
saveRDS(do.call(rbind,res), "./study3b.rds")
