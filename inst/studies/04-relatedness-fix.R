suppressMessages(invisible(lapply(list.files("../R", full.names=TRUE), source)))
source("./study3_fns.R")
lod_set_ibd <- function(S, err, ibd) {
  L <- matrix(0, nrow(S$O1), S$ncand)
  for (l in seq_len(S$nloc)) { p <- S$freqs[[l]]
    for (j in seq_len(S$ncand)) {
      v <- lod_locus(S$O1[,l],S$O2[,l],S$C1[,j,l],S$C2[,j,l],p=p,error=err,
                     type="both_unknown", ibd=ibd)
      v[!is.finite(v)] <- 0; L[,j] <- L[,j] + v } }
  L
}
# mixture of products:  LOD_mix = -log( sum_c w_c exp(-LOD_c) )
combine <- function(Ls, w) {
  w <- w/sum(w)
  A <- array(unlist(lapply(seq_along(Ls), function(i) log(w[i]) - Ls[[i]])),
             c(dim(Ls[[1]]), length(Ls)))
  m <- apply(A, c(1,2), max)
  -(m + log(apply(exp(sweep(A, c(1,2), m)), c(1,2), sum)))
}
err <- 0.01; n_off <- 800
cat(sprintf("%-11s %4s %6s | %-26s %9s %8s %9s\n","panel","nFS","piSire","model","asgn>=.95","stated","observed"))
for (pn in list(c(15,10,"15 msat"), c(96,2,"96 SNP"))) {
 nloc<-as.integer(pn[1]); nall<-as.integer(pn[2]); lab<-pn[3]
 for (nfs in c(3,10)) for (pi_ in c(1,0.5)) {
  S <- make_set(n_off,nloc,nall,nfs,0,199-nfs,err, 2000+nloc*7+nfs+round(10*pi_))
  set.seed(31+nloc+nfs+round(100*pi_)); present <- runif(n_off) < pi_
  frac <- nfs/(S$ncand-1)
  Lu <- lod_set_ibd(S, err, c(1,0,0))
  La <- lod_set_ibd(S, err, c(0.5,0.5,0))
  models <- list("naive (unrelated only)"       = Lu,
                 "product-of-mixtures"          = lod_set_ibd(S, err, ibd_mixture(unrelated=1-frac, avuncular=frac)),
                 "mixture-of-products"          = combine(list(Lu,La), c(1-frac, frac)),
                 "mixture-of-products (2x rel)" = combine(list(Lu,La), c(1-2*frac, 2*frac)))
  for (nm in names(models)) {
    L <- models[[nm]]; L[!present,1] <- NA
    rownames(L)<-paste0("O",seq_len(n_off)); colnames(L)<-paste0("C",seq_len(S$ncand))
    mm <- structure(list(lod=L, type="both_unknown"), class="parentage_lod_matrix")
    po <- parentage_posterior(mm, prop_sampled=pi_, n_candidates=S$ncand)
    topc <- as.integer(sub("C","",po$candidate)); correct <- present & topc==1L
    hi <- po$posterior >= 0.95
    cat(sprintf("%-11s %4d %6.1f | %-26s %8.1f%% %8.3f %9.3f\n", lab,nfs,pi_,nm,
      100*mean(hi), if(sum(hi)) mean(po$posterior[hi]) else NA, if(sum(hi)) mean(correct[hi]) else NA))
  }
  cat("\n")
 }
}
