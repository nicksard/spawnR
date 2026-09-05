suppressMessages(invisible(lapply(list.files("../R", full.names=TRUE), source)))
# Pools where the true sire AND dam each have full sibs present, and both true
# parents are jointly withheld from a fraction of offspring.
build <- function(n_off, nloc, nall, nfs, nun, err, seed) {
  set.seed(seed)
  fr <- lapply(seq_len(nloc), function(l){p<-(0.75^(seq_len(nall)-1))*runif(nall,.5,1.5)
    p<-p/sum(p); names(p)<-paste0("A",seq_len(nall)); p})
  dr <- function(n,l) sample.int(nall,n,TRUE,prob=fr[[l]])
  kid <- function(a1,a2,b1,b2) list(ifelse(runif(length(a1))<.5,a1,a2),
                                    ifelse(runif(length(b1))<.5,b1,b2))
  nM <- 1L+nfs+nun; nF <- nM
  O1<-O2<-matrix(NA_integer_,n_off,nloc)
  M1<-M2<-array(NA_integer_,c(n_off,nM,nloc)); F1<-F2<-array(NA_integer_,c(n_off,nF,nloc))
  for (l in seq_len(nloc)) {
    gs1<-dr(n_off,l);gs2<-dr(n_off,l);gd1<-dr(n_off,l);gd2<-dr(n_off,l)
    hs1<-dr(n_off,l);hs2<-dr(n_off,l);hd1<-dr(n_off,l);hd2<-dr(n_off,l)
    S<-kid(gs1,gs2,gd1,gd2); D<-kid(hs1,hs2,hd1,hd2)
    o<-kid(S[[1]],S[[2]],D[[1]],D[[2]]); O1[,l]<-o[[1]]; O2[,l]<-o[[2]]
    M1[,1,l]<-S[[1]]; M2[,1,l]<-S[[2]]; F1[,1,l]<-D[[1]]; F2[,1,l]<-D[[2]]
    for (j in seq_len(nfs)) {
      x<-kid(gs1,gs2,gd1,gd2); M1[,1+j,l]<-x[[1]]; M2[,1+j,l]<-x[[2]]
      y<-kid(hs1,hs2,hd1,hd2); F1[,1+j,l]<-y[[1]]; F2[,1+j,l]<-y[[2]] }
    for (j in seq_len(nun)) {
      M1[,1+nfs+j,l]<-dr(n_off,l); M2[,1+nfs+j,l]<-dr(n_off,l)
      F1[,1+nfs+j,l]<-dr(n_off,l); F2[,1+nfs+j,l]<-dr(n_off,l) }
    ef<-function(A1,A2){h<-runif(length(A1))<err
      A1[h]<-dr(sum(h),l);A2[h]<-dr(sum(h),l);list(A1,A2)}
    z<-ef(O1[,l],O2[,l]);O1[,l]<-z[[1]];O2[,l]<-z[[2]]
    z<-ef(as.vector(M1[,,l]),as.vector(M2[,,l]));M1[,,l]<-matrix(z[[1]],n_off,nM);M2[,,l]<-matrix(z[[2]],n_off,nM)
    z<-ef(as.vector(F1[,,l]),as.vector(F2[,,l]));F1[,,l]<-matrix(z[[1]],n_off,nF);F2[,,l]<-matrix(z[[2]],n_off,nF)
  }
  list(O1=O1,O2=O2,M1=M1,M2=M2,F1=F1,F2=F2,freqs=fr,nM=nM,nF=nF,nloc=nloc)
}
pair_lod <- function(S, err, ibd) {
  L <- array(0, c(nrow(S$O1), S$nM, S$nF))
  for (l in seq_len(S$nloc)) { p<-S$freqs[[l]]
    for (m in seq_len(S$nM)) for (f in seq_len(S$nF)) {
      v <- lod_locus(S$O1[,l],S$O2[,l], S$F1[,f,l],S$F2[,f,l], p=p, error=err,
                     k1=S$M1[,m,l], k2=S$M2[,m,l], type="pair", ibd=ibd)
      v[!is.finite(v)]<-0; L[,m,f]<-L[,m,f]+v } }
  L
}
combine <- function(Ls, w) { w<-w/sum(w)
  A <- array(unlist(lapply(seq_along(Ls), function(i) log(w[i])-Ls[[i]])), c(dim(Ls[[1]]), length(Ls)))
  mx <- apply(A, 1:3, max); -(mx + log(apply(exp(sweep(A,1:3,mx)),1:3,sum))) }

err<-0.01; n_off<-300; nloc<-15; nall<-10; nfs<-3; nun<-16
S <- build(n_off,nloc,nall,nfs,nun,err,555)
frac <- nfs/(S$nM-1)
Lu <- pair_lod(S, err, c(1,0,0))
La <- pair_lod(S, err, c(0.5,0.5,0))    # avuncular on both alleged parents
Lmix <- combine(list(Lu, La), c(1-frac, frac))
cat(sprintf("%6s | %-24s %9s %8s %9s\n","piPar","model","asgn>=.95","stated","observed"))
for (pi_ in c(1, 0.5)) {
  set.seed(9); present <- runif(n_off) < pi_
  for (nm in c("naive (unrelated)","relatedness-aware")) {
    L <- if (nm=="naive (unrelated)") Lu else Lmix
    M <- matrix(NA_real_, n_off, S$nM*S$nF)
    for (i in seq_len(n_off)) { v <- as.vector(L[i,,]); 
      if (!present[i]) { v[1] <- NA }   # true pair (m=1,f=1) withheld
      M[i,] <- v }
    rownames(M)<-paste0("O",seq_len(n_off)); colnames(M)<-paste0("P",seq_len(ncol(M)))
    mm <- structure(list(lod=M, type="pair"), class="parentage_lod_matrix")
    po <- parentage_posterior(mm, prop_sampled=pi_, n_candidates=ncol(M))
    top <- as.integer(sub("P","",po$candidate)); correct <- present & top==1L
    hi <- po$posterior>=0.95
    cat(sprintf("%6.1f | %-24s %8.1f%% %8.3f %9.3f\n", pi_, nm, 100*mean(hi),
        if(sum(hi)) mean(po$posterior[hi]) else NA, if(sum(hi)) mean(correct[hi]) else NA))
  }
}
