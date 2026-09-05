make_set <- function(n_off, nloc, nall, n_fs, n_hs, n_unrel, err, seed) {
  set.seed(seed)
  fr <- lapply(seq_len(nloc), function(l) {
    p <- (0.75^(seq_len(nall)-1)) * runif(nall, .5, 1.5); p <- p/sum(p)
    names(p) <- paste0("A", seq_len(nall)); p })
  draw <- function(n, l) sample.int(nall, n, TRUE, prob = fr[[l]])
  kid  <- function(a1,a2,b1,b2,l) list(ifelse(runif(length(a1))<.5,a1,a2),
                                       ifelse(runif(length(b1))<.5,b1,b2))
  ncand <- 1L + n_fs + n_hs + n_unrel
  O1<-O2<-matrix(NA_integer_,n_off,nloc); C1<-C2<-array(NA_integer_,c(n_off,ncand,nloc))
  for (l in seq_len(nloc)) {
    gs1<-draw(n_off,l); gs2<-draw(n_off,l); gd1<-draw(n_off,l); gd2<-draw(n_off,l)
    s <- kid(gs1,gs2,gd1,gd2,l)                       # true sire
    d1<-draw(n_off,l); d2<-draw(n_off,l)              # dam of focal offspring
    o <- kid(s[[1]],s[[2]],d1,d2,l)
    O1[,l]<-o[[1]]; O2[,l]<-o[[2]]
    C1[,1,l]<-s[[1]]; C2[,1,l]<-s[[2]]
    if (n_fs) for (j in seq_len(n_fs)) {
      f <- kid(gs1,gs2,gd1,gd2,l); C1[,1+j,l]<-f[[1]]; C2[,1+j,l]<-f[[2]] }
    if (n_hs) for (j in seq_len(n_hs)) {
      od1<-draw(n_off,l); od2<-draw(n_off,l)
      h <- kid(gs1,gs2,od1,od2,l); C1[,1+n_fs+j,l]<-h[[1]]; C2[,1+n_fs+j,l]<-h[[2]] }
    if (n_unrel) for (j in seq_len(n_unrel)) {
      C1[,1+n_fs+n_hs+j,l]<-draw(n_off,l); C2[,1+n_fs+n_hs+j,l]<-draw(n_off,l) }
  }
  # genotyping error: random genotype replacement
  errify <- function(A1,A2,l) { h <- runif(length(A1))<err
    A1[h]<-sample.int(nall,sum(h),TRUE,prob=fr[[l]]); A2[h]<-sample.int(nall,sum(h),TRUE,prob=fr[[l]])
    list(A1,A2) }
  for (l in seq_len(nloc)) {
    z<-errify(O1[,l],O2[,l],l); O1[,l]<-z[[1]]; O2[,l]<-z[[2]]
    z<-errify(as.vector(C1[,,l]),as.vector(C2[,,l]),l)
    C1[,,l]<-matrix(z[[1]],n_off,ncand); C2[,,l]<-matrix(z[[2]],n_off,ncand)
  }
  list(O1=O1,O2=O2,C1=C1,C2=C2,freqs=fr,ncand=ncand,nloc=nloc)
}
lod_set <- function(S, err) {
  L <- matrix(0, nrow(S$O1), S$ncand)
  for (l in seq_len(S$nloc)) {
    p <- S$freqs[[l]]
    for (j in seq_len(S$ncand)) {
      v <- lod_locus(S$O1[,l], S$O2[,l], S$C1[,j,l], S$C2[,j,l], p = p,
                     error = err, type = "both_unknown")
      v[!is.finite(v)] <- 0; L[,j] <- L[,j] + v
    }
  }
  L
}
