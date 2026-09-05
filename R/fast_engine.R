# Internal fast engine ------------------------------------------------------
#
# The per-locus LOD contribution depends on the candidate only through the
# candidate's genotype, and a locus with k alleles has only G = k(k+1)/2
# genotypes. So instead of evaluating the likelihood once per candidate, the
# engine evaluates it once per distinct genotype - G times per locus - and then
# reduces the offspring-by-candidate score to a table lookup.
#
# Writing V[i, l, a] for the LOD contribution at locus l when the candidate
# carries genotype a, the total is
#
#     LOD[i, j] = sum_l V[i, l, code(j, l)]
#
# which can be accumulated two ways:
#
#   "blas"   LOD = sum_a V[, , a] %*% C_a, where C_a is the nl-by-nc indicator
#            of candidate genotype a. G dense matrix products, so this wins when
#            G is small - biallelic SNP panels above all.
#   "lookup" one gather of an nl-by-nc slice per locus. Independent of G, so it
#            wins for microsatellites, where G is large but loci are few.
#
# The switch is made on G. Both paths are exercised against the reference
# implementation in the test suite.

.geno_index <- function(g) {
  nl <- n_loci(g)
  code <- matrix(NA_integer_, n_ind(g), nl)
  tabs <- vector("list", nl)
  for (l in seq_len(nl)) {
    k <- length(g$alleles[[l]])
    gi <- unlist(lapply(seq_len(k), function(i) rep(i, k - i + 1L)))
    gj <- unlist(lapply(seq_len(k), function(i) seq.int(i, k)))
    a1 <- pmin(g$a1[, l], g$a2[, l])
    a2 <- pmax(g$a1[, l], g$a2[, l])
    code[, l] <- match(a1 + (a2 - 1L) * k, gi + (gj - 1L) * k)
    tabs[[l]] <- list(gi = gi, gj = gj, G = length(gi))
  }
  list(code = code, tabs = tabs, G = max(vapply(tabs, `[[`, integer(1), "G")) + 1L)
}

# V arrays of dimension (n_off, n_loci, G) for the LOD and the mismatch flag.
# Column G is the untyped-candidate slot and stays zero. The transition
# probabilities are computed once and used for both quantities.
.build_V <- function(gi, ocode, kcode, freqs, e, type, ibd = NULL) {
  ibd <- .norm_ibd(ibd, type)
  no <- nrow(ocode); nl <- ncol(ocode); G <- gi$G
  Vl <- array(0, c(no, nl, G))
  Vm <- array(0, c(no, nl, G))
  for (l in seq_len(nl)) {
    tb <- gi$tabs[[l]]; p <- freqs[[l]]; ee <- e[l]
    o_code <- ocode[, l]
    o1 <- tb$gi[o_code]; o2 <- tb$gj[o_code]
    Po <- geno_freq(o1, o2, p)
    if (type == "both_unknown") {
      u0 <- (1 - ee)^2; u1 <- ee * (1 - ee); u2 <- ee^2
      k1 <- k2 <- NULL; T_k <- NULL
    } else {
      kc <- kcode[, l]; k1 <- tb$gi[kc]; k2 <- tb$gj[kc]
      w0 <- (1 - ee)^3; w1 <- ee * (1 - ee)^2
      w23 <- 3 * ee^2 * (1 - ee) + ee^3
      T_k <- trans_prob(o1, o2, k1, k2, p)
    }
    for (a in seq_len(tb$G)) {
      c1 <- tb$gi[a]; c2 <- tb$gj[a]
      T_c <- trans_prob(o1, o2, c1, c2, p)
      if (type == "both_unknown") {
        num <- u0 * T_c + (2 * u1 + u2) * Po
        den <- if (is.null(ibd)) Po else
          u0 * (ibd$m[1] * Po + ibd$m[2] * T_c + ibd$m[3] * as.numeric(o_code == a)) +
          (2 * u1 + u2) * Po
        bad <- T_c == 0
      } else {
        T_kc <- trans_prob_pair(o1, o2, k1, k2, c1, c2, p)
        num <- w0 * T_kc + w1 * (T_k + T_c + Po) + w23 * Po
        den <- if (type == "pair") {
                 if (is.null(ibd)) Po else {
                   aa <- ibd$m[2]; bb <- ibd$f[2]
                   w0 * (aa * bb * T_kc + aa * (1 - bb) * T_k +
                         (1 - aa) * bb * T_c + (1 - aa) * (1 - bb) * Po) +
                   w1 * (Po + aa * T_k + (1 - aa) * Po + bb * T_c + (1 - bb) * Po) +
                   w23 * Po
                 }
               } else if (is.null(ibd)) w0 * T_k + w1 * (T_k + 2 * Po) + w23 * Po
               else w0 * (ibd$m[2] * T_kc + (1 - ibd$m[2]) * T_k) +
                    w1 * (Po + T_k + ibd$m[2] * T_c + (1 - ibd$m[2]) * Po) + w23 * Po
        bad <- T_kc == 0
      }
      v <- .safe_log_ratio(num, den)
      cmp <- !is.na(v)
      v[!cmp] <- 0
      inf <- !is.finite(v)
      if (any(inf)) v[inf] <- sign(v[inf]) * 700
      bad[is.na(bad)] <- FALSE
      Vl[, l, a] <- v
      Vm[, l, a] <- as.numeric(bad)
    }
  }
  list(lod = Vl, mism = Vm, G = G)
}

.accumulate <- function(V, cc, G, method, no, nl, nc) {
  if (method == "blas") {
    acc <- matrix(0, no, nc)
    for (a in seq_len(G)) {
      Ca <- (cc == a) * 1.0
      if (!any(Ca != 0)) next
      acc <- acc + array(V[, , a], c(no, nl)) %*% Ca
    }
    acc
  } else {
    acc <- matrix(0, no, nc)
    idx_i <- rep(seq_len(no), times = nc)
    for (l in seq_len(nl)) {
      Sl <- array(V[, l, ], c(no, G))
      acc <- acc + matrix(Sl[cbind(idx_i, rep(cc[l, ], each = no))], no, nc)
    }
    acc
  }
}

.pick_method <- function(G, nl, nc) if (G <= 48L) "blas" else "lookup"

# Offspring-by-candidate LOD, mismatch and comparability matrices.
.lod_matrix <- function(gi, off_idx, cand_idx, known_idx, freqs, e, type,
                        chunk = 4e6, ibd = NULL) {
  nl <- ncol(gi$code); G <- gi$G
  nc <- length(cand_idx); no_all <- length(off_idx)
  ccode <- gi$code[cand_idx, , drop = FALSE]
  cc <- t(ccode); cc[is.na(cc)] <- G
  method <- .pick_method(G, nl, nc)

  # Comparability needs no per-genotype table: a locus counts when the
  # offspring, the candidate and (where used) the known parent are all typed.
  ok_c <- matrix(as.numeric(!is.na(t(ccode))), nl, nc)
  ok_o <- !is.na(gi$code[off_idx, , drop = FALSE])
  if (type != "both_unknown") ok_o <- ok_o & !is.na(gi$code[known_idx, , drop = FALSE])
  ncm <- (ok_o * 1.0) %*% ok_c

  per <- max(1L, floor(chunk / (nl * G)))
  out_l <- out_m <- matrix(0, no_all, nc)
  for (s in seq(1L, no_all, by = per)) {
    ii <- s:min(s + per - 1L, no_all)
    ocode <- gi$code[off_idx[ii], , drop = FALSE]
    kcode <- if (type == "both_unknown") NULL else gi$code[known_idx[ii], , drop = FALSE]
    V <- .build_V(gi, ocode, kcode, freqs, e, type, ibd)
    out_l[ii, ] <- .accumulate(V$lod,  cc, G, method, length(ii), nl, nc)
    out_m[ii, ] <- .accumulate(V$mism, cc, G, method, length(ii), nl, nc)
  }
  list(lod = out_l, mism = round(out_m), ncmp = round(ncm), method = method)
}
