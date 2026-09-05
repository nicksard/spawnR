# Clean-room derivation of the likelihood equations in `certusR`

This note records how every equation implemented in `certusR` was obtained.
The package was written **without reading any existing implementation's source
code**. In particular, the source of `github.com/mahmood225/PairwisePaternity`
(non-commercial licence, GPL-incompatible) was never fetched, cloned or read.
The only inputs were the published equations and prose of:

* **Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)**
  Statistical confidence for likelihood-based paternity inference in natural
  populations. *Molecular Ecology* **7**, 639-655. doi:10.1046/j.1365-294x.1998.00374.x
* **Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007)**
  Revising how the computer program CERVUS accommodates genotyping error
  increases success in paternity assignment. *Molecular Ecology* **16**,
  1099-1106. doi:10.1111/j.1365-294X.2007.03089.x
  (see also the corrigendum, *Molecular Ecology* **19**, 1512,
  doi:10.1111/j.1365-294x.2010.04544.x)
* **Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
  (2025)** Pairwise paternity assignment with forward-backward simulations:
  refining CERVUS using trio-based likelihood and locus-specific error rates.
  *Ecology and Evolution* **15**(10), e72230. doi:10.1002/ece3.72230


## 1. Why the published coefficients were re-derived rather than transcribed

The error-expansion coefficients printed in these papers are demonstrably
mis-typeset, and this is acknowledged in the literature: the 2010 corrigendum to
Kalinowski et al. (2007) states that "the likelihood equations printed in the
Appendix of Kalinowski et al. (2007) contained typesetting errors", while noting
that the errors never affected the CERVUS software itself.

The same damage is visible in the equations of Amiri Roudbar et al. (2025).
Their equations 1, 2, 5 and 6 print a third term with coefficient
`eps^2 (1-eps)^3`, and their equations 3 and 4 print a second term with
coefficient `eps (1-eps)^2`. Neither set of coefficients sums to one over the
mutually exclusive error configurations, so neither can be a probability
decomposition as printed. The most parsimonious reading is that the
multinomial multipliers `3` and `2` were absorbed into the exponents during
typesetting: `3 eps^2 (1-eps)` became `eps^2 (1-eps)^3`, and
`2 eps (1-eps)` became `eps (1-eps)^2`.

Rather than propagate a defect, the expansion below is derived from first
principles out of the *error model prose*, which all three papers state
unambiguously and consistently. The derivation reproduces the published term
**structure** exactly - including the distinctive three-part bracket
`T(o|m) + T(o|a) + P(o)` of equation 1 - while yielding coefficients that are a
proper probability decomposition. Where this implementation departs from a
printed coefficient, it is flagged in `?likelihood_notes` and in the
`LOD_expansion_sums_to_one` regression test.


## 2. The error model (random genotype replacement)

Marshall et al. (1998), restated by Kalinowski et al. (2007) and adopted with a
locus-specific rate by Amiri Roudbar et al. (2025):

> when a genotyping error occurs, the true genotype is replaced by a genotype
> selected at random under Hardy-Weinberg assumptions.

Writing `G` for a true genotype, `g` for the observed genotype, `e_l` for the
genotyping error rate at locus `l`, and `P(g)` for the Hardy-Weinberg frequency
of `g`:

```
Pr(observe g | true G) = (1 - e_l) * 1{G = g}  +  e_l * P(g)                (E1)
```

Errors are independent across individuals and across loci. Note the identity

```
sum_G P(G) Pr(g | G) = (1 - e_l) P(g) + e_l P(g) = P(g)                     (E2)
```

so genotyping error under this model leaves the marginal genotype frequencies
unchanged. This is why the denominator of the "both parents unknown" case is
*not* modified by error (section 4.2 below).


## 3. Mendelian transition probabilities

Let `p_i` be the frequency of allele `A_i` at the locus.

**Both parents' genotypes conditioned on** (`m = (m1,m2)`, `f = (f1,f2)`):

```
T(o | m, f) = (1/4) * sum_{a in m} sum_{b in f} 1{ {a,b} = o }              (T2)
```

**One parent conditioned on, the other drawn from the population** (`f = (f1,f2)`):
the known parent transmits each of its two alleles with probability 1/2, and the
unsampled parent contributes an allele drawn at Hardy-Weinberg frequency. For
`o = (x, y)`:

```
T(o | f) = (1/2) * [ h(f1) + h(f2) ]                                        (T1)

  where, if x = y :  h(a) = p_x  if a = x, else 0
        if x != y:  h(a) = p_y  if a = x
                    h(a) = p_x  if a = y
                    h(a) = 0    otherwise
```

Expanding (T1) gives the familiar table:

| parent | offspring | T(o \| f)      |
|--------|-----------|----------------|
| AiAi   | AiAi      | p_i            |
| AiAi   | AiAj      | p_j            |
| AiAj   | AiAj      | (p_i + p_j)/2  |
| AiAj   | AiAi      | p_i / 2        |
| AiAj   | AiAk      | p_k / 2        |
| any    | no shared allele | 0       |

Hardy-Weinberg genotype frequency:

```
P(o) = p_x^2        if x = y
P(o) = 2 p_x p_y    if x != y                                              (P0)
```

Two identities used repeatedly below, both consequences of Hardy-Weinberg:

```
sum_o T(o | m, f) = 1 ,   sum_o T(o | f) = 1 ,   sum_o P(o) = 1            (I1)
sum_F P(F) T(o | m, F) = T(o | m)                                          (I2)
sum_M sum_F P(M) P(F) T(o | M, F) = P(o)                                   (I3)
```

(I2) says that marginalising the unsampled father out of the trio transition
probability returns the single-parent transition probability; (I3) says that the
offspring of two random parents is itself a random individual. Both are checked
numerically in the test suite.


## 4. Marginalising the error model over unobserved true genotypes

### 4.1 One parent known (mother `m` known, alleged father `a` tested)

The likelihood of the three *observed* genotypes under H1 (the alleged father is
the true father) marginalises over the three *true* genotypes:

```
L(H1) = sum_{M,A,O} P(M) P(A) T(O|M,A) Pr(g_m|M) Pr(g_a|A) Pr(g_o|O)
```

Substituting (E1) expands this into 2^3 = 8 terms, one per subset of individuals
that was mistyped. Applying (I1)-(I3) to collapse each term:

| mistyped | weight        | collapses to                     |
|----------|---------------|----------------------------------|
| none     | (1-e)^3       | T(g_o \| g_m, g_a)               |
| offspring| e(1-e)^2      | P(g_o)          (by I1)          |
| father   | e(1-e)^2      | T(g_o \| g_m)   (by I2)          |
| mother   | e(1-e)^2      | T(g_o \| g_a)   (by I2)          |
| any two  | e^2(1-e), x3  | P(g_o)          (by I1, I3)      |
| all three| e^3           | P(g_o)                           |

giving, after dividing out the common factor `P(g_m) P(g_a)`:

```
L(H1) ~ (1-e)^3 T(o|m,a)
      + e(1-e)^2 [ T(o|m) + T(o|a) + P(o) ]
      + 3 e^2 (1-e) P(o)
      + e^3 P(o)                                                           (L1)
```

This is *structurally identical* to equation 1 of Amiri Roudbar et al. (2025),
including the three-part bracket, and differs only in the third coefficient
(`3 e^2 (1-e)` here, `e^2 (1-e)^3` as printed).

Under H2 (the alleged father is unrelated, so the true father is a random male)
the alleged father's genotype drops out of the transition probability entirely:

```
L(H2) ~ (1-e)^3 T(o|m)
      + e(1-e)^2 [ T(o|m) + 2 P(o) ]
      + 3 e^2 (1-e) P(o)
      + e^3 P(o)                                                           (L2)
```

The `2 P(o)` arises because mistyping the *offspring* and mistyping the *mother*
both collapse to `P(o)`, while mistyping the irrelevant alleged father leaves
`T(o|m)` intact. Equation 2 of Amiri Roudbar et al. (2025) prints only
`T(o|m) + P(o)` here; one `P(o)` appears to have been dropped.

Both (L1) and (L2) sum to one over `o`, since
`(1-e)^3 + 3e(1-e)^2 + 3e^2(1-e) + e^3 = 1`.

### 4.2 Both parents unknown (alleged father `a` tested alone)

Only two individuals can be mistyped, so with `u0 = (1-e)^2`, `u1 = e(1-e)`,
`u2 = e^2`:

```
L(H1) ~ (1-e)^2 T(o|a) + 2 e(1-e) P(o) + e^2 P(o)                          (L3)
L(H2) ~ (1-e)^2 P(o)   + 2 e(1-e) P(o) + e^2 P(o)  =  P(o)                 (L4)
```

(L4) collapses exactly, which is the (E2) invariance again: the denominator of
this case is untouched by genotyping error. Equations 3 and 4 of Amiri Roudbar
et al. (2025) print `e(1-e)^2` for the middle coefficient; `2 e(1-e)` is required
for the weights to sum to one, and matches the independently published
restatement of the Kalinowski model by Huang et al. (2018, *Genetics* 210:1467).

### 4.3 Both parents tested jointly

H1 is (L1) with the alleged mother in place of the known mother. Under H2 both
alleged parents are unrelated to the offspring, so every term collapses to
`P(o)` by (I1) and (I3):

```
L(H1) ~ (1-e)^3 T(o|am,af)
      + e(1-e)^2 [ T(o|am) + T(o|af) + P(o) ]
      + 3 e^2 (1-e) P(o) + e^3 P(o)                                        (L5)
L(H2) ~ P(o)                                                               (L6)
```

### 4.4 LOD

```
LOD_l = ln( L(H1)_l / L(H2)_l ),      LOD = sum_l LOD_l                     (L7)
```

Loci at which any required genotype is missing contribute exactly zero.

### 4.5 The zero-error limit

Setting `e = 0` recovers the classical Marshall et al. (1998) ratios, which is
the primary structural check on the expansion:

```
one parent known : LOD_l = ln[ T(o|m,a) / T(o|m) ]
both unknown     : LOD_l = ln[ T(o|a)   / P(o)   ]
joint pair       : LOD_l = ln[ T(o|am,af) / P(o) ]
```

Setting `e = 1` gives `LOD_l = 0` at every locus, as it must.


## 5. Trio-specific null distributions

Amiri Roudbar et al. (2025) replace CERVUS's single population-level Delta
distribution with a reference distribution conditioned on the actual genotypes
of the trio being tested, on the grounds that the alleles carried by the
particular parent and offspring - and the pattern of missing data - change the
distribution of the statistic.

**Forward simulation** (their section 2.3; used for *paternity testing*: is this
specific male the father?). Holding the observed candidate genotype fixed, and
the known parent's genotype fixed if available: draw an offspring genotype by
Mendelian sampling from the candidate, taking the second allele from the known
parent if present and otherwise from population allele frequencies; apply
genotyping error and missing data to the simulated offspring; compute the LOD of
the candidate. Repeated `nsim` times, this is the distribution of LOD **when the
candidate really is the parent**.

**Backward simulation** (their section 2.4; used for *assignment*: which of the
candidates is the father?). Holding the observed offspring genotype fixed: draw
a random true father by sampling his first allele from the offspring's alleles -
excluding the allele traceable to the known mother, where it is traceable - and
his second allele from population allele frequencies; apply genotyping error and
missing data to the simulated father; compute the LOD. Repeated `nsim` times,
this is the distribution of LOD **over the true fathers compatible with this
offspring**.

Both distributions describe *true* parents, so the decision rule is a **lower**
tail: the critical value at significance `alpha` is the `alpha` quantile, and a
candidate is accepted when its observed LOD is at least that large. The paper
uses `alpha = 0.05` and `alpha = 0.01` and `nsim = 10000`.

Note that this is the opposite tail from the classical CERVUS criterion, which
places a threshold on the *upper* tail of a distribution built from *unrelated*
candidates. Both are provided here: `lod_null()` / `lod_critical()` implement
the trio-specific true-parent rule, and `sim_delta()` / `delta_critical()`
implement the Marshall population-level Delta rule.


## 6. Delta and its critical values

Marshall et al. (1998) define, over the candidate parents not excluded for a
given offspring,

```
Delta = LOD(1) - LOD(2)   if two or more candidates remain
Delta = LOD(1)            if exactly one candidate remains                  (D1)
```

Confidence is defined as the proportion of assignments, among all simulated
assignments with Delta at or above a threshold, in which the top candidate is
the true parent. The critical Delta at confidence level `c` is therefore the
smallest threshold whose implied proportion reaches `c`. Simulation draws a true
parent and offspring under Hardy-Weinberg, adds `n_candidates - 1` unrelated
candidates, includes the true parent in the candidate set with probability
`prop_sampled`, and applies typing error and missing data at the stated rates.
CERVUS's conventional levels, relaxed 80% and strict 95%, are the defaults here.


## 7. Locus-specific error rates

Amiri Roudbar et al. (2025), section 2.6, estimate a per-locus rate by counting
mismatching parent-offspring pairs at each locus. The paper does not state the
denominator used to turn those counts into a rate, so `estimate_error_rates()`
offers two documented estimators:

* `method = "naive"`: `e_l = m_l / n_l`, the raw proportion of scored
  parent-offspring pairs that mismatch at locus `l`. This under-estimates the
  error rate, because many errors produce a genotype that is still Mendelian-
  compatible and therefore invisible.

* `method = "corrected"` (default): inverts the probability that an error is
  actually *detectable*. Under the replacement model an error somewhere in a
  parent-offspring pair occurs with probability `2 e (1 - e)` for exactly one
  member, and is visible only if the replacement genotype is incompatible with
  the other member's. Writing `d_l` for that detection probability,

  ```
  d_l = 1 - sum_{gp} sum_{go} P(gp) P(go) 1{ T(go | gp) > 0 }
  ```

  the expected mismatch rate is `m_l / n_l ~ 2 e_l (1 - e_l) d_l`, which is
  solved for `e_l` by taking the smaller root. `d_l` is computed exactly from
  the estimated allele frequencies.

The second estimator is this package's own derivation, not a published
equation, and is documented as such.
