# certusR

Likelihood-ratio ("LOD score") parentage inference for codominant markers, with
locus-specific genotyping error and simulation-based confidence. Base R only —
no dependencies.

## The name

Latin *certus* — "determined, resolved, certain" — is the past participle of
*cernere*, to sift apart, to distinguish, to decide; the root that also gives
*discern*, *criterion* and *crisis*. It is one letter from *cervus*, the red
deer of Rum whose pedigree motivated Marshall et al. (1998) and named CERVUS.
That distinction is the point of the package: CERVUS returns a score, and what
is added here is how far that score can be trusted for the particular trio in
front of you.

## What it does

LOD scores in the three configurations of Kalinowski et al. (2007):

| configuration     | `type`           | hypothesis pair |
|-------------------|------------------|-----------------|
| one parent known  | `"one_known"`    | candidate is the second parent, vs. unrelated |
| both unknown      | `"both_unknown"` | candidate is a parent, vs. unrelated |
| both jointly      | `"pair"`         | candidate pair are the parents, vs. both unrelated |

Two ways to attach confidence:

* **Trio-specific** (Amiri Roudbar et al. 2025) — `lod_null()` builds a
  reference distribution by forward or backward simulation conditioned on the
  actual genotypes and missing-data pattern of the case being tested;
  `lod_critical()` takes the lower-tail quantile.
* **Population-level** (Marshall et al. 1998) — `sim_delta()` and
  `delta_critical()` reproduce the Delta criterion that CERVUS reports, at the
  conventional relaxed (80%) and strict (95%) levels.

Per-locus genotyping error rates can be estimated from known parent–offspring
pairs with `estimate_error_rates()`.

## Install

```r
install.packages("certusR_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Usage

```r
library(certusR)

g <- genotypes(my_wide_table, id_col = "id")   # id, then 2 columns per locus
p <- allele_freqs(g)

# per-locus error rates from a known pedigree, if you have one
er <- estimate_error_rates(g, known_pairs, p)

# LOD scores for every offspring against every candidate father
tab <- parentage_lod(g, offspring = offs, candidates = males,
                     freqs = p, error = er$error, type = "both_unknown")
delta_stat(tab)

# assignment with trio-specific criteria
assign_parentage(g, offs, males, p, error = er$error,
                 criterion = "pairwise", nsim = 10000)
```

`simulate_population()` generates a population with a known pedigree for
calibration and testing.

## Performance

Base R throughout, no compiled code. The engine evaluates the likelihood once
per distinct genotype per locus - a locus with *k* alleles has only *k(k+1)/2*
of them - and reduces the offspring-by-candidate score to a table lookup,
accumulated either with dense matrix products (SNP panels) or per-locus gathers
(microsatellites).

Measured on a single core:

| task                                                    | time   |
|---------------------------------------------------------|--------|
| `parentage_lod`, 300 offspring x 300 cand x 300 SNPs     | 0.5 s  |
| `parentage_lod`, 1000 x 2000 x 300 SNPs (matrix output)  | 3.3 s  |
| `parentage_lod`, 500 x 500 x 15 microsatellites          | 0.7 s  |
| `lod_null`, 10000 replicates, 300 SNPs                   | 0.19 s |
| `assign_parentage` pairwise, 500 x 500 x 15 msats        | 4.3 s  |
| `assign_parentage` pairwise, 500 x 500 x 300 SNPs        | 48 s   |

For anything larger, `output = "matrix"` avoids building a long-format table of
millions of rows, which by then costs more than the arithmetic.

An Rcpp version of the accumulation kernel is in `inst/optional-cpp/` with
benchmarks. It is roughly 4x faster than the matrix-product path on SNP panels,
which is well under 2x end to end, and it would cost every user a compiler. The
package stays dependency-free unless profiling on real data says otherwise.

## Choosing a confidence criterion

The two criteria control different things. The trio-specific `"pairwise"` rule
bounds the chance of *rejecting a true parent*; it does not bound the chance
that an unrelated candidate clears the same threshold, and that grows with the
pool. At a nominal 99% level with 15 microsatellites and the true sire always
present, the share of assignments that were wrong was 0% with 100 or 500
candidates, 3.1% with 2000, and 5.6% with 5000. The Marshall `"delta"` rule
defines confidence as the proportion of assignments that are correct and
simulates the whole pool, so it is the one to report when screening broadly.
See `?assign_parentage`.


## Provenance

This is a clean-room implementation. The likelihood formulation was
reconstructed from the published equations and prose of the three papers below.
No existing implementation's source code was consulted; in particular the source
of `PairwisePaternity` (non-commercial licence, GPL-incompatible) was never
fetched or read.

The genotyping-error expansion is **re-derived from the stated
random-genotype-replacement model rather than transcribed**, because the
coefficients printed in the sources do not form a probability decomposition. A
corrigendum to Kalinowski et al. (2007) records typesetting errors in that
paper's appendix, and the coefficients of Amiri Roudbar et al. (2025) show the
same damage — as printed, their equation 1 sums to 0.924 rather than 1 at an
error rate of 0.2. The re-derivation reproduces the published term *structure*
exactly, including the distinctive `T(o|m) + T(o|a) + P(o)` bracket, while
restoring the multinomial multipliers so that every hypothesis is a proper
distribution. That property is asserted in the test suite, and the zero-error
limit is checked against the classical Marshall ratios.

The full term-by-term derivation is in
[`inst/notes/derivation.md`](inst/notes/derivation.md), also reachable at
runtime:

```r
file.show(system.file("notes", "derivation.md", package = "certusR"))
```

Two places where the implementation departs from a published procedure, both
documented in the help and in the derivation note:

* The error-expansion coefficients, as described above.
* `estimate_error_rates()` defaults to using every supplied pair. Retaining only
  pairs with a single typing error, as described in the source, discards exactly
  the pairs carrying the most errors; at a true rate of 10% that procedure
  returns about 4.8% while the default returns about 10.5%.

## References

Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998) Statistical
confidence for likelihood-based paternity inference in natural populations.
*Molecular Ecology* **7**, 639–655. doi:10.1046/j.1365-294x.1998.00374.x

Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the computer
program CERVUS accommodates genotyping error increases success in paternity
assignment. *Molecular Ecology* **16**, 1099–1106.
doi:10.1111/j.1365-294X.2007.03089.x (corrigendum: *Molecular Ecology* **19**,
1512, doi:10.1111/j.1365-294x.2010.04544.x)

Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
(2025) Pairwise paternity assignment with forward–backward simulations: refining
CERVUS using trio-based likelihood and locus-specific error rates. *Ecology and
Evolution* **15**(10), e72230. doi:10.1002/ece3.72230

## Licence

MIT.
