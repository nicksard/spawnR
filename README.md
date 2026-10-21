# spawnR

Likelihood-ratio ("LOD score") parentage inference for codominant markers, with
locus-specific genotyping error and simulation-based confidence. Base R only —
no dependencies.

## The name

To *spawn*, of an aquatic animal, is to shed eggs or milt into open water —
from Latin *expandere*, to spread out. Spawning is reproduction you cannot
watch: gametes are broadcast, several males may fertilise one female's clutch
at once, and no pair is visible from the bank. That unobservability is exactly
why genetic parentage assignment exists in these systems.

It reads the other way too. To `spawn` is what a program does when it starts a
child process. Both readings are true, which is the point.

## What it does

LOD scores in the three configurations of Kalinowski et al. (2007):

| configuration     | `type`           | hypothesis pair |
|-------------------|------------------|-----------------|
| one parent known  | `"one_known"`    | candidate is the second parent, vs. unrelated |
| both unknown      | `"both_unknown"` | candidate is a parent, vs. unrelated |
| both jointly      | `"pair"`         | candidate pair are the parents, vs. both unrelated |

Three ways to attach confidence:

* **Trio-specific** (Amiri Roudbar et al. 2025) — `lod_null()` builds a
  reference distribution by forward or backward simulation conditioned on the
  actual genotypes and missing-data pattern of the case being tested;
  `lod_critical()` takes the lower-tail quantile.
* **Population-level** (Marshall et al. 1998) — `sim_delta()` and
  `delta_critical()` reproduce the Delta criterion that CERVUS reports, at the
  conventional relaxed (80%) and strict (95%) levels, for single parents and
  (`type = "pair"`) for mother-father pairs, optionally under a known mating
  design.
* **Posterior** — `parentage_posterior()` gives the probability that the named
  candidate is the parent, given the pool size and the proportion of parents
  sampled; `parentage_fdr()` picks a cutoff for a target false discovery rate.

`assign_parentage()` runs a complete single-parent analysis under any of the
three. `assign_trios()` does mothers, fathers and mother-father pairs together,
and with known crosses scores only the crosses actually made.

Per-locus genotyping error rates can be estimated from known parent–offspring
pairs with `estimate_error_rates()`.

## Install

```r
# install.packages("remotes")
remotes::install_github("nicksard/spawnR", build_vignettes = TRUE)
```

`build_vignettes = TRUE` is needed for `vignette()` to find the vignettes below
(it needs the knitr and rmarkdown packages).

Base R only (R >= 4.0); nothing else is installed.

## Input

`genotypes()` takes one row per individual: an identifier column, then two
adjacent columns per locus holding the two alleles. Allele labels can be any
text or numbers.

```
id    L1_a L1_b L2_a L2_b
F01   152  156  203  203
M07   152  160  NA   NA
```

By default the codes `""`, `"0"`, `"000"`, `"NA"`, `"*"` and `"-9"` mean
missing, and a locus with either allele missing is treated as untyped. **SNP
data coded 0/1, as in a VCF, must be recoded first** (for example 0 to `"A"`,
1 to `"B"`), or every reference allele will be read as missing.

## Usage

```r
library(spawnR)

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

### Mothers, fathers and trios together

When candidate mothers and fathers are sampled separately, `assign_trios()`
gives one row per offspring: the best mother, the best father and the best
mother-father pair (scored jointly), each with LOD, Delta, the Delta confidence
attained, a significance flag, and the number of loci at which that parent is
excluded outright.

```r
a <- assign_trios(g, offspring = offs, mothers = females, fathers = males,
                  freqs = p, error = er$error,
                  prop_sampled_mothers = 0.8, prop_sampled_fathers = 0.6)
a            # best mother / father / trio per offspring
summary(a)   # counts by confidence level, critical Delta values

# every candidate pair for every offspring, ranked
assign_trios(g, offs, females, males, p, er$error, output = "all",
             max_mismatch = 3)
```

`max_mismatch` and `min_loci` screen out candidates excluded at too many loci
or compared at too few, in the data and in the significance simulation alike.

### Known crosses

In a hatchery or breeding programme the crosses are written down. Give them to
`assign_trios()` and the trio question becomes *which of the crosses actually
made* produced each offspring:

```r
crosses <- data.frame(mother = c("D1", "D1", "D2", "D2"),
                      father = c("S1", "S2", "S2", "S3"))
ac <- assign_trios(g, offs, freqs = p, error = er$error, crosses = crosses)
ac[, c("offspring", "trio_mother", "trio_father", "trio_confidence",
       "trio_posterior", "trio_set")]
attr(ac, "families")   # expected offspring per cross
```

The significance simulation reproduces the real design, each offspring gets a
posterior probability for each cross, and `trio_set` lists the smallest set of
crosses holding 95% of it (one cross when the data decide, usually a few
half-sib crosses when they do not). Summed over offspring, the posteriors
estimate family sizes. See `vignette("trios")`.

`simulate_population()` generates a population with a known pedigree for
calibration and testing.

## Performance

Base R throughout, no compiled code. The engine evaluates the likelihood once
per distinct genotype per locus - a locus with *k* alleles has only *k(k+1)/2*
of them - and reduces the offspring-by-candidate score to a table lookup,
accumulated either with dense matrix products (SNP panels) or per-locus gathers
(microsatellites).

Measured on a single core of a cloud container, spawnR 0.8.0, default
simulation sizes:

| task                                                         | time   |
|--------------------------------------------------------------|--------|
| `parentage_lod`, 300 offspring x 300 cand x 300 SNPs          | 0.5 s  |
| `parentage_lod`, 1000 x 2000 x 300 SNPs (matrix output)       | 4.9 s  |
| `parentage_lod`, 500 x 500 x 15 microsatellites               | 0.5 s  |
| `lod_null`, 10000 replicates, 300 SNPs                        | 0.14 s |
| `assign_parentage` pairwise, 500 x 500 x 15 msats             | 7.8 s  |
| `assign_trios`, 500 offspring, 50 dams x 50 sires, 15 msats   | 40 s   |
| `assign_trios`, 500 offspring, 150 known crosses, 15 msats    | 35 s   |

Most of the `assign_trios` time is the significance simulation; `nsim` and
`nsim_trio` trade precision of the critical values for speed.

For anything larger, `output = "matrix"` avoids building a long-format table of
millions of rows, which by then costs more than the arithmetic.

An Rcpp version of the accumulation kernel is in `inst/optional-cpp/` with
benchmarks. It is roughly 4x faster than the matrix-product path on SNP panels,
which is well under 2x end to end, and it would cost every user a compiler. The
package stays dependency-free unless profiling on real data says otherwise.

## Choosing a confidence criterion

The criteria control different things. The trio-specific `"pairwise"` rule
bounds the chance of *rejecting a true parent*; it does not bound the chance
that an unrelated candidate clears the same threshold, and that grows with the
pool. At a nominal 99% level with 15 microsatellites and the true sire always
present, the share of assignments that were wrong was 0% with 100 or 500
candidates, 3.1% with 2000, and 5.6% with 5000. The Marshall `"delta"` rule
defines confidence as the proportion of assignments that are correct and
simulates the whole pool, so it is the one to report when screening broadly.
The posterior answers the question most readers actually ask, the probability
that this assignment is right, and stayed calibrated across pool sizes from 50
to 5000. See `?assign_parentage` and `vignette("criteria")`.

When the Delta simulation contains no wrong assignments at all, as it often
does with a strong panel, the critical Delta is 0 and `delta_critical()` marks
it `resolved = FALSE`: any assignment with a positive-LOD candidate passes,
because the simulation found no threshold below which assignments fail.

## Vignettes

* `vignette("trios")` — mothers, fathers and known crosses
* `vignette("criteria")` — what the confidence criteria actually control
* `vignette("posterior")` — posterior probabilities and false discovery rates
* `vignette("relatedness")` — relatives in the candidate pool
* `vignette("validation")` — replicating Kalinowski et al. (2007), and comparing
  with CERVUS

The results tables in the last four come from the scripts in `inst/studies/`,
last re-run with spawnR 0.8.0.


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
file.show(system.file("notes", "derivation.md", package = "spawnR"))
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
