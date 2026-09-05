# spawnR 0.4.0

Handles relatives in the candidate pool, the last known way the posterior could
be confidently wrong.

* New `relatedness` argument to `parentage_lod()`: the alternative hypothesis
  becomes a mixture over relationship classes rather than "unrelated stranger".
  Classes are mixed over *multilocus* likelihoods, `LOD = -log(sum_c w_c
  exp(-L_c))`, because a candidate either is or is not a sib for its whole
  genome. Mixing inside the per-locus denominator instead - the obvious
  approach, and exact algebra - shifts parents and their sibs almost equally and
  recovers almost none of the lost calibration. Both are measured in the new
  vignette.
* `ibd_mixture()` builds the coefficients from named relationship classes.
* `lod_locus()` and `parentage_lod()` gain a low-level `ibd` argument for a
  single relationship class.
* `pool_relatedness()` screens a pool for excess allele sharing, reported as a
  lower bound rather than an estimate.
* `vignette("relatedness")` shows that the failure needs two conditions at once,
  relatives present *and* the true parent unsampled: either alone is harmless,
  together they take a stated 0.993 down to 0.838 observed. The multilocus
  mixture restores calibration to within a point or two everywhere except the
  hardest cell, where it improves 0.838 to 0.907 and cuts the assignment rate to
  match what the markers can actually support.
* Not implemented for `type = "pair"`.

# spawnR 0.3.0

Adds a posterior criterion, which resolves the trade-off between the two
threshold rules, and two vignettes documenting the problem and the fix.

* `parentage_posterior()` returns the probability that each candidate is the
  parent, given the LOD scores, the size of the pool searched, and a prior
  `prop_sampled` that the true parent is in it. Trio-specific and pool-aware at
  once. Computed on the log scale.
* `parentage_fdr()` chooses a posterior cutoff controlling the expected false
  discovery rate, using the fact that calibrated posteriors make the expected
  number of errors a sum of `1 - posterior`.
* `assign_parentage(criterion = "posterior")` is now the default. It requires
  `prop_sampled`, which is a modelling assumption rather than a default.
* `vignette("criteria")` documents, by factorial simulation, that the
  trio-specific rule holds recall near-constant while precision falls to 54%
  over pools of 50 to 5000 and sampling fractions of 1 to 0.5, and that the
  Delta rule fails in exactly the mirror direction, holding precision while
  recall falls to 54%.
* `vignette("posterior")` shows the posterior is calibrated across that whole
  design, holds precision between 98.5% and 100% in every cell, costs at most
  2.3 points of precision under a prior misspecified by 0.4, and controls FDR
  to its target.
* Simulation scripts behind both vignettes ship in `inst/studies/`.

# spawnR

Renamed from `parentageLR` (briefly `certusR`). To spawn is to shed gametes
into open water — reproduction that cannot be observed, which is why genetic
parentage assignment exists. It also reads as spawning an R process. No code
changed in the rename.

# spawnR 0.2.0

Performance. Results are unchanged; the test suite holds the new engine against
a locus-by-locus reference implementation and the new null samplers against the
allele-level ones they replace.

* `parentage_lod()` now evaluates the likelihood once per distinct genotype at
  each locus rather than once per candidate, and reduces the
  offspring-by-candidate score to a table lookup. Accumulation switches between
  dense matrix products (few distinct genotypes, i.e. SNP panels) and per-locus
  gathers (microsatellites). A 300-offspring by 300-candidate by 300-SNP problem
  went from 17.2 s to 0.5 s.
* `parentage_lod(output = "matrix")` returns the score matrices directly. For
  large candidate sets, assembling and sorting the long-format table costs more
  than computing the scores.
* `lod_null()` draws each locus contribution from its exact discrete
  distribution instead of simulating alleles: 3-10x faster, same distribution.
* `sim_delta()` uses the same per-genotype reduction, roughly 2x faster.
* New: "Choosing a criterion" in `?assign_parentage`, documenting that the
  trio-specific criterion controls the risk of rejecting a true parent but not
  the false-positive rate, which grows with the candidate pool.
* `inst/optional-cpp/` carries an Rcpp version of the accumulation kernel, with
  benchmarks, for anyone who later needs it. The package itself stays
  dependency-free and has no compiled code.

# spawnR 0.1.0

First release.

* LOD scores for one-parent-known, both-parents-unknown and joint parent-pair
  configurations, under a random-genotype-replacement error model with a
  separate rate per locus.
* Trio-specific null distributions by forward and backward simulation, with
  lower-tail critical values.
* Population-level Delta criterion with relaxed and strict confidence levels.
* Locus-specific error-rate estimation from known parent-offspring pairs.
* Clean-room reconstruction from the primary literature; see
  `inst/notes/derivation.md`.
