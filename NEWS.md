# spawnR 0.9.0

Marker informativeness and non-exclusion probabilities, as CERVUS reports them.

* New `marker_summary()`. For each locus it gives the number of alleles,
  number typed, observed and unbiased expected heterozygosity, PIC and F_IS,
  plus five non-exclusion probabilities. The `combined` attribute holds the
  product over loci. The five probabilities are:
  - NE-1P: first parent, other parent unknown
  - NE-2P: second parent, first parent known
  - NE-PP: an unrelated parent pair
  - NE-I: identity of two unrelated individuals
  - NE-SI: identity of two full sibs
* The exclusion formulas are those of Jamieson & Taylor (1997). Their eq. 2a is
  NE-1P, eq. 1a is NE-2P and eq. 3a is NE-PP; the paper numbers the first two
  the other way round from CERVUS. The test suite checks every probability
  against direct enumeration over all genotypes. It also checks the paper's
  worked values: the equal-frequency maxima (eqs 1b, 2b and 3b, and Table 1)
  and the HTG4 row of Table 2 (exclusion 0.54, 0.36 and 0.73).
* Frequencies are counted with no floor by default, as CERVUS does. Pass
  `reference` to count them in the candidate parents only.
* F_IS is reported so that loci with excess heterozygosity stand out. In
  salmonids these are often paralogous and cause spurious mismatches.

# spawnR 0.8.0

Known crosses, and documentation brought up to date.

* `assign_trios()` gains `crosses`: the mating design actually used, one row
  per cross (`mother`, `father`, and optionally `weight`, e.g. eggs
  fertilised). The trio question becomes "which of the crosses made is this
  offspring from", not "which of every mother-by-father combination".
  `mothers` and `fathers` default to the parents in `crosses`.
  - **Significance from the real design.** `sim_delta(type = "pair")` gains
    `crosses` and `cross_weights`. The simulation uses the same mothers,
    fathers and pattern of who was crossed with whom, and draws the true cross
    from the list, so the runner-up is typically a half-sib cross sharing a
    parent, as in the data. `prop_sampled_crosses` is the probability that an
    offspring comes from a listed cross at all.
  - **Posterior over crosses.** Pair LODs share a denominator, so their
    differences are exact log-likelihood ratios between crosses. Each
    offspring gets `trio_posterior` (the best cross), `trio_post_none` (from
    no listed cross), and `trio_set`, the smallest set of crosses reaching
    `set_level` (0.95) posterior. The set is one cross when the data decide,
    several (usually half-sib crosses) when they do not. Weights act as a
    prior.
  - **Family sizes.** Attribute `families` gives each cross its expected
    number of offspring (summed posteriors) and the number assigned
    significantly, and `summary()` prints it.
  - **Calibration check.** 10 dams x 10 sires, each dam crossed with 3 sires,
    120 offspring, 6 replicates. With 4 four-allele loci, 82% of top crosses
    were correct (58% when scoring all combinations), and the share assigned at
    95% rose from 20% to 63% at 93.8% precision. 97.9% of 95% cross sets
    contained the true cross, with a mean set size of 1.9. With 6 five-allele
    loci, 98% were correct against 92%.
* Printing a column subset of an `assign_trios()` result no longer fails.
* New `vignette("trios")`: mothers, fathers and known crosses, computed when
  the vignette is built.
* The study scripts in `inst/studies/` now run against the installed package.
  Five had sourced `../R`, a path that does not exist from that directory, and
  two sourced a helper file under a wrong name. All six were re-run with 0.8.0.
  Every results table in the vignettes reproduced except one cell of the
  criteria table, which changed as a result of the 0.7.1 fix (Delta recall with
  50 candidates, all sampled: 99.8% to 100%). Each vignette now names the
  script behind its tables.
* README: install from GitHub, the input format (including the 0/1 SNP coding
  trap), known crosses, the posterior criterion, what a critical Delta of 0
  means, the vignette list, and re-measured timings.

# spawnR 0.7.2

The significance simulations now follow the same rules as the analysis they
calibrate.

* **Screens are applied inside the simulation.** `max_mismatch` dropped
  candidates from the real data before Delta was computed, but the simulation
  in `sim_delta()` never dropped anyone, so observed and simulated Deltas were
  computed over different pools. `sim_delta()` gains `max_mismatch` and
  `min_loci`. `assign_trios()` and `assign_parentage(criterion = "delta")` pass
  their screens through, so a true parent removed by the screen counts against
  the criterion, as it does in the real analysis. In a check with 3% error and
  `max_mismatch = 1`, the 0.7.1 results were already close to calibrated and
  this changed them little. The point is that the two now agree by
  construction rather than by luck.
* **New `min_loci` screen** in `parentage_lod()`, `assign_parentage()`,
  `assign_trios()` and `sim_delta()`. A comparison scored at fewer than
  `min_loci` loci is dropped, so an offspring typed at too few loci is left
  unassigned instead of being judged against a threshold simulated at average
  typing. CERVUS applies a rule of this kind; a common choice is half the
  loci. The default `NULL` applies no screen.
* **`prop_typed` defaults to the individuals analysed.** It was the typing rate
  over every individual in the genotype object, including ones not in the
  analysis. `assign_trios()` now uses the offspring, mothers and fathers;
  `assign_parentage()` uses the offspring, candidates and known parents. This
  shifts results only when the genotype object holds individuals with typing
  rates different from those analysed.
* **Exclusions at zero error are handled the same way in simulation and
  scoring.** With `error = 0`, an outright Mendelian exclusion was scored as
  -700 in real data, but the simulations in `sim_delta()`, `lod_null()` and
  `pool_relatedness()` treated it as a locus not compared. All now use -700.
  Results with `error > 0` (and no allele frequency of exactly 0) are
  unchanged.

# spawnR 0.7.1

Fixes a calibration bug in the Delta criterion.

* **`delta_critical()` no longer turns a perfect simulation into a high bar.**
  When every simulated assignment was correct, the critical Delta came back as
  the *smallest simulated Delta* rather than 0. The simulation had seen no wrong
  assignment at any Delta, yet any real assignment below that minimum was called
  non-significant. The tell-tale sign was an identical `critical_delta` at the
  80% and 95% levels. The critical value is now 0 whenever the simulated
  assignments taken all together attain the level, matching the 0.00 CERVUS
  reports. Otherwise the threshold is unchanged: the smallest Delta at which
  the cumulative proportion correct still meets the level, subject to `min_n`.
  The result is still `NA` when the level is never reached.

  **Who is affected:** anyone using `delta_critical()` directly, or through
  `assign_trios()` or `assign_parentage(criterion = "delta")`, with a panel
  strong enough that the simulation contains no wrong assignments. Typically
  that means hundreds of SNPs or many microsatellites, especially with small
  candidate pools. In these cases the 0.7.0 results understated significance,
  sometimes badly. On a 2-dam x 4-sire test where every offspring was assigned
  correctly in simulation, only 79-85% were called significant, and the trio
  critical Delta moved between about 14 and 65 depending on the locus set.
  Re-run those analyses. Calls made at a positive, located threshold are
  unaffected.

* `delta_critical()` gains a `resolved` column: `TRUE` for a positive threshold
  located within the simulated Deltas, `FALSE` when the critical value is 0
  because the whole simulation already attains the level, and `NA` when the
  level is not reached. A `FALSE` means "any assignment with a positive-LOD
  contender passes", not that a threshold was found. The column is carried
  into the `criteria` attribute of `assign_trios()` and `assign_parentage()`.
* `allele_freqs()`: the default `min_freq` is now `min(1/(2n), 0.01)` rather
  than `1/(2n)`, for n typed individuals at a locus. With small reference sets
  the old floor was large: with six reference parents every allele was floored
  at about 0.083, so an allele never seen in the reference counted as common as
  one seen once. That distorted LOD scores and the Delta simulation, which
  draws genotypes from these frequencies. Results change only where fewer than
  50 individuals are typed at a locus; pass `min_freq = 1 / (2 * n)` to
  reproduce 0.7.0. A message now flags loci with fewer than 25 typed reference
  individuals, since frequencies from so few allele copies are noisy whatever
  the floor; an explicit `min_freq` silences it.
* Unit tests now use testthat (`tests/testthat/`), alongside the existing
  `tests/test-spawnR.R` script. testthat is added to Suggests.

# spawnR 0.7.0

Mothers, fathers and trios in one call, each with a significance test.

* New `assign_trios()`. Given offspring, candidate mothers and candidate
  fathers, it returns one row per offspring with the most likely mother, the
  most likely father, and the most likely mother-father pair chosen on the joint
  trio LOD. For each of the three it reports the LOD, Delta, loci compared, the
  Delta confidence level attained and whether it is significant. It also counts
  mismatching loci (exclusion alone, ignoring error) for each parent separately,
  and for the trio jointly.
  `output = "all"` returns every mother-by-father combination per offspring,
  ranked, with both single-parent LODs and all three mismatch counts.
* `sim_delta()` gains `type = "pair"`, with `n_candidates2` and
  `prop_sampled2` for the mother pool. Until now, Delta significance existed only
  for single parents. A simulated pair counts as correct only when both parents
  are right, so the critical Delta is a statement about the whole trio.
  Replicates are chunked to bound memory, and 1000 replicates against 60 x 60
  candidate pairs at 12 loci take about 8 s.
* Calibration check: 6 replicate populations, 7 six-allele loci, 70% of each
  sex sampled. Observed precision at the stated levels was 81.0% / 95.5% for
  mothers, 81.5% / 98.1% for fathers, and 80.5% / 94.4% for trios, at 80% / 95%
  respectively.

# spawnR 0.6.0

External validation against the published literature and against CERVUS itself.

* Reproduces the headline simulation result of Kalinowski et al. (2007). From
  their stated parameters alone - an 11-allele frequency spectrum, six loci, 100
  candidate males, 1% error, Delta 0.99 - this package assigns **72.1%** of
  offspring against their reported **73%**, with precision at exactly 99.0%. The
  re-derived error coefficients therefore recover what CERVUS implements, which
  is what the paper's corrigendum says the printed equations failed to convey.
  Script in `inst/studies/06-kalinowski-replication.R`.
* `read_cervus_results()` reads a CERVUS parentage export. Columns are detected
  by pattern and the mapping is printed, since CERVUS's header text varies by
  analysis type and version and no authoritative list is published.
* `compare_cervus()` joins this package's scores to CERVUS's and reports matched
  rows, correlation, mean and maximum absolute difference, rows over tolerance,
  and whether the mismatch and loci-compared counts agree.
* `vignette("validation")` documents both, including the four inputs that must
  match before a difference means anything, how to read a disagreement, and what
  is expected *not* to match.
* `genotypes()` now also treats `*` and `-9` as missing allele codes.
* Honest gap: the 51% figure Kalinowski et al. report for Marshall et al.'s
  *original* equations is not reproduced. That needs the 1998 error formula as
  printed, which is closed-access and could not be obtained; only its defect is
  documented, not its form. Guessing at it would not be a replication.

# spawnR 0.5.0

Extends the relatedness-aware alternative to the joint parent-pair
configuration, which turns out to have been the most exposed of the three.

* `parentage_lod(type = "pair")` now accepts `relatedness` and `ibd`. Under H2
  neither alleged parent is a parent, but each may stand in for a relative of
  the one it is tested as, and the offspring has only two alleles to share, so
  the two sides enter jointly:
  `Pr(o | am, af) = ab T(o|am,af) + a(1-b) T(o|am) + (1-a)b T(o|af) + (1-a)(1-b) P(o)`.
  Setting `a = b = 0` recovers the classical `P(o)` and `a = b = 1` recovers the
  numerator, so the LOD goes to zero as it must.
* `ibd` for a pair takes one vector for both alleged parents, or a two-row
  matrix or list of two to treat the sexes differently. `relatedness` lets each
  alleged parent draw its class independently, mixing over every pair of classes.
* Simulation with full sibs of both the true sire and the true dam: with both
  true parents present, nothing is wrong (1.000 observed). With them withheld
  half the time, a stated 0.997 became an observed 0.653, about twice the worst
  single-parent gap, because a sib-by-sib pair mimics the true pair on both
  sides at once. The correction recovers that to 0.946, dropping the assignment
  rate from 75% to 50%.
* New cross-configuration test: a pair whose alleged mother is fixed as the true
  parent and whose alleged father is unrelated now provably collapses onto the
  one-known formula.

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
