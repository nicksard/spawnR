# certusR

Renamed from `parentageLR`. Latin *certus*, "decided, certain", from *cernere*
"to distinguish, to decide" — and one letter from *cervus*. No code changed in
the rename.

# certusR 0.2.0

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

# certusR 0.1.0

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
