# Optional compiled accumulation kernel

`spawnR` is deliberately dependency-free and contains no compiled code, so
it installs anywhere R does, with no toolchain and no Rtools on Windows. The
pure-R engine reduces the per-candidate likelihood to a table lookup and then
accumulates it either with dense matrix products (small numbers of distinct
genotypes, i.e. SNP panels) or with per-locus gathers (microsatellites). That
gets within a small factor of hand-written C++.

`lod_accumulate.cpp` is the same accumulation written as an Rcpp kernel, kept
here in case the pure-R engine ever stops being fast enough. Measured against
the shipped engine on this machine, accumulating a LOD matrix only:

| problem                                  | R (BLAS path) | Rcpp  | gain  |
|------------------------------------------|---------------|-------|-------|
| 1000 offspring x 2000 cand x 300 SNPs    | 1.38 s        | 0.36 s| 3.9x  |
| 500 x 5000 x 500 SNPs                    | 2.70 s        | 0.66 s| 4.1x  |
| 2000 x 5000 x 96 SNPs                    | 2.92 s        | 0.63 s| 4.6x  |
| 500 x 500 x 20 microsatellites           | 0.19 s        | 0.004 s| 48x  |

Note that in the shipped package this accumulation is no longer the dominant
cost - assembling the long-format result table and sorting it costs as much
again - so a 4x kernel speedup is well under 2x end to end. The microsatellite
row is a large ratio on an already negligible absolute time.

To try it:

```r
Rcpp::sourceCpp(system.file("optional-cpp", "lod_accumulate.cpp",
                            package = "spawnR"))
# lod_acc2(ocode, ccode, as.numeric(LT), G)
```

Adopting it properly would mean adding `Rcpp` to `LinkingTo`/`Imports`, a
`src/` directory, and a compiler requirement for every user installing from
source. That trade is only worth making if profiling on real data shows the
accumulation dominating.
