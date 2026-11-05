# netLOO 0.1.0

This release fixes several bugs found in a code review (November 2026).
**Results from earlier versions for H2', NODF_c, null-model tests,
generality and vulnerability should be recomputed.** The leave-one-out
contribution analysis is unchanged except for the edge cases noted below.

## Bug fixes that change results

* `H2_fixed()` returned 1 - H2' (the formula was inverted) on a web
  multiplied by an arbitrary `weight_scale`. `bipartite::H2fun()` had no
  bug: H2' = 0 is the expected value for binary webs, where H2' carries no
  information. `H2_fixed()` is deprecated; the new `calc_H2()` uses
  `H2fun()` unchanged and returns NA (with a warning) for binary webs.
  `calculate_network_metrics()` now reports `H2 = NA` for binary webs.
* `calc_NODFc()` / `nodf_c()` divided a 0-100 NODF by a 0-1 maximum NODF and
  used the natural log, inflating NODF_c about 43-fold. It now matches
  `maxnodf::NODFc()` (0-1 scales, log10). Rankings among networks are
  unchanged.
* `enhanced_null_models()` used `bipartite::nullmodel(method = 1)`
  (Patefield) by default. On binary webs this generates weighted matrices
  with fewer links than the observed web (about 14% fewer in the bundled
  2015 example), biasing nestedness and modularity z-scores. The default is
  now `null_method = "quasiswap"` (vegan), which keeps row/column totals and
  the number of links fixed. `p_value` is now a two-sided empirical p-value;
  the normal approximation is kept as `p_value_z`.
* `calculate_network_metrics()` always returned NA for `generality` and
  `vulnerability` (networklevel() names them `generality.HL` and
  `vulnerability.LL`).
* `build_metaweb()` returned, for each pair, the number of networks in which
  it occurred instead of presence/absence. It is now binary by default
  (`binary = FALSE` restores counts). LOO results are unaffected because
  they only test `> 0`.

## Other changes

* New `seed` argument in `pipeline_config()` (default 1) and
  `enhanced_null_models()`. Null models draw one seed per null web, so
  results no longer depend on the number of cores. `n_cores = 1` runs
  without a cluster.
* New `nodfc_quality` argument (default 1). The previous hard-coded
  `quality = 2` (simulated annealing) has crashed R on some matrices.
* New `generate_null_webs()`, `check_metaweb()` and `is_binary_web()`.
* A metaweb supplied through `pipeline_config(metaweb = )` is checked
  against the union of the networks, with a warning if species or links
  differ (usually a sign of inconsistent name cleaning).
* `read_network_matrix()` trims whitespace in species names and warns about
  duplicated names.
* `calc_contribution()` tags species by level, so a name used in both levels
  counts as two species, and builds interaction keys that cannot collide
  when names contain `_`.
* `temporal_trend_tests()` reports `n` and `min_possible_p` (the smallest
  attainable p-value; 0.083 with 4 time points) and skips constant metrics.
* `pipeline_config(weight_scale = )` is deprecated and ignored.
* `plot_beta_pair()` documents that its returned metrics use the original
  betalink partition, which is not comparable with `results$beta`.
