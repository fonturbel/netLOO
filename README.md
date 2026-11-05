# netLOO

[![R-CMD-check](https://github.com/fonturbel/netLOO/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fonturbel/netLOO/actions/workflows/R-CMD-check.yaml)

### netLOO package version 0.1.0

> **Upgrading from 0.0.3?** Version 0.1.0 fixes bugs that change H2′, NODF_c,
> null-model tests, generality and vulnerability. Please recompute those results;
> see [NEWS.md](NEWS.md).

This package contains a dataset-agnostic toolkit for **multi-temporal bipartite network analysis** that apply the full pipeline to any set of bipartite networks sampled at three or more time points (plant–pollinator webs across years, host–parasite
series, networks across an invasion gradient, etc.).

This package was based on the ideas developed by Rodrigo Medel & Francisco E. Fontúrbel's to assess changes in plant-pollinator networks in the Chilean blooming desert.  Software development was made using Claude Code (model Opus 5).

## What this package computes

- **Per network**: connectance, NODF, connectance-corrected NODF (NODF_c,
  same formula as `maxnodf::NODFc()`), modularity, generality,
  vulnerability, web asymmetry, links/species, and robustness. H2′ is
  computed for weighted webs only (it is NA for binary webs, where it
  carries no information).
- **Temporal**: Spearman trend of each metric over time, with the smallest
  p-value attainable for that many time points (with 4 points no trend can
  reach p < 0.05); species persistence, Jaccard turnover, gains/losses,
  core species.
- **Beta diversity**: `betalinkr` partitioning (S / OS / WN / ST) for every
  pair of time points.
- **Leave-one-out contribution**: for each time point, expected contribution
  to the pooled metaweb vs. observed-unique contribution relative to the
  leave-one-out metaweb, plus an efficiency index. The conceptual key aspect here is that we compute the metaweb leaving out the network to be compared in order to avoid redundancy in the analyses that may cause non-independence.
- **Null models**: observed vs null nestedness and modularity, with z-scores
  and two-sided empirical p-values. The default `"quasiswap"` null keeps
  row/column totals and the number of links fixed; seeds make results
  reproducible regardless of the number of cores.
- **Rarefaction**: species/interaction accumulation curves.
- **Plots**: temporal trajectories, richness, turnover, LOO contribution,
  taxonomic composition, bipartite web diagrams, and betalink overlays.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fonturbel/netLOO")
```

Two optional dependencies are not on CRAN: `maxnodf` (needed only for
`do_nodfc = TRUE`) and `betalink` (needed only for `plot_beta_pair()`;
archived from CRAN):

```r
remotes::install_github("christophhoeppke/maxnodf")
remotes::install_version("betalink", version = "2.2.1")
```

## Dependencies

Base R (`stats`, `utils`, `graphics`, `grDevices`) covers the core; the
following packages are needed for full functionality (network metrics,
beta diversity, plots, parallel null models):

```r
install.packages(c("bipartite", "vegan", "igraph", "ggplot2", "tidyr",
                    "dplyr", "foreach", "doParallel"))

# not on CRAN:
remotes::install_github("christophhoeppke/maxnodf")
remotes::install_version("betalink", version = "2.2.1")
```

| Package                 | Used for                                                                  |
| ----------------------- | ------------------------------------------------------------------------- |
| `bipartite`             | core bipartite network metrics (connectance, NODF, modularity, H2′, etc.) and `betalinkr()` beta diversity partitioning |
| `vegan`                 | diversity/ordination utilities                                            |
| `igraph`                | network representation/algorithms                                         |
| `ggplot2`               | plotting functions (`plot_richness()`, `plot_loo_contribution()`, ...)    |
| `tidyr`, `dplyr`        | data wrangling                                                            |
| `foreach`, `doParallel` | parallel null models                                                      |
| `betalink` (archived)   | shared vs. exclusive species/links plot (`plot_beta_pair()`)              |
| `maxnodf` (GitHub)      | connectance-corrected NODF (`nodf_c()`, `do_nodfc = TRUE`)                |
| `rmarkdown`, `testthat` | vignettes / package tests (dev only)                                      |

## Input format

- One CSV matrix **per time point**: rows = lower trophic level (e.g.
  plants), columns = higher trophic level (e.g. pollinators), `;`-separated,
  first column = row names. Binary or weighted. Species sets need **not**
  match across time points.
- A **manifest** CSV with columns `label,file` listing the time points in
  order (or omit it to auto-discover every CSV in `data_dir`).
- Optional **taxonomy** CSV with columns `species,group` to classify the
  higher-level species into groups.

## Quick start

```r
library(netLOO)

cfg <- pipeline_config(
  manifest = system.file("extdata", "manifest.csv", package = "netLOO"),
  data_dir = system.file("extdata", package = "netLOO"),
  n_null      = 1000,
  null_method = "quasiswap",   # fixed-fill null for binary webs
  do_nodfc    = TRUE,
  seed        = 1
)
results <- run_network_pipeline(cfg)

results$metrics_df          # per-period structural metrics
results$persistence$turnover
results$beta                # betalinkr S/OS/WN/ST, all pairs
results$loo                 # leave-one-out contribution + efficiency
results$null_models         # z-score / empirical p-value per metric

plot_richness(results)
plot_loo_contribution(results)
```

## Running a single component

Skip the full pipeline (and its slow null models) when you only need one
piece. `nodf_c()` and `loo_only()` take either a `bd_pipeline_config`, a
named list of matrices, a single matrix, or (for `nodf_c()`) a file path:

```r
# Just NODF_c, on one matrix or a whole folder of them
nodf_c(read_network_matrix("2020.csv"))
nodf_c(cfg)                    # named list, one result per time point
nodf_c(cfg, label = "2020")    # single time point

# Just the leave-one-out contribution analysis, no null models
loo_only(cfg)
loo_only(list(t1 = net1, t2 = net2, t3 = net3))
```

Both still accept a `bd_pipeline_config` so existing `cfg`-based scripts
keep working unchanged. The lower-level building blocks
(`load_networks()`, `build_metaweb()`, `calc_NODFc()`,
`loo_contribution_analysis()`, `calculate_network_metrics()`, ...) are
exported individually too, if you want finer control.

## Applying to a new dataset

1. Drop the per-period matrices in a folder.
2. Write a `manifest.csv` (`label,file`) in time order.
3. (Optional) write a `taxonomy_lookup.csv` (`species,group`).
4. Point `pipeline_config()` at them. Done.

Rename trophic levels with `lower_name` / `higher_name` (e.g. `"Hosts"` /
`"Parasites"`). Set `do_nodfc = FALSE` to skip the slow `maxnodf` step.

> **Known issue**: `maxnodf::maxnodf()` can segfault in its simulated
> annealing step (`quality = 2`) on some matrices. netLOO now defaults to
> `nodfc_quality = 1`; if you still hit a crash, use `nodfc_quality = 0` or
> `do_nodfc = FALSE`.

## Citation

Fontúrbel FE, Medel, R (2026) netLOO: a comprehensive package to compare multiple networks across time periods. <https://github.com/fonturbel/netLOO>

## License

GPL (>= 3)
