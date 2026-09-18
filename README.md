# netLOO

### netLOO package version 0.0-1

This package contains a dataset-agnostic toolkit for **multi-temporal bipartite network analysis** that apply the full pipeline to any set of bipartite networks sampled at three or more time points (plant–pollinator webs across years, host–parasite
series, networks across an invasion gradient, etc.).

This package was based on the ideas developed by Rodrigo Medel & Francisco E. Fontúrbel's to assess changes in plant-pollinator networls in the Chilean blooming desert.  Software development was made using Claude Code (model Opus 5).

## What this package computes

- **Per network**: connectance, NODF, connectance-corrected NODF (NODF_c),
  modularity, H2′ (implementing a bug fix that caused an unexpectwed error in bipartite 2.23 bug fix), generality, vulnerability, web asymmetry, links/species, and robustness.
- **Temporal**: Spearman trend of each metric over time; species
  persistence, Jaccard turnover, gains/losses, core species.
- **Beta diversity**: `betalinkr` partitioning (S / OS / WN / ST) for every
  pair of time points.
- **Leave-one-out contribution**: for each time point, expected contribution
  to the pooled metaweb vs. observed-unique contribution relative to the
  leave-one-out metaweb, plus an efficiency index. The conceptual key aspect here is that we compute the metaweb leaving out the network to be compared in order to avoid redundancy in the analyses that may cause non-independence.
- **Null models**: parallel z-score / p-value for connectance, nestedness,
  modularity.
- **Rarefaction**: species/interaction accumulation curves.
- **Plots**: temporal trajectories, richness, turnover, LOO contribution,
  taxonomic composition, bipartite web diagrams, and betalink overlays.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fonturbel/netLOO")
```

`maxnodf` (needed only for `do_nodfc = TRUE`) is not on CRAN:

```r
remotes::install_github("christophhoeppke/maxnodf")
```

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
  n_null   = 1000,
  do_nodfc = TRUE
)
results <- run_network_pipeline(cfg)

results$metrics_df          # per-period structural metrics
results$persistence$turnover
results$beta                # betalinkr S/OS/WN/ST, all pairs
results$loo                 # leave-one-out contribution + efficiency
results$null_models         # z-score / p-value per metric

plot_richness(results)
plot_loo_contribution(results)
```

Or just call the leave-one-out pieces directly on any list of matrices:

```r
networks <- load_networks(cfg)
metaweb  <- build_metaweb(networks)
loo      <- loo_contribution_analysis(networks, metaweb)
```

## Applying to a new dataset

1. Drop the per-period matrices in a folder.
2. Write a `manifest.csv` (`label,file`) in time order.
3. (Optional) write a `taxonomy_lookup.csv` (`species,group`).
4. Point `pipeline_config()` at them. Done.

Rename trophic levels with `lower_name` / `higher_name` (e.g. `"Hosts"` /
`"Parasites"`). Set `do_nodfc = FALSE` to skip the slow `maxnodf` step.

## Citation

Fontúrbel FE, Medel, R (2026) netLOO: a comprehensive package to compare multiple networks across time periods. <https://github.com/fonturbel/netLOO>

## License

GPL (>= 3)
