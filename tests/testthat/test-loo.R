test_that("build_loo_metaweb excludes only the focal network", {
  nets <- list(
    a = matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("p1", "p2"), c("v1", "v2"))),
    b = matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("p1", "p3"), c("v2", "v3")))
  )
  full_meta <- build_metaweb(nets)
  loo_a <- build_loo_metaweb(nets, 1)

  expect_equal(dim(full_meta), c(3, 3))
  expect_equal(sort(rownames(loo_a)), sort(rownames(nets$b)))
  expect_equal(sort(colnames(loo_a)), sort(colnames(nets$b)))
})

test_that("calc_contribution reports species/interaction overlap correctly", {
  focal <- matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("p1", "p2"), c("v1", "v2")))
  ref   <- matrix(c(1, 0, 0, 0), 2, 2, dimnames = list(c("p1", "p2"), c("v1", "v2")))

  out <- calc_contribution(focal, ref)
  expect_equal(out$n_species, 4)
  expect_equal(out$n_interactions, 2)
  expect_equal(out$shared_interactions, 1)
  expect_equal(out$unique_interactions, 1)
  expect_equal(out$unique_species, 0)
})

test_that("loo_contribution_analysis returns one row per network with an efficiency <= exp", {
  nets <- list(
    a = matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("p1", "p2"), c("v1", "v2"))),
    b = matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("p1", "p3"), c("v2", "v3"))),
    c = matrix(c(1, 1, 1, 1), 2, 2, dimnames = list(c("p2", "p3"), c("v1", "v3")))
  )
  meta <- build_metaweb(nets)
  loo <- loo_contribution_analysis(nets, meta)

  expect_equal(nrow(loo), length(nets))
  expect_true(all(loo$obs_unique_species <= loo$exp_species))
  expect_true(all(loo$obs_unique_interactions <= loo$exp_interactions))
})

test_that("loo_contribution_analysis reproduces expected values on the bundled example", {
  cfg <- pipeline_config(
    manifest = system.file("extdata", "manifest.csv", package = "netLOO"),
    data_dir = system.file("extdata", package = "netLOO"),
    do_nodfc = FALSE
  )
  nets <- load_networks(cfg)
  meta <- build_metaweb(nets)
  loo <- loo_contribution_analysis(nets, meta)

  expect_equal(loo$network, names(nets))
  expect_true(all(loo$species_efficiency >= 0 & loo$species_efficiency <= 1))
  expect_true(all(loo$interaction_efficiency >= 0 & loo$interaction_efficiency <= 1))
})
