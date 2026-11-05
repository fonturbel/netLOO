# Regression tests for the fixes in netLOO 0.1.0

ex_web <- function(y) read_network_matrix(
  system.file("extdata", sprintf("data%s_clean.csv", y), package = "netLOO"))

test_that("calc_H2 matches bipartite on a weighted web and is NA for binary", {
  skip_if_not_installed("bipartite")
  utils::data("Safariland", package = "bipartite", envir = environment())
  expect_equal(unname(calc_H2(Safariland)["H2"]), 0.8537, tolerance = 1e-3)
  expect_warning(h <- calc_H2(ex_web("2015")), "not defined for binary")
  expect_true(is.na(h["H2"]))
})

test_that("H2_fixed is deprecated and no longer returns 1 - H2'", {
  skip_if_not_installed("bipartite")
  utils::data("Safariland", package = "bipartite", envir = environment())
  expect_warning(h <- H2_fixed(Safariland), "deprecated")
  expect_equal(unname(h["H2"]), 0.8537, tolerance = 1e-3)
})

test_that("calc_NODFc is on the right scale and matches maxnodf::NODFc", {
  skip_if_not_installed("maxnodf"); skip_if_not_installed("vegan")
  w <- ex_web("2017")
  nc <- calc_NODFc(w, quality = 0)
  expect_true(nc$NODF_n > 0 && nc$NODF_n <= 1)
  ref <- suppressWarnings(maxnodf::NODFc(w, quality = 0))
  expect_equal(nc$NODF_c, ref, tolerance = 1e-8)
})

test_that("calculate_network_metrics fills generality/vulnerability and NODFc", {
  skip_if_not_installed("bipartite"); skip_if_not_installed("maxnodf")
  w <- ex_web("2015")
  m <- calculate_network_metrics(w, "2015", robustness_reps = 5, nodfc_quality = 0)
  expect_false(is.na(m$generality))
  expect_false(is.na(m$vulnerability))
  expect_false(is.na(m$NODFc))
  expect_true(is.na(m$H2))
  expect_equal(m$links_per_species, sum(w) / sum(dim(w)))
})

test_that("binary null webs keep links and marginal totals", {
  skip_if_not_installed("vegan")
  w <- ex_web("2015")
  nw <- generate_null_webs(w, 20, "quasiswap")
  expect_true(all(vapply(nw, function(x) sum(x > 0), numeric(1)) == sum(w > 0)))
  expect_true(all(vapply(nw, function(x) all(rowSums(x) == rowSums(w)) &&
                           all(colSums(x) == colSums(w)), logical(1))))
  expect_warning(generate_null_webs(w, 2, 1), "weighted")
})

test_that("null models are reproducible and give empirical p-values", {
  skip_if_not_installed("bipartite"); skip_if_not_installed("vegan")
  w <- ex_web("2017")
  r1 <- enhanced_null_models(w, n_null = 30, n_cores = 1, seed = 42)
  r2 <- enhanced_null_models(w, n_null = 30, n_cores = 1, seed = 42)
  expect_identical(r1, r2)
  expect_equal(r1$connectance$note, "fixed by null model")
  expect_true(is.na(r1$connectance$p_value))
  expect_true(r1$nestedness$p_value > 0 && r1$nestedness$p_value <= 1)
})

test_that("build_metaweb is binary by default and check_metaweb detects drift", {
  nets <- list(a = ex_web("2015"), b = ex_web("2017"))
  mw <- build_metaweb(nets)
  expect_true(is_binary_web(mw))
  expect_true(max(build_metaweb(list(nets$a, nets$a), binary = FALSE)) == 2)
  expect_null(check_metaweb(mw, nets))
  mw2 <- mw; mw2[1, 1] <- 1 - mw2[1, 1]
  expect_equal(check_metaweb(mw2, nets)$mismatched_cells, 1)
})

test_that("contribution keys do not collide and keep levels apart", {
  a <- matrix(1, 1, 1, dimnames = list("x_y", "z"))
  b <- matrix(1, 1, 1, dimnames = list("x", "y_z"))
  expect_equal(calc_contribution(a, b)$shared_interactions, 0)
  s <- matrix(1, 1, 1, dimnames = list("Same", "Same"))
  expect_equal(calc_contribution(s, s)$n_species, 2)
})

test_that("read_network_matrix trims names and warns on duplicates", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(";v1;v1 ", "p1;1;0", " p2;0;1"), f)
  expect_warning(m <- read_network_matrix(f), "duplicated column")
  expect_equal(rownames(m), c("p1", "p2"))
})

test_that("trend tests report the minimum attainable p-value", {
  expect_message(tr <- temporal_trend_tests(data.frame(x = c(1, 3, 2, 4))),
                 "descriptive")
  expect_equal(tr$x$min_possible_p, 2 / 24)
})

test_that("weight_scale is deprecated in pipeline_config", {
  expect_warning(pipeline_config(weight_scale = 5), "deprecated")
})
