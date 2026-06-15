## Concurrency / fork-safety regression tests for the C++ layer (see misc/proposal.md §12).

test_that("get_nearest_neighbors is deterministic across n_cores (fork-safe threading)", {
  ## The hitting/commute-time kNN must give identical results regardless of thread count:
  ## each parallel iteration writes a disjoint index, so dropping the `omp critical` and using
  ## `if(n_cores>1) num_threads(n_cores)` must not change the output.
  set.seed(123)
  nv <- 60L; k <- 6L
  adj <- lapply(seq_len(nv), function(i) as.integer(sort(sample(setdiff(0:(nv - 1L), i - 1L), k))))
  trans <- lapply(adj, function(nb) rep(1 / length(nb), length(nb)))

  r1 <- conos:::get_nearest_neighbors(adj, trans, n_cores = 1L, verbose = FALSE)
  r2 <- conos:::get_nearest_neighbors(adj, trans, n_cores = 2L, verbose = FALSE)
  expect_identical(r1, r2)
})

test_that("adjustedRandcpp returns correct values for known partitions", {
  rand <- function(cl1, cl2, flag) {
    cl1u <- unique(cl1); cl2u <- unique(cl2)
    conos:::adjustedRandcpp(as.integer(cl1), as.integer(cl1u), as.integer(cl2), as.integer(cl2u),
                            length(cl1u), length(cl2u), length(cl1), as.integer(flag))
  }
  cl <- rep(1:3, each = 20)
  ## identical partitions -> 1 for Rand (flag 1), Hubert-Arabie (2), Jaccard (5)
  expect_equal(rand(cl, cl, 1L), 1)
  expect_equal(rand(cl, cl, 2L), 1)
  expect_equal(rand(cl, cl, 5L), 1)
  ## a perturbed partition is less than perfectly concordant
  cl2 <- cl; cl2[1:10] <- 2L
  expect_lt(rand(cl, cl2, 2L), 1)
})
