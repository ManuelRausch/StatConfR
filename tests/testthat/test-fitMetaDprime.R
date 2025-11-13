test_that("fitMetaDprime handles non-numeric participant IDs", {
  dummy_data <- data.frame(
    rating = factor(c("low", "high", "low", "high")),
    stimulus = factor(c("S1", "S1", "S2", "S2")),
    correct = c(1, 0, 1, 0),
    participant = factor(c("p1", "p1", "p2", "p2"))
  )

  mocked_int_fit <- function(...) {
    data.frame(
      dprime = 1,
      c = 0,
      metaD = 1,
      Ratio = 1
    )
  }

  res <- with_mocked_bindings(
    fitMetaDprime(dummy_data, model = "ML", nInits = 1, nRestart = 1, .parallel = FALSE),
    int_fitMetaDprime = mocked_int_fit
  )

  expect_s3_class(res, "data.frame")
  expect_equal(unique(res$model), "ML")
  expect_equal(nrow(res), 2L)
  expect_equal(sort(unique(res$participant)), c("p1", "p2"))
})
