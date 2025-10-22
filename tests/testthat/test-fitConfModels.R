test_that("fitConfModels handles non-numeric participant IDs", {
  dummy_data <- data.frame(
    diffCond = factor(c("hard", "hard", "easy", "easy")),
    rating = factor(c("low", "high", "low", "high")),
    stimulus = factor(c("S1", "S2", "S1", "S2")),
    correct = c(1, 0, 1, 0),
    participant = factor(c("subA", "subA", "subB", "subB"))
  )

  mocked_fit <- function(...) {
    data.frame(
      negLogLik = 0,
      N = 10,
      k = 2,
      BIC = 0,
      AICc = 0,
      AIC = 0,
      d_1 = 1,
      c = 0,
      theta_minus.1 = -1,
      theta_plus.1 = 1
    )
  }

  res <- with_mocked_bindings(
    fitConfModels(dummy_data, models = "SDT", nInits = 1, nRestart = 1, .parallel = FALSE),
    fitConf = mocked_fit
  )

  expect_s3_class(res, "data.frame")
  expect_equal(unique(res$model), "SDT")
  expect_equal(nrow(res), 2L)
  expect_equal(sort(unique(res$participant)), c("subA", "subB"))
})
