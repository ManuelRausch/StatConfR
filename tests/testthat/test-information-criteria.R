test_that("int_fitSDT reports AICc with +1 correction", {
  N_SA_RA <- matrix(c(20, 15), nrow = 1)
  N_SA_RB <- matrix(c(10, 5), nrow = 1)
  N_SB_RA <- matrix(c(12, 8), nrow = 1)
  N_SB_RB <- matrix(c(18, 22), nrow = 1)

  total_trials <- sum(N_SA_RA + N_SA_RB + N_SB_RA + N_SB_RB)

  fit_sdt <- getFromNamespace("fitSDT", "statConfR")

  res <- fit_sdt(
    N_SA_RA = N_SA_RA,
    N_SA_RB = N_SA_RB,
    N_SB_RA = N_SB_RA,
    N_SB_RB = N_SB_RB,
    nInits = 1,
    nRestart = 1,
    nRatings = 2,
    nCond = 1,
    nTrials = total_trials
  )

  expect_s3_class(res, "data.frame")
  expect_equal(res$AIC, 2 * res$negLogLik + 2 * res$k)

  denom <- res$N - res$k - 1
  skip_if(denom <= 0, "AICc undefined when N <= k + 1")
  expected_AICc <- res$AIC + (2 * res$k * (res$k + 1)) / denom
  expect_equal(res$AICc, expected_AICc)
})
