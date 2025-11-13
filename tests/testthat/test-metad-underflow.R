test_that("meta-d likelihood clamps probabilities before logging", {
  nRatings <- 3
  parameters <- log(c(1.5, 1e-9, 1e-9, 1e-9, 1e-9))
  nC_rS1 <- c(20, 0, 0)
  nI_rS1 <- c(0, 5, 0)
  nC_rS2 <- c(4, 25, 0)
  nI_rS2 <- c(0, 6, 0)

  value_ml <- getFromNamespace("negLoglMetaD", "statConfR")(
    parameters = parameters,
    nC_rS1 = nC_rS1,
    nI_rS1 = nI_rS1,
    nC_rS2 = nC_rS2,
    nI_rS2 = nI_rS2,
    nRatings = nRatings,
    cprime = 0
  )

  value_f <- getFromNamespace("negLoglFleming", "statConfR")(
    parameters = parameters,
    nC_rS1 = nC_rS1,
    nI_rS1 = nI_rS1,
    nC_rS2 = nC_rS2,
    nI_rS2 = nI_rS2,
    nRatings = nRatings,
    type1_c = 0
  )

  expect_true(is.finite(value_ml))
  expect_true(is.finite(value_f))
})
