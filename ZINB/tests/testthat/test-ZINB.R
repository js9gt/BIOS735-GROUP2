test_that("errors for bad input", {
  expect_error(fit.zinb(1, c(1, 2, 3, 4, 0, 0, 0, 0)))
  expect_error(LRT1D(1, c(1, 2, 3, 4, 0, 0, 0, 0)))
})
