test_that("test udf", {

  udf <- function(x, y, ...) {
    x + y
  }

  expect_error(rxode2parse("b <- udf(x, y)"))

  udf <- function(x, y) {
    x + y
  }

  expect_error(rxode2parse("b <- udf(x, y)"), NA)

  expect_error(rxode2parse("b <- udf(x, y, z)"))


})
