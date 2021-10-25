

test_that("Detection under alternative hypothesis, one DE", {
  expect_equal(DArand(matrix(c(rep(100,6), rep(1,6), rep(1,12*(10-1))), byrow = FALSE, nrow=12)+
         matrix(rnorm(120,mean=0, sd=0.1), byrow = FALSE, nrow=12),n1=6, r=100,set.seed = 19811101),1)
})
test_that("Detection under alternative hypothesis, two DE", {
  expect_equal(DArand(matrix(c(rep(c(rep(150,6), rep(1,6)),2), rep(1,12*(10-2))), byrow = FALSE, nrow=12)+
                        matrix(rnorm(120,mean=0, sd=0.1), byrow = FALSE, nrow=12),n1=6, r=100, set.seed = 19811101),1:2)
})

test_that("Detection under null hypothesis, zero DE", {
expect_equal(DArand(matrix(rnorm(132,mean=1, sd=0.1), byrow = FALSE, nrow=12),n1=6, r=100,set.seed = 19811101),integer(0))
})
