
test_that("nullcat works for all algos", {

      set.seed(123)
      m1 <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)

      for(method in nullcat_methods()){

            # no error for default "category" or "index" modes
            expect_no_error(nullcat(m1, method, 1000, output = "index"))
            expect_no_error(m2 <- nullcat(m1, method, 1000))

            # check that matrix is actually perturbed
            expect_failure(expect_equal(m1, m2))
      }
})


