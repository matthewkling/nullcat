test_that("curvecat works", {

      m1 <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)
      m2 <- curvecat(m1, 1000)

      # check that matrix is actually perturbed
      expect_failure(expect_equal(m1, m2))

      # check that row multisets are preserved
      expect_equal(as.vector(apply(m1, 1, sort)),
                   as.vector(apply(m2, 1, sort)))

      # check that col multisets are preserved
      expect_equal(as.vector(apply(m1, 2, sort)),
                   as.vector(apply(m2, 2, sort)))
})
