test_that("curvecat works", {

      set.seed(123)
      m1 <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)

      expect_no_error(m2 <- curvecat(m1, 1000, output = "index"))
      expect_no_error(m2 <- curvecat(m1, 1000))

      # check that matrix is actually perturbed
      expect_failure(expect_equal(m1, m2))

      # check that row multisets are preserved
      expect_equal(as.vector(apply(m1, 1, sort)),
                   as.vector(apply(m2, 1, sort)))

      # check that col multisets are preserved
      expect_equal(as.vector(apply(m1, 2, sort)),
                   as.vector(apply(m2, 2, sort)))
})



test_that("category and index modes follow same randomization (given a shared seed)", {

      m <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)

      set.seed(123)
      m_rand <- curvecat(m, n_iter = 1000, output = "category")

      set.seed(123)
      idx <- curvecat(m, n_iter = 1000, output = "index")
      m_vec <- as.integer(m)
      m_rand2 <- matrix(m_vec[idx], nrow = nrow(m), ncol = ncol(m))

      expect_equal(m_rand, m_rand2)
})



test_that("curvecat stationary distribution matches vegan curveball given binary data", {
      skip_on_cran()

      # a small binary matrix
      set.seed(1234)
      b <- matrix(sample(0:1, 16, replace = T), nrow = 4)

      sim <- function(nullmod){
            apply(as.array(simulate(nullmod, burnin = 1e4, nsim = 1e5)), c(1, 2), mean)
      }

      cb <- sim(vegan::nullmodel(b, method = "curveball"))
      cc <- sim(vegan::nullmodel(b, method = commsim_curvecat_seq()))

      # test cell-wise MAE
      expect_lt(mean(abs(cb - cc)), .01)
})
