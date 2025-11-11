
test_that("nullcat returns a result for all methods", {

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



test_that("`category` and `index` modes agree (given a shared seed)", {

      set.seed(123)
      m <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)

      for(method in nullcat_methods()){

            set.seed(123)
            m_rand <- nullcat(m, method = method, n_iter = 1000, output = "category")

            set.seed(123)
            idx <- nullcat(m, method = method, n_iter = 1000, output = "index")
            m_vec <- as.integer(m)
            m_rand2 <- matrix(m_vec[idx], nrow = nrow(m), ncol = ncol(m))

            expect_equal(m_rand, m_rand2)
      }


})



test_that("row and column multisets are fixed as indended", {

      compare_margin_sets <- function(a, b, margin){
            expect_equal(as.vector(apply(a, margin, sort)),
                         as.vector(apply(a, margin, sort)))
      }

      set.seed(123)
      m0 <- matrix(sample(1:4, 20*30, replace = T), nrow = 30)

      # methods with fixed rows and cols
      for(method in c("curvecat", "swapcat", "tswapcat")){
            m1 <- nullcat(m0, method = method, 1000)
            compare_margin_sets(m0, m1, 1)
            compare_margin_sets(m0, m1, 2)
      }

      # methods with fixed rows only
      for(method in c("r0cat")){
            m1 <- nullcat(m0, method = method, 1000)
            compare_margin_sets(m0, m1, 1)
      }

      # methods with fixed cols only
      for(method in c("c0cat")){
            m1 <- nullcat(m0, method = method, 1000)
            compare_margin_sets(m0, m1, 2)
      }

})


test_that("stationary distributions match vegan versions, given binary data", {
      skip_on_cran()

      # a small binary matrix
      set.seed(1234)
      b <- matrix(sample(0:1, 16, replace = T), nrow = 4)


      bins <- c("curveball", "swap", "tswap")
      cats <- c("curvecat", "swapcat", "tswapcat")

      for(i in seq_len(length(bins))){

            sim <- function(nullmod){
                  nsim <- 1e5
                  apply(as.array(simulate(nullmod, burnin = nsim/2, nsim = nsim)), c(1, 2), mean)
            }

            cb <- sim(vegan::nullmodel(b, method = bins[i]))
            cc <- sim(vegan::nullmodel(b, method = commsim_cat_seq(method = cats[i])))

            # test proportional MAE
            expect_lt(mean(abs(cb - cc) / cc), .1)
      }


})


