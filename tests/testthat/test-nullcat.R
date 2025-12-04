
test_that("nullcat returns a result for all methods", {

      set.seed(123)
      m1 <- cat_mat(500)

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
      m <- cat_mat(500)

      for(method in nullcat_methods()){

            seed <- 123
            m_rand <- nullcat(m, method = method, n_iter = 1000,
                              output = "category", swaps = "alternating", seed = seed)

            idx <- nullcat(m, method = method, n_iter = 1000,
                           output = "index", swaps = "alternating", seed = seed)
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
      m0 <- cat_mat(500)

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

      require_vegan()

      # a small binary matrix
      set.seed(1234)
      b <- matrix(sample(0:1, 16, replace = T), nrow = 4)

      bins <- c("curveball", "swap", "tswap")
      cats <- c("curvecat", "swapcat", "tswapcat")

      for(i in seq_len(length(bins))){

            sim <- function(nullmod){
                  nsim <- 1e4
                  apply(as.array(simulate(nullmod, burnin = nsim/2, nsim = nsim)), c(1, 2), mean)
            }

            cb <- sim(vegan::nullmodel(b, method = bins[i]))
            cc <- sim(vegan::nullmodel(b, method = nullcat_commsim_seq(method = cats[i])))

            # test proportional MAE
            expect_lt(mean(abs(cb - cc) / cc), .1)
      }
})


test_that("swap direction operates as expected for curvecat, swapcat, and tswapcat", {

      set.seed(123)
      x <- matrix(sample(1:4, 2500, replace = TRUE), 50)
      idx <- matrix(1:length(x), nrow(x))

      for(fun in list(curvecat_cpp, swapcat_cpp, tswapcat_cpp)){

            # vertical swaps: preserve column token sums, not row sums
            # (with token sums as a proxy for token sets)
            ver <- fun(x, 1000, swaps = "vertical", output = "index")
            expect_equal(colSums(idx), colSums(ver))
            expect_false(identical(rowSums(idx), rowSums(ver)))

            # horizontal swaps: preserve row token sums, not column sums
            hor <- fun(x, 1000, swaps = "horizontal", output = "index")
            expect_equal(rowSums(idx), rowSums(hor))
            expect_false(identical(colSums(idx), colSums(hor)))

            # alternating swaps: preserve neither margin
            alt <- fun(x, 1000, swaps = "alternating", output = "index")
            expect_false(identical(colSums(idx), colSums(alt)))
            expect_false(identical(rowSums(idx), rowSums(alt)))
      }

})
