test_that("quantize commsim functions run, and can be used, without erroring", {

      require_vegan()

      x <- matrix(runif(100), 10)

      # quantize_commsim
      expect_no_error({
            cs <- quantize_commsim(n_iter = 100, method = "curvecat", n_strata = 3)
            nm <- vegan::nullmodel(x, cs)
            sims <- simulate(nm, nsim = 10)
      })

      # quantize_commsim_seq
      expect_no_error({
            cs <- quantize_commsim_seq(method = "curvecat", transform = sqrt)
            nm <- vegan::nullmodel(x, cs)
            sims <- simulate(nm, nsim = 10, thin = 5, burnin = 40)
      })

})
