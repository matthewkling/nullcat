test_that("nullcat commsim functions run awithout error", {

      require_vegan()

      x <- cat_mat(100)

      # nullcat_commsim
      expect_no_error({
            cs <- nullcat_commsim(n_iter = 100, method = "curvecat")
            nm <- vegan::nullmodel(x, cs)
            sims <- simulate(nm, nsim = 10)
      })

      # nullcat_commsim_seq
      expect_no_error({
            cs <- nullcat_commsim_seq(method = "curvecat")
            nm <- vegan::nullmodel(x, cs)
            sims <- simulate(nm, nsim = 10, thin = 5, burnin = 40)
      })

})
