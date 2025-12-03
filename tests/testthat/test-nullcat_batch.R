test_that("`nullcat_batch` works", {
      set.seed(123)
      m <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

      # single core
      expect_no_error(nullcat_batch(m, n_reps = 10, n_iter = 50, method = "curvecat"))

      # # multicore (skipping because devtools::check() doesn't like multicore)
      # expect_no_error(nullcat_batch(m, n_reps = 10, n_iter = 50, n_cores = 4, method = "curvecat"))
})
