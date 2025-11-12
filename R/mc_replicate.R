# Internal helper: replicate fun() n_reps times, optionally in parallel.
# fun() should return something that can be simplified via simplify2array()
# (e.g. a vector, matrix, or array).
mc_replicate <- function(n_reps,
                         fun,
                         n_cores = 1L) {

      n_reps  <- as.integer(n_reps)
      n_cores <- as.integer(n_cores)

      if (n_reps <= 0L) {
            stop("n_reps must be positive")
      }

      # Sequential path: just use base::replicate to match semantics
      if (n_cores <= 1L) {
            return(replicate(n_reps, fun(), simplify = "array"))
      }

      # Simple, reproducible seeding per replicate
      seeds <- sample.int(.Machine$integer.max, n_reps)

      # Worker: i is the replicate index
      worker_fun <- function(i, fun, seeds) {
            set.seed(seeds[i])
            fun()
      }

      # Choose backend by OS
      if (.Platform$OS.type == "unix") {
            # Forked processes (not on Windows)
            sims_list <- parallel::mclapply(
                  X        = seq_len(n_reps),
                  FUN      = worker_fun,
                  fun      = fun,
                  seeds    = seeds,
                  mc.cores = n_cores
            )
      } else {
            # PSOCK clusters (works everywhere, including Windows)
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            sims_list <- parallel::parLapply(
                  cl,
                  X     = seq_len(n_reps),
                  fun   = worker_fun,
                  fun   = fun,
                  seeds = seeds
            )
      }

      # Simplify list -> array in the same way replicate(..., simplify="array") would
      simplify2array(sims_list)
}
