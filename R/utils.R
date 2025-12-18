

# internal helper -- simulate a categorical matrix for use in testing.
cat_mat <- function(n = NULL, nr = NULL, nc = NULL, ncat = 4, prob = NULL){
      if(!is.null(n)){
            nr <- round(sqrt(n))
            nc <- round(sqrt(n))
      }
      matrix(sample(1:ncat, nr * nc, replace = T, prob = prob), nrow = nr)
}


# internal helper -- require vegan package.
require_vegan <- function() {
      if (!requireNamespace("vegan", quietly = TRUE)) {
            stop("Package 'vegan' is required for this operation. ",
                 "Please install it with install.packages('vegan').",
                 call. = FALSE)
      }
}


# internal helper -- set seed, run expression, reinstate prior seed.
with_seed <- function(seed, expr) {
      if (is.null(seed)) return(eval.parent(substitute(expr)))

      # Work with parent environment instead of .GlobalEnv
      parent_env <- parent.frame()
      old_exists <- exists(".Random.seed", envir = parent_env, inherits = FALSE)
      if (old_exists) old_seed <- get(".Random.seed", envir = parent_env, inherits = FALSE)

      set.seed(as.integer(seed))

      on.exit({
            if (old_exists) {
                  assign(".Random.seed", old_seed, envir = parent_env)
            } else if (exists(".Random.seed", envir = parent_env, inherits = FALSE)) {
                  rm(".Random.seed", envir = parent_env)
            }
      }, add = TRUE)

      eval.parent(substitute(expr))
}


# internal helper -- replicate fun() n_reps times, optionally in parallel.
# fun() should return something that can be simplified via simplify2array()
# (e.g. a vector, matrix, or array).
# Updated signature to accept additional arguments
mc_replicate <- function(n_reps, fun, n_cores = 1L) {

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

      # Choose backend by OS
      if (.Platform$OS.type == "unix") {
            # Forked processes (closures work automatically)
            sims_list <- parallel::mclapply(
                  seq_len(n_reps),
                  function(i) {
                        set.seed(seeds[i])
                        fun()
                  },
                  mc.cores = n_cores
            )
      } else {
            # PSOCK clusters - need to export function's environment
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            # Export everything from fun's closure
            parallel::clusterExport(cl,
                                    varlist = ls(environment(fun), all.names = TRUE),
                                    envir = environment(fun))

            sims_list <- parallel::parLapply(
                  cl,
                  seq_len(n_reps),
                  function(i) {
                        set.seed(seeds[i])
                        fun()
                  }
            )
      }

      # Simplify list -> array
      simplify2array(sims_list)
}
