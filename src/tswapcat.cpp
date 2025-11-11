#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <unordered_map>  // needed by inline helpers in cat.hpp

#include "cat.hpp"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// tswapcat algorithm internals
//------------------------------------------------------------------------------

// Attempt one categorical 2x2 trial swap.
// Returns true if a swap was performed, false otherwise.
//
// We look for a 2x2 block:
//
//   a  b
//   c  d
//
// where a != b, a == d, b == c.
// This is the categorical analog of the binary 1-0 / 0-1 pattern,
// and we flip it to the opposite orientation (AB/BA <-> BA/AB).
static bool tswapcat_try_block(CatState &S, IntegerMatrix *idx,
                               int r1, int r2, int c1, int c2) {
      // Guard against degenerate cases
      if (r1 == r2 || c1 == c2) return false;

      int a = S.get(r1, c1);
      int b = S.get(r1, c2);
      int c = S.get(r2, c1);
      int d = S.get(r2, c2);

      // Check for AB/BA or BA/AB pattern in category ids
      if (!(a != b && a == d && b == c)) {
            return false; // no valid swap
      }

      // If tracking indices, pull out the 2x2 block of indices
      int ia, ib, ic, id;
      if (idx) {
            ia = (*idx)(r1, c1);
            ib = (*idx)(r1, c2);
            ic = (*idx)(r2, c1);
            id = (*idx)(r2, c2);
      }

      // Flip orientation:
      //   a  b        b  a
      //   c  d   ->   d  c
      S.set(r1, c1, b);
      S.set(r1, c2, a);
      S.set(r2, c1, d);
      S.set(r2, c2, c);

      if (idx) {
            (*idx)(r1, c1) = ib;
            (*idx)(r1, c2) = ia;
            (*idx)(r2, c1) = id;
            (*idx)(r2, c2) = ic;
      }

      return true;
}

// Run n_iter trial swaps on CatState, optionally tracking an index matrix.
//
// For each "iteration", we try up to max_trials randomly chosen 2x2 blocks.
// If any trial succeeds (i.e., we find a valid AB/BA block), we perform
// exactly one swap and move on to the next iteration.
static void tswapcat_engine(CatState &S, IntegerMatrix *idx, int n_iter) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      // Simple default: number of trials per iteration proportional to matrix size.
      const int max_trials = nrow * ncol;

      for (int it = 0; it < n_iter; ++it) {
            bool swapped = false;

            for (int trial = 0; trial < max_trials; ++trial) {
                  // Sample two distinct rows
                  int r1 = static_cast<int>(std::floor(R::runif(0.0, (double)nrow)));
                  int r2 = static_cast<int>(std::floor(R::runif(0.0, (double)(nrow - 1))));
                  if (r2 >= r1) r2++;

                  // Sample two distinct columns
                  int c1 = static_cast<int>(std::floor(R::runif(0.0, (double)ncol)));
                  int c2 = static_cast<int>(std::floor(R::runif(0.0, (double)(ncol - 1))));
                  if (c2 >= c1) c2++;

                  if (tswapcat_try_block(S, idx, r1, r2, c1, c2)) {
                        swapped = true;
                        break; // proceed to next "iteration"
                  }
            }

            // If no swap happened in max_trials attempts, this iteration
            // is effectively a no-op and we move on to the next one.
            (void)swapped; // silence unused-variable warning if not inspected
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix tswapcat_cpp(IntegerMatrix mat,
                           int n_iter,
                           std::string output = "category") {
      const int nrow = mat.nrow();
      const int ncol = mat.ncol();

      if (nrow < 2 || ncol < 2 || n_iter <= 0) {
            // Nothing to do; return input matrix unchanged
            return mat;
      }

      const bool return_index = (output == "index");

      RNGScope scope;  // hook into R's RNG

      // Build shared categorical state
      CatState S = make_cat_state_from_matrix(mat);

      // Optional index matrix for token tracking
      IntegerMatrix idx;
      IntegerMatrix *idx_ptr = nullptr;
      if (return_index) {
            idx = IntegerMatrix(nrow, ncol);
            int counter = 1;
            // COLUMN-major order to match R's as.vector()
            for (int j = 0; j < ncol; ++j) {
                  for (int i = 0; i < nrow; ++i) {
                        idx(i, j) = counter++;
                  }
            }
            idx_ptr = &idx;
      }

      // Run the engine
      tswapcat_engine(S, idx_ptr, n_iter);

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
