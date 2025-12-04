#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "cat.h"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// Swapcat algorithm internals
//------------------------------------------------------------------------------

// Perform one categorical swap between two rows and two columns with horizontal flip.
// Horizontal flip swaps columns within the 2x2 block (tokens move between columns).
//
// We look for a 2x2 block:
//   a  b
//   c  d
// where a != b, a == d, b == c, and flip it horizontally to:
//   b  a
//   d  c
static void swapcat_block_horizontal(CatState &S, IntegerMatrix *idx,
                                     int r1, int r2, int c1, int c2) {
      // Guard against degenerate cases
      if (r1 == r2 || c1 == c2) return;

      int a = S.get(r1, c1);
      int b = S.get(r1, c2);
      int c = S.get(r2, c1);
      int d = S.get(r2, c2);

      // Check for AB/BA or BA/AB pattern in category ids
      if (!(a != b && a == d && b == c)) {
            return; // no valid swap
      }

      // If tracking indices, pull out the 2x2 block of indices
      int ia, ib, ic, id;
      if (idx) {
            ia = (*idx)(r1, c1);
            ib = (*idx)(r1, c2);
            ic = (*idx)(r2, c1);
            id = (*idx)(r2, c2);
      }

      // Flip horizontally:
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
}

// Perform one categorical swap between two rows and two columns with vertical flip.
// Vertical flip swaps rows within the 2x2 block (tokens move between rows).
//
// We look for a 2x2 block:
//   a  b
//   c  d
// where a != b, a == d, b == c, and flip it vertically to:
//   c  d
//   a  b
static void swapcat_block_vertical(CatState &S, IntegerMatrix *idx,
                                   int r1, int r2, int c1, int c2) {
      // Guard against degenerate cases
      if (r1 == r2 || c1 == c2) return;

      int a = S.get(r1, c1);
      int b = S.get(r1, c2);
      int c = S.get(r2, c1);
      int d = S.get(r2, c2);

      // Check for AB/BA or BA/AB pattern in category ids
      if (!(a != b && a == d && b == c)) {
            return; // no valid swap
      }

      // If tracking indices, pull out the 2x2 block of indices
      int ia, ib, ic, id;
      if (idx) {
            ia = (*idx)(r1, c1);
            ib = (*idx)(r1, c2);
            ic = (*idx)(r2, c1);
            id = (*idx)(r2, c2);
      }

      // Flip vertically:
      //   a  b        c  d
      //   c  d   ->   a  b
      S.set(r1, c1, c);
      S.set(r1, c2, d);
      S.set(r2, c1, a);
      S.set(r2, c2, b);

      if (idx) {
            (*idx)(r1, c1) = ic;
            (*idx)(r1, c2) = id;
            (*idx)(r2, c1) = ia;
            (*idx)(r2, c2) = ib;
      }
}

// Run n_iter swap attempts with horizontal flips.
static void swapcat_engine_horizontal(CatState &S, IntegerMatrix *idx, int n_iter) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows
            int r1 = static_cast<int>(std::floor(R::runif(0.0, (double)nrow)));
            int r2 = static_cast<int>(std::floor(R::runif(0.0, (double)(nrow - 1))));
            if (r2 >= r1) r2++;

            // Sample two distinct columns
            int c1 = static_cast<int>(std::floor(R::runif(0.0, (double)ncol)));
            int c2 = static_cast<int>(std::floor(R::runif(0.0, (double)(ncol - 1))));
            if (c2 >= c1) c2++;

            swapcat_block_horizontal(S, idx, r1, r2, c1, c2);
      }
}

// Run n_iter swap attempts with vertical flips.
static void swapcat_engine_vertical(CatState &S, IntegerMatrix *idx, int n_iter) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows
            int r1 = static_cast<int>(std::floor(R::runif(0.0, (double)nrow)));
            int r2 = static_cast<int>(std::floor(R::runif(0.0, (double)(nrow - 1))));
            if (r2 >= r1) r2++;

            // Sample two distinct columns
            int c1 = static_cast<int>(std::floor(R::runif(0.0, (double)ncol)));
            int c2 = static_cast<int>(std::floor(R::runif(0.0, (double)(ncol - 1))));
            if (c2 >= c1) c2++;

            swapcat_block_vertical(S, idx, r1, r2, c1, c2);
      }
}

// Run n_iter swap attempts alternating between horizontal and vertical flips.
static void swapcat_engine_alternating(CatState &S, IntegerMatrix *idx, int n_iter) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows
            int r1 = static_cast<int>(std::floor(R::runif(0.0, (double)nrow)));
            int r2 = static_cast<int>(std::floor(R::runif(0.0, (double)(nrow - 1))));
            if (r2 >= r1) r2++;

            // Sample two distinct columns
            int c1 = static_cast<int>(std::floor(R::runif(0.0, (double)ncol)));
            int c2 = static_cast<int>(std::floor(R::runif(0.0, (double)(ncol - 1))));
            if (c2 >= c1) c2++;

            // Alternate between horizontal and vertical
            if (it % 2 == 0) {
                  swapcat_block_horizontal(S, idx, r1, r2, c1, c2);
            } else {
                  swapcat_block_vertical(S, idx, r1, r2, c1, c2);
            }
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix swapcat_cpp(IntegerMatrix mat,
                          int n_iter,
                          std::string swaps = "vertical",
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

      // Run the appropriate engine based on swaps parameter
      if (swaps == "horizontal") {
            swapcat_engine_horizontal(S, idx_ptr, n_iter);
      } else if (swaps == "vertical") {
            swapcat_engine_vertical(S, idx_ptr, n_iter);
      } else if (swaps == "alternating") {
            swapcat_engine_alternating(S, idx_ptr, n_iter);
      } else {
            Rcpp::stop("Invalid swaps parameter. Must be 'vertical', 'horizontal', or 'alternating'.");
      }

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
