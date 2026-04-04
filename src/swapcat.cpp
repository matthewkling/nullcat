#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "cat.h"
#include "weighted_pair.h"

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
// For horizontal swaps, tokens move between columns, so the relevant
// weighted margin (if any) is the column margin.
static void swapcat_engine_horizontal(CatState &S, IntegerMatrix *idx, int n_iter,
                                      const WeightedPairSampler &col_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      // Row pairs are always sampled uniformly in swap (both rows and cols
      // are needed to identify the 2x2 block). The wt_row weights apply
      // to the margin corresponding to the swap direction.
      WeightedPairSampler row_sampler = make_pair_sampler(nrow, std::vector<double>());

      for (int it = 0; it < n_iter; ++it) {
            int r1, r2, c1, c2;
            row_sampler.sample(r1, r2);
            col_sampler.sample(c1, c2);
            swapcat_block_horizontal(S, idx, r1, r2, c1, c2);
      }
}

// Run n_iter swap attempts with vertical flips.
// For vertical swaps, tokens move between rows, so the relevant
// weighted margin (if any) is the row margin.
static void swapcat_engine_vertical(CatState &S, IntegerMatrix *idx, int n_iter,
                                    const WeightedPairSampler &row_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      // Column pairs always sampled uniformly
      WeightedPairSampler col_sampler = make_pair_sampler(ncol, std::vector<double>());

      for (int it = 0; it < n_iter; ++it) {
            int r1, r2, c1, c2;
            row_sampler.sample(r1, r2);
            col_sampler.sample(c1, c2);
            swapcat_block_vertical(S, idx, r1, r2, c1, c2);
      }
}

// Run n_iter swap attempts alternating between horizontal and vertical flips.
static void swapcat_engine_alternating(CatState &S, IntegerMatrix *idx, int n_iter,
                                       const WeightedPairSampler &row_sampler,
                                       const WeightedPairSampler &col_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      for (int it = 0; it < n_iter; ++it) {
            int r1, r2, c1, c2;
            row_sampler.sample(r1, r2);
            col_sampler.sample(c1, c2);

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
                          std::string output = "category",
                          Rcpp::Nullable<Rcpp::NumericMatrix> wt_row = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericMatrix> wt_col = R_NilValue) {
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
            WeightedPairSampler col_sampler = wt_col.isNotNull()
            ? make_pair_sampler_from_matrix(ncol, Rcpp::as<Rcpp::NumericMatrix>(wt_col))
                  : make_pair_sampler(ncol, std::vector<double>());
            swapcat_engine_horizontal(S, idx_ptr, n_iter, col_sampler);
      } else if (swaps == "vertical") {
            WeightedPairSampler row_sampler = wt_row.isNotNull()
            ? make_pair_sampler_from_matrix(nrow, Rcpp::as<Rcpp::NumericMatrix>(wt_row))
                  : make_pair_sampler(nrow, std::vector<double>());
            swapcat_engine_vertical(S, idx_ptr, n_iter, row_sampler);
      } else if (swaps == "alternating") {
            WeightedPairSampler row_sampler = wt_row.isNotNull()
            ? make_pair_sampler_from_matrix(nrow, Rcpp::as<Rcpp::NumericMatrix>(wt_row))
                  : make_pair_sampler(nrow, std::vector<double>());
            WeightedPairSampler col_sampler = wt_col.isNotNull()
                  ? make_pair_sampler_from_matrix(ncol, Rcpp::as<Rcpp::NumericMatrix>(wt_col))
                        : make_pair_sampler(ncol, std::vector<double>());

            swapcat_engine_alternating(S, idx_ptr, n_iter, row_sampler, col_sampler);
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
