#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <unordered_map>

#include "cat.h"
#include "weighted_pair.h"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// tswapcat algorithm internals
//------------------------------------------------------------------------------

// Attempt one categorical 2x2 trial swap with horizontal flip.
// Returns true if a swap was performed, false otherwise.
static bool tswapcat_try_block_horizontal(CatState &S, IntegerMatrix *idx,
                                          int r1, int r2, int c1, int c2) {
      if (r1 == r2 || c1 == c2) return false;

      int a = S.get(r1, c1);
      int b = S.get(r1, c2);
      int c = S.get(r2, c1);
      int d = S.get(r2, c2);

      if (!(a != b && a == d && b == c)) {
            return false;
      }

      int ia, ib, ic, id;
      if (idx) {
            ia = (*idx)(r1, c1);
            ib = (*idx)(r1, c2);
            ic = (*idx)(r2, c1);
            id = (*idx)(r2, c2);
      }

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

// Attempt one categorical 2x2 trial swap with vertical flip.
// Returns true if a swap was performed, false otherwise.
static bool tswapcat_try_block_vertical(CatState &S, IntegerMatrix *idx,
                                        int r1, int r2, int c1, int c2) {
      if (r1 == r2 || c1 == c2) return false;

      int a = S.get(r1, c1);
      int b = S.get(r1, c2);
      int c = S.get(r2, c1);
      int d = S.get(r2, c2);

      if (!(a != b && a == d && b == c)) {
            return false;
      }

      int ia, ib, ic, id;
      if (idx) {
            ia = (*idx)(r1, c1);
            ib = (*idx)(r1, c2);
            ic = (*idx)(r2, c1);
            id = (*idx)(r2, c2);
      }

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

      return true;
}

// Run n_iter trial swaps with horizontal flips.
static void tswapcat_engine_horizontal(CatState &S, IntegerMatrix *idx, int n_iter,
                                       const WeightedPairSampler &col_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      const int max_trials = nrow * ncol;
      WeightedPairSampler row_sampler = make_pair_sampler(nrow, std::vector<double>());

      for (int it = 0; it < n_iter; ++it) {
            bool swapped = false;

            for (int trial = 0; trial < max_trials; ++trial) {
                  int r1, r2, c1, c2;
                  row_sampler.sample(r1, r2);
                  col_sampler.sample(c1, c2);

                  if (tswapcat_try_block_horizontal(S, idx, r1, r2, c1, c2)) {
                        swapped = true;
                        break;
                  }
            }

            (void)swapped;
      }
}

// Run n_iter trial swaps with vertical flips.
static void tswapcat_engine_vertical(CatState &S, IntegerMatrix *idx, int n_iter,
                                     const WeightedPairSampler &row_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      const int max_trials = nrow * ncol;
      WeightedPairSampler col_sampler = make_pair_sampler(ncol, std::vector<double>());

      for (int it = 0; it < n_iter; ++it) {
            bool swapped = false;

            for (int trial = 0; trial < max_trials; ++trial) {
                  int r1, r2, c1, c2;
                  row_sampler.sample(r1, r2);
                  col_sampler.sample(c1, c2);

                  if (tswapcat_try_block_vertical(S, idx, r1, r2, c1, c2)) {
                        swapped = true;
                        break;
                  }
            }

            (void)swapped;
      }
}

// Run n_iter trial swaps alternating between horizontal and vertical flips.
static void tswapcat_engine_alternating(CatState &S, IntegerMatrix *idx, int n_iter,
                                        const WeightedPairSampler &row_sampler,
                                        const WeightedPairSampler &col_sampler) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;

      if (nrow < 2 || ncol < 2 || n_iter <= 0) return;

      const int max_trials = nrow * ncol;

      for (int it = 0; it < n_iter; ++it) {
            bool swapped = false;

            for (int trial = 0; trial < max_trials; ++trial) {
                  int r1, r2, c1, c2;
                  row_sampler.sample(r1, r2);
                  col_sampler.sample(c1, c2);

                  bool success;
                  if (it % 2 == 0) {
                        success = tswapcat_try_block_horizontal(S, idx, r1, r2, c1, c2);
                  } else {
                        success = tswapcat_try_block_vertical(S, idx, r1, r2, c1, c2);
                  }

                  if (success) {
                        swapped = true;
                        break;
                  }
            }

            (void)swapped;
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix tswapcat_cpp(IntegerMatrix mat,
                           int n_iter,
                           std::string swaps = "vertical",
                           std::string output = "category",
                           Rcpp::Nullable<Rcpp::NumericMatrix> wt_row = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericMatrix> wt_col = R_NilValue) {
      const int nrow = mat.nrow();
      const int ncol = mat.ncol();

      if (nrow < 2 || ncol < 2 || n_iter <= 0) {
            return mat;
      }

      const bool return_index = (output == "index");

      RNGScope scope;

      CatState S = make_cat_state_from_matrix(mat);

      IntegerMatrix idx;
      IntegerMatrix *idx_ptr = nullptr;
      if (return_index) {
            idx = IntegerMatrix(nrow, ncol);
            int counter = 1;
            for (int j = 0; j < ncol; ++j) {
                  for (int i = 0; i < nrow; ++i) {
                        idx(i, j) = counter++;
                  }
            }
            idx_ptr = &idx;
      }

      if (swaps == "horizontal") {
            WeightedPairSampler col_sampler = wt_col.isNotNull()
            ? make_pair_sampler_from_matrix(ncol, Rcpp::as<Rcpp::NumericMatrix>(wt_col))
                  : make_pair_sampler(ncol, std::vector<double>());
            tswapcat_engine_horizontal(S, idx_ptr, n_iter, col_sampler);
      } else if (swaps == "vertical") {
            WeightedPairSampler row_sampler = wt_row.isNotNull()
            ? make_pair_sampler_from_matrix(nrow, Rcpp::as<Rcpp::NumericMatrix>(wt_row))
                  : make_pair_sampler(nrow, std::vector<double>());
            tswapcat_engine_vertical(S, idx_ptr, n_iter, row_sampler);
      } else if (swaps == "alternating") {
            WeightedPairSampler row_sampler = wt_row.isNotNull()
            ? make_pair_sampler_from_matrix(nrow, Rcpp::as<Rcpp::NumericMatrix>(wt_row))
                  : make_pair_sampler(nrow, std::vector<double>());
            WeightedPairSampler col_sampler = wt_col.isNotNull()
                  ? make_pair_sampler_from_matrix(ncol, Rcpp::as<Rcpp::NumericMatrix>(wt_col))
                        : make_pair_sampler(ncol, std::vector<double>());

            tswapcat_engine_alternating(S, idx_ptr, n_iter, row_sampler, col_sampler);
      } else {
            Rcpp::stop("Invalid swaps parameter. Must be 'vertical', 'horizontal', or 'alternating'.");
      }

      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
