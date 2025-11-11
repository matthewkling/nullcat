#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "cat.hpp"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// c0cat: column-wise categorical permutation (column margin fixed, rows free)
//------------------------------------------------------------------------------

// Randomly permute the categories within each column, preserving the multiset
// of categories in that column. If idx is non-null, apply the same permutation
// to the index matrix so tokens move with their cells.
static void c0cat_engine(CatState &S, IntegerMatrix *idx) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;
      if (nrow == 0 || ncol == 0) return;

      std::vector<int> order;
      order.reserve(nrow);
      std::vector<int> col_vals;
      col_vals.reserve(nrow);

      for (int j = 0; j < ncol; ++j) {
            // Build a permutation of row indices 0..nrow-1
            order.clear();
            for (int i = 0; i < nrow; ++i) {
                  order.push_back(i);
            }
            shuffle_in_place(order); // shared helper from cat.hpp

            // Copy original column categories
            col_vals.clear();
            for (int i = 0; i < nrow; ++i) {
                  col_vals.push_back(S.get(i, j));
            }

            // Apply permutation: new row i gets category from old row order[i]
            for (int i = 0; i < nrow; ++i) {
                  int old_i = order[i];
                  S.set(i, j, col_vals[old_i]);
            }

            // If tracking indices, apply the same permutation to idx
            if (idx) {
                  std::vector<int> col_idx;
                  col_idx.reserve(nrow);
                  for (int i = 0; i < nrow; ++i) {
                        col_idx.push_back((*idx)(i, j));
                  }
                  for (int i = 0; i < nrow; ++i) {
                        int old_i = order[i];
                        (*idx)(i, j) = col_idx[old_i];
                  }
            }
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix c0cat_cpp(IntegerMatrix mat,
                        int n_iter,
                        std::string output = "category") {
      const int nrow = mat.nrow();
      const int ncol = mat.ncol();

      // For c0cat, n_iter is effectively ignored beyond checking > 0:
      // any n_iter >= 1 gives a single random column-wise permutation.
      if (nrow == 0 || ncol == 0 || n_iter <= 0) {
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

      // Run the engine (one randomization)
      c0cat_engine(S, idx_ptr);

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
