#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "cat.h"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// r0cat: row-wise categorical permutation (row margin fixed, columns free)
//------------------------------------------------------------------------------

// Randomly permute the categories within each row, preserving the multiset
// of categories in that row. If idx is non-null, apply the same permutation
// to the index matrix so tokens move with their cells.
static void r0cat_engine(CatState &S, IntegerMatrix *idx) {
      const int nrow = S.n_row;
      const int ncol = S.n_col;
      if (nrow == 0 || ncol == 0) return;

      std::vector<int> order;
      order.reserve(ncol);
      std::vector<int> row_vals;
      row_vals.reserve(ncol);

      for (int i = 0; i < nrow; ++i) {
            // Build a permutation of column indices 0..ncol-1
            order.clear();
            for (int j = 0; j < ncol; ++j) {
                  order.push_back(j);
            }
            shuffle_in_place(order); // shared helper from cat.hpp

            // Copy original row categories
            row_vals.clear();
            for (int j = 0; j < ncol; ++j) {
                  row_vals.push_back(S.get(i, j));
            }

            // Apply permutation: new column j gets category from old column order[j]
            for (int j = 0; j < ncol; ++j) {
                  int old_j = order[j];
                  S.set(i, j, row_vals[old_j]);
            }

            // If tracking indices, apply the same permutation to idx
            if (idx) {
                  std::vector<int> row_idx;
                  row_idx.reserve(ncol);
                  for (int j = 0; j < ncol; ++j) {
                        row_idx.push_back((*idx)(i, j));
                  }
                  for (int j = 0; j < ncol; ++j) {
                        int old_j = order[j];
                        (*idx)(i, j) = row_idx[old_j];
                  }
            }
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix r0cat_cpp(IntegerMatrix mat,
                        int n_iter,
                        std::string output = "category") {
      const int nrow = mat.nrow();
      const int ncol = mat.ncol();

      // For r0cat, n_iter is effectively ignored beyond checking > 0:
      // any n_iter >= 1 gives a single random row-wise permutation.
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
      r0cat_engine(S, idx_ptr);

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
