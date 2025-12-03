#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>

#include "cat.h"

using namespace Rcpp;
using namespace nullcat;

//------------------------------------------------------------------------------
// Curvecat algorithm internals
//------------------------------------------------------------------------------

struct CurvecatGroup {
      std::vector<int>  cols;    // column indices in this (a,b) multiset
      std::vector<char> orient;  // 0 = (a,b), 1 = (b,a) in the ORIGINAL pair
      int a;                     // smaller label (category id)
      int b;                     // larger label (category id)
      int n_ab;                  // count of (a,b) orientation in ORIGINAL pair

      CurvecatGroup() : a(0), b(0), n_ab(0) {}
};

// Perform one categorical curveball trade between two rows of CatState.
// If idx is non-null, we also update an index matrix (R-style) in lockstep.
static void curvecat_pair(CatState &S, IntegerMatrix *idx, int r1, int r2) {
      const int ncol = S.n_col;
      if (ncol == 0 || r1 == r2) return;

      // Group differing columns by unordered pair (a,b),
      // while simultaneously counting (a,b) orientations.
      std::unordered_map<long long, CurvecatGroup> groups;
      groups.reserve(std::min(ncol, 32));

      for (int j = 0; j < ncol; ++j) {
            int v1 = S.get(r1, j);
            int v2 = S.get(r2, j);
            if (v1 == v2) continue;  // identical => no trade

            int a = v1;
            int b = v2;
            bool is_ab = true; // orientation in ORIGINAL pair

            if (a > b) {
                  std::swap(a, b);
                  is_ab = false;  // originally (b,a)
            }

            // encode unordered pair into 64-bit key
            long long key = (((long long)a) << 32) ^ (unsigned int)b;

            auto it = groups.find(key);
            if (it == groups.end()) {
                  CurvecatGroup g;
                  g.a = a;
                  g.b = b;
                  auto res = groups.emplace(key, std::move(g));
                  it = res.first;
            }

            CurvecatGroup &g = it->second;
            g.cols.push_back(j);
            g.orient.push_back(is_ab ? 0 : 1);

            // ORIGINAL orientation: (a,b) vs (b,a)
            if (is_ab) {
                  ++g.n_ab;
            }
      }

      if (groups.empty()) return;

      // For each multiset group, shuffle and reassign
      for (auto &kv : groups) {
            CurvecatGroup &g = kv.second;
            std::vector<int>  &cols   = g.cols;
            std::vector<char> &orient = g.orient;
            const int m = static_cast<int>(cols.size());
            if (m == 0) continue;

            const int a    = g.a;
            const int b    = g.b;
            const int n_ab = g.n_ab;

            // Shuffle columns using the shared helper, ensuring orient
            // gets the same permutation.
            shuffle_two(cols, orient);

            // Assign first n_ab to (a,b), rest to (b,a).
            // While doing so, if we have an index matrix, we swap indices
            // for columns whose orientation flips relative to the original.
            for (int i = 0; i < m; ++i) {
                  int col        = cols[i];
                  bool new_is_ab = (i < n_ab);        // TRUE if we assign (a,b)
                  bool old_is_ab = (orient[i] == 0);  // TRUE if originally (a,b)

                  // If we are tracking indices and the orientation flips,
                  // swap the indices in this column.
                  if (idx && (new_is_ab != old_is_ab)) {
                        int tmp = (*idx)(r1, col);
                        (*idx)(r1, col) = (*idx)(r2, col);
                        (*idx)(r2, col) = tmp;
                  }

                  // Now assign categories in S
                  if (new_is_ab) {
                        S.set(r1, col, a);
                        S.set(r2, col, b);
                  } else {
                        S.set(r1, col, b);
                        S.set(r2, col, a);
                  }
            }
      }
}

// Run n_iter curveball trades on CatState, optionally tracking an index matrix.
static void curvecat_engine(CatState &S, IntegerMatrix *idx, int n_iter) {
      const int nrow = S.n_row;
      if (nrow < 2 || n_iter <= 0) return;

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows uniformly
            int r1 = static_cast<int>(std::floor(R::runif(0.0, (double)nrow)));
            int r2 = static_cast<int>(std::floor(R::runif(0.0, (double)(nrow - 1))));
            if (r2 >= r1) r2++;  // ensure r2 != r1

            curvecat_pair(S, idx, r1, r2);
      }
}

//------------------------------------------------------------------------------
// Rcpp export
//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix curvecat_cpp(IntegerMatrix mat,
                           int n_iter,
                           std::string output = "category") {
      const int nrow = mat.nrow();
      const int ncol = mat.ncol();

      if (nrow < 2 || n_iter <= 0) {
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
      curvecat_engine(S, idx_ptr, n_iter);

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return cat_state_to_matrix(S);
      }
}
