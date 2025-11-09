#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>

using namespace Rcpp;

struct CurvecatGroup {
      std::vector<int> cols;     // column indices in this (a,b) multiset
      std::vector<char> orient;  // 0 = (a,b), 1 = (b,a) in the ORIGINAL pair
      int a;                     // smaller label
      int b;                     // larger label
      int n_ab;                  // count of (a,b) orientation in ORIGINAL pair

      CurvecatGroup() : a(0), b(0), n_ab(0) {}
};

// Helper: perform one categorical curveball trade between two rows.
// If idx is non-null, we also update an index matrix in lockstep.
// - mat: category matrix
// - idx: optional index matrix (same dimensions as mat); may be nullptr
// - r1, r2: row indices to trade
void curvecat_pair(IntegerMatrix &mat, IntegerMatrix *idx, int r1, int r2) {
      const int ncol = mat.ncol();
      if (ncol == 0 || r1 == r2) return;

      // Group differing columns by unordered pair (a,b),
      // while simultaneously counting (a,b) orientations.
      std::unordered_map<long long, CurvecatGroup> groups;
      groups.reserve(std::min(ncol, 32));

      for (int j = 0; j < ncol; ++j) {
            int v1 = mat(r1, j);
            int v2 = mat(r2, j);
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
            std::vector<int> &cols = g.cols;
            std::vector<char> &orient = g.orient;
            const int m = (int)cols.size();
            if (m == 0) continue;

            const int a = g.a;
            const int b = g.b;
            const int n_ab = g.n_ab;

            // Shuffle columns (Fisher–Yates) using R RNG, and apply
            // the same permutation to the orientation vector.
            if (m > 1) {
                  for (int i = m - 1; i > 0; --i) {
                        int j = (int)std::floor(R::runif(0.0, 1.0) * (i + 1));
                        if (j < 0) j = 0;
                        if (j > i) j = i;
                        std::swap(cols[i], cols[j]);
                        std::swap(orient[i], orient[j]);
                  }
            }

            // Assign first n_ab to (a,b), rest to (b,a).
            // While doing so, if we have an index matrix, we swap indices
            // for columns whose orientation flips relative to the original.
            for (int i = 0; i < m; ++i) {
                  int col = cols[i];
                  bool new_is_ab = (i < n_ab);       // TRUE if we assign (a,b)
                  bool old_is_ab = (orient[i] == 0); // TRUE if originally (a,b)

                  // If we are tracking indices and the orientation flips,
                  // swap the indices in this column.
                  if (idx && (new_is_ab != old_is_ab)) {
                        int tmp = (*idx)(r1, col);
                        (*idx)(r1, col) = (*idx)(r2, col);
                        (*idx)(r2, col) = tmp;
                  }

                  // Now assign categories in mat
                  if (new_is_ab) {
                        mat(r1, col) = a;
                        mat(r2, col) = b;
                  } else {
                        mat(r1, col) = b;
                        mat(r2, col) = a;
                  }
            }
      }
}

// [[Rcpp::export]]
IntegerMatrix curvecat_cpp(IntegerMatrix mat,
                           int n_iter,
                           std::string output = "category") {

      const int nrow = mat.nrow();
      const int ncol = mat.ncol();
      if (nrow < 2 || n_iter <= 0) return mat;

      bool return_index = (output == "index");

      RNGScope scope;  // hook into R's RNG

      // Optional index matrix for token tracking
      IntegerMatrix idx;
      if (return_index) {
            idx = IntegerMatrix(nrow, ncol);
            int counter = 1;
            // COLUMN-major order to match R's as.vector()
            for (int j = 0; j < ncol; ++j) {
                  for (int i = 0; i < nrow; ++i) {
                        idx(i, j) = counter++;
                  }
            }
      }

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows uniformly
            int r1 = (int)std::floor(R::runif(0.0, (double)nrow));
            int r2 = (int)std::floor(R::runif(0.0, (double)(nrow - 1)));
            if (r2 >= r1) r2++;  // ensure r2 != r1

            if (return_index) {
                  curvecat_pair(mat, &idx, r1, r2);
            } else {
                  curvecat_pair(mat, nullptr, r1, r2);
            }
      }

      // Return either the randomized categories or the index mapping
      if (return_index) {
            return idx;
      } else {
            return mat;
      }
}
