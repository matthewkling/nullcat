#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>

using namespace Rcpp;

struct CurvecatGroup {
      std::vector<int> cols;  // column indices in this (a,b) multiset
      int a;                  // smaller label
      int b;                  // larger label
      int n_ab;               // count of (a,b) orientation in ORIGINAL pair

      CurvecatGroup() : a(0), b(0), n_ab(0) {}
};

// Helper: perform one categorical curveball trade between two rows
void curvecat_pair(IntegerMatrix &mat, int r1, int r2) {
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
            if (a > b) std::swap(a, b);

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

            // ORIGINAL orientation: (a,b) vs (b,a)
            if (v1 == a && v2 == b) {
                  ++g.n_ab;
            }
            // otherwise it must be (b,a), so we don't increment
      }

      if (groups.empty()) return;

      // For each multiset group, shuffle and reassign
      for (auto &kv : groups) {
            CurvecatGroup &g = kv.second;
            std::vector<int> &cols = g.cols;
            const int m = (int)cols.size();
            if (m == 0) continue;

            const int a = g.a;
            const int b = g.b;
            const int n_ab = g.n_ab;

            // Shuffle columns (Fisher–Yates) using R RNG
            if (m > 1) {
                  for (int i = m - 1; i > 0; --i) {
                        int j = (int)std::floor(R::runif(0.0, 1.0) * (i + 1));
                        if (j < 0) j = 0;
                        if (j > i) j = i;
                        std::swap(cols[i], cols[j]);
                  }
            }

            // Assign first n_ab to (a,b), rest to (b,a)
            for (int i = 0; i < n_ab; ++i) {
                  int col = cols[i];
                  mat(r1, col) = a;
                  mat(r2, col) = b;
            }
            for (int i = n_ab; i < m; ++i) {
                  int col = cols[i];
                  mat(r1, col) = b;
                  mat(r2, col) = a;
            }
      }
}

// [[Rcpp::export]]
IntegerMatrix curvecat_cpp(IntegerMatrix mat, int n_iter) {
      const int nrow = mat.nrow();
      if (nrow < 2 || n_iter <= 0) return mat;

      RNGScope scope;  // hook into R's RNG

      for (int it = 0; it < n_iter; ++it) {
            // Sample two distinct rows uniformly
            int r1 = (int)std::floor(R::runif(0.0, (double)nrow));
            int r2 = (int)std::floor(R::runif(0.0, (double)(nrow - 1)));
            if (r2 >= r1) r2++;  // ensure r2 != r1

            curvecat_pair(mat, r1, r2);
      }

      return mat;
}
