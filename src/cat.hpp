#ifndef NULLCAT_CAT_HPP
#define NULLCAT_CAT_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <Rcpp.h>

// Shared categorical matrix state and utilities for nullcat
namespace nullcat {

struct CatState {
      int n_row = 0;
      int n_col = 0;
      int n_cat = 0;

      // Row-major storage: cells[i * n_col + j]
      std::vector<int> cells;
      // Mapping from category id -> original integer value
      std::vector<int> cat_values;

      inline int index(int i, int j) const {
            return i * n_col + j;
      }

      inline int get(int i, int j) const {
            return cells[index(i, j)];
      }

      inline void set(int i, int j, int value) {
            cells[index(i, j)] = value;
      }
};

// Fisher–Yates shuffle on a single container with operator[]
template <typename Vec>
inline void shuffle_in_place(Vec &x) {
      int m = static_cast<int>(x.size());
      if (m <= 1) return;

      for (int i = m - 1; i > 0; --i) {
            int j = static_cast<int>(std::floor(R::runif(0.0, 1.0) * (i + 1)));
            if (j < 0) j = 0;
            if (j > i) j = i;
            std::swap(x[i], x[j]);
      }
}

// Fisher–Yates shuffle applied to two parallel containers of equal length
template <typename Vec1, typename Vec2>
inline void shuffle_two(Vec1 &a, Vec2 &b) {
      int m = static_cast<int>(a.size());
      if (static_cast<int>(b.size()) != m || m <= 1) return;

      for (int i = m - 1; i > 0; --i) {
            int j = static_cast<int>(std::floor(R::runif(0.0, 1.0) * (i + 1)));
            if (j < 0) j = 0;
            if (j > i) j = i;
            std::swap(a[i], a[j]);
            std::swap(b[i], b[j]);
      }
}

// Build a CatState from an R integer matrix.
// - All distinct integer values (including NA) are mapped to ids 0..n_cat-1.
// - Cell storage is row-major for C++ convenience.
inline CatState make_cat_state_from_matrix(const Rcpp::IntegerMatrix &mat) {
      const int n_row = mat.nrow();
      const int n_col = mat.ncol();

      CatState S;
      S.n_row = n_row;
      S.n_col = n_col;
      S.n_cat = 0;

      std::unordered_map<int, int> value_to_id;
      value_to_id.reserve(std::max(16, n_row * n_col / 4));

      S.cells.resize(static_cast<size_t>(n_row) * static_cast<size_t>(n_col));

      for (int i = 0; i < n_row; ++i) {
            for (int j = 0; j < n_col; ++j) {
                  int v = mat(i, j);
                  auto it = value_to_id.find(v);
                  int id;
                  if (it == value_to_id.end()) {
                        id = S.n_cat;
                        value_to_id.emplace(v, id);
                        S.cat_values.push_back(v);
                        ++S.n_cat;
                  } else {
                        id = it->second;
                  }
                  S.cells[S.index(i, j)] = id;
            }
      }

      return S;
}

// Convert CatState back to an R IntegerMatrix, restoring original values
// from cat_values.
inline Rcpp::IntegerMatrix cat_state_to_matrix(const CatState &S) {
      Rcpp::IntegerMatrix mat(S.n_row, S.n_col);

      for (int i = 0; i < S.n_row; ++i) {
            for (int j = 0; j < S.n_col; ++j) {
                  int id = S.get(i, j);
                  mat(i, j) = S.cat_values[id];
            }
      }

      return mat;
}

} // namespace nullcat

#endif // NULLCAT_CAT_HPP
