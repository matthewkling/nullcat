#ifndef NULLCAT_WEIGHTED_PAIR_HPP
#define NULLCAT_WEIGHTED_PAIR_HPP

#include <vector>
#include <cmath>
#include <Rcpp.h>

// Weighted pair sampler for nullcat.
//
// Given an n x n weight matrix (symmetric, non-negative, diagonal ignored),
// builds an alias table over all n*(n-1)/2 unique pairs for O(1) sampling.
// When weights are NULL/empty, falls back to uniform pair sampling.

namespace nullcat {

struct WeightedPairSampler {
      int n;                    // dimension (number of rows or columns)
      int n_pairs;              // n*(n-1)/2
      bool weighted;            // false = uniform sampling (no alias table)

      // Alias table arrays (only populated when weighted = true)
      std::vector<double> prob;   // probability table
      std::vector<int>    alias;  // alias table

      // Precomputed pair lookup (avoids sqrt decode at sample time)
      std::vector<int> pair_i;    // first element of pair k
      std::vector<int> pair_j;    // second element of pair k

      // Sample a pair. Returns (a, b) with a != b (not necessarily a < b).
      inline void sample(int &a, int &b) const {
            if (!weighted) {
                  // Uniform sampling of two distinct elements from 0..n-1
                  a = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n))));
                  b = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n - 1))));
                  if (b >= a) b++;
                  return;
            }

            // Alias method: O(1) sampling from the pair distribution
            double u = R::runif(0.0, static_cast<double>(n_pairs));
            int k = static_cast<int>(std::floor(u));
            if (k >= n_pairs) k = n_pairs - 1; // guard against edge case
            double frac = u - static_cast<double>(k);

            int idx = (frac < prob[k]) ? k : alias[k];

            // Direct lookup — no sqrt needed
            a = pair_i[idx];
            b = pair_j[idx];

            // Randomly flip order so we don't always have a < b
            if (R::runif(0.0, 1.0) < 0.5) {
                  std::swap(a, b);
            }
      }
};

// Build a WeightedPairSampler from a flat upper-triangle weight vector.
// If weights is empty, returns an unweighted (uniform) sampler.
inline WeightedPairSampler make_pair_sampler(int n,
                                             const std::vector<double> &weights) {
      WeightedPairSampler S;
      S.n = n;
      S.n_pairs = n * (n - 1) / 2;
      S.weighted = !weights.empty();

      if (!S.weighted || S.n_pairs == 0) {
            return S;
      }

      // Build alias table using Vose's algorithm
      // Reference: Vose, M.D. (1991). A linear algorithm for generating random
      // numbers with a given distribution. IEEE Trans. Software Eng. 17(9).

      S.prob.resize(S.n_pairs);
      S.alias.resize(S.n_pairs);

      // Precompute pair lookup table
      S.pair_i.resize(S.n_pairs);
      S.pair_j.resize(S.n_pairs);
      {
            int k = 0;
            for (int i = 0; i < n; ++i) {
                  for (int j = i + 1; j < n; ++j) {
                        S.pair_i[k] = i;
                        S.pair_j[k] = j;
                        ++k;
                  }
            }
      }

      // Normalize weights
      double total = 0.0;
      for (int k = 0; k < S.n_pairs; ++k) {
            total += weights[k];
      }

      if (total <= 0.0) {
            // All weights zero: fall back to uniform
            S.weighted = false;
            S.prob.clear();
            S.alias.clear();
            return S;
      }

      // Scaled probabilities: p[k] * n_pairs
      std::vector<double> scaled(S.n_pairs);
      for (int k = 0; k < S.n_pairs; ++k) {
            scaled[k] = (weights[k] / total) * static_cast<double>(S.n_pairs);
      }

      // Partition into small and large
      std::vector<int> small, large;
      small.reserve(S.n_pairs);
      large.reserve(S.n_pairs);

      for (int k = 0; k < S.n_pairs; ++k) {
            if (scaled[k] < 1.0) {
                  small.push_back(k);
            } else {
                  large.push_back(k);
            }
      }

      while (!small.empty() && !large.empty()) {
            int s = small.back(); small.pop_back();
            int l = large.back(); large.pop_back();

            S.prob[s] = scaled[s];
            S.alias[s] = l;

            scaled[l] = (scaled[l] + scaled[s]) - 1.0;

            if (scaled[l] < 1.0) {
                  small.push_back(l);
            } else {
                  large.push_back(l);
            }
      }

      // Remaining entries get probability 1.0
      while (!large.empty()) {
            int l = large.back(); large.pop_back();
            S.prob[l] = 1.0;
            S.alias[l] = l; // self-alias (never used)
      }
      while (!small.empty()) {
            // Can happen due to floating point
            int s = small.back(); small.pop_back();
            S.prob[s] = 1.0;
            S.alias[s] = s;
      }

      return S;
}

// Build a WeightedPairSampler from an R NumericMatrix (full square matrix).
// Extracts upper triangle weights[i,j] for i < j.
// If wt_mat has 0 rows (i.e. is R_NilValue / empty), returns uniform sampler.
inline WeightedPairSampler make_pair_sampler_from_matrix(int n,
                                                         const Rcpp::NumericMatrix &wt_mat) {
      if (wt_mat.nrow() == 0) {
            // No weights supplied: uniform sampler
            std::vector<double> empty;
            return make_pair_sampler(n, empty);
      }

      int n_pairs = n * (n - 1) / 2;
      std::vector<double> weights(n_pairs);

      for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                  int k = i * n - i * (i + 1) / 2 + (j - i - 1);
                  double w = wt_mat(i, j);
                  if (w < 0.0) {
                        Rcpp::stop("Weight matrix contains negative values.");
                  }
                  weights[k] = w;
            }
      }

      return make_pair_sampler(n, weights);
}

} // namespace nullcat

#endif // NULLCAT_WEIGHTED_PAIR_HPP
