#pragma once

#include <bits/stdc++.h>

template <typename T>
using Matrix = std::vector<std::vector<T>>;

class dynamic_bitset {
  static constexpr size_t SIZE = 64;
  size_t _size;

 public:
  std::bitset<SIZE> impl;

  dynamic_bitset() = default;
  explicit dynamic_bitset(size_t size) : _size(size), impl(0) {
    assert(size <= SIZE);
  }

  std::bitset<SIZE>::reference operator[](size_t n) { return impl[n]; }

  bool operator[](size_t n) const { return impl[n]; }

  void flip(size_t position) { impl.flip(position); }

  void reset() { impl.reset(); }

  size_t count() const { return impl.count(); }

  size_t size() const { return _size; }

  dynamic_bitset &operator^=(const dynamic_bitset &rhs) {
    assert(impl.size() == rhs.impl.size());
    impl ^= rhs.impl;
    return *this;
  }
};

class Mat {
 public:
  std::vector<dynamic_bitset> m;

  Mat() = default;
  explicit Mat(size_t size) { m.resize(size, dynamic_bitset(size)); }
  Mat(const Matrix<int> &mat, int size = -1)
      : Mat(size == -1 ? mat.size() : size) {
    if (size == -1) {
      size = mat.size();
    }
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        m[i][j] = mat[i][j];
      }
    }
  }

  dynamic_bitset &operator[](size_t n) { return m[n]; }

  const dynamic_bitset &operator[](size_t n) const { return m[n]; }

  size_t count() const {
    size_t cnt = 0;
    for (const auto &l : m) cnt += l.count();
    return cnt;
  }

  size_t size() const { return m.size(); }

  Mat inv() const;
  Mat trans() const;
  Mat row_i2j(int i, int j) const;
  Mat col_i2j(int i, int j) const;
  bool can_depthone() const;
};

template <>
struct std::formatter<Mat> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const Mat &mat, std::format_context &ctx) const {
    size_t SIZE = mat.size();
    for (unsigned i = 0; i < SIZE; ++i) {
      const auto &layer = mat[i];
      for (unsigned j = 0; j < SIZE; ++j) {
        std::format_to(ctx.out(), "{} ", (int)layer[j]);
      }
      std::format_to(ctx.out(), "\n");
    }
    return ctx.out();
  }
};

typedef struct {
  int src;
  int dst;
  int type;  // 0: row; 1: col.
} Op;