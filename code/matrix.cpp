#include "matrix.h"

Mat Mat::inv() const {
  int SIZE = m.size();
  Mat ret(SIZE);
  Mat original(SIZE);
  Mat tmp(SIZE);
  for (int i = 0; i < SIZE; i++) {
    ret[i][i] = 1;
    original[i] = m[i];
  }
  for (int i = 0; i < SIZE; i++) {
    if (original[i][i] == 0) {
      int k = i;
      for (int j = i + 1; j < SIZE; j++)
        if (original[j][i] == 1) {
          k = j;
          break;
        }
      std::swap(original[i], original[k]);
      std::swap(ret[i], ret[k]);
    }
    for (int j = i + 1; j < SIZE; j++) {
      if (original[j][i] == 1) {
        original[j] ^= original[i];
        ret[j] ^= ret[i];
      }
    }
  }
  for (int i = SIZE - 1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      if (original[j][i] == 1) {
        original[j] ^= original[i];
        ret[j] ^= ret[i];
      }
    }
  }
  return ret;
}

Mat Mat::trans() const {
  int SIZE = m.size();
  Mat ret(SIZE);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      ret[i][j] = this->m[j][i];
    }
  }
  return ret;
}

Mat Mat::row_i2j(int i, int j) const {
  Mat ret = *this;
  ret.m[j] ^= ret.m[i];
  return ret;
}

Mat Mat::col_i2j(int i, int j) const {
  int SIZE = m.size();
  Mat ret = *this;
  for (int k = 0; k < SIZE; k++) {
    if (ret.m[k][i]) ret.m[k].flip(j);
  }
  return ret;
}

bool Mat::can_depthone() const {
  int SIZE = m.size();
  for (int i = 0; i < SIZE; i++)
    if (m[i].count() > 2) return false;
  std::vector<int> cnt_col(SIZE, 0);
  for (int i = 0; i < SIZE; i++) {
    int cnt = 0;
    for (int j = 0; j < SIZE; j++)
      if (m[j][i]) cnt++;
    if (cnt > 2) return false;
    cnt_col[i] = cnt;
  }

  std::vector<int> cnt_r(SIZE, 0);
  std::vector<int> cnt_c(SIZE, 0);
  for (int i = 0; i < SIZE; i++)
    for (int j = 0; j < SIZE; j++) {
      if (m[i][j] && m[i].count() == 2 && cnt_col[j] == 2)
        cnt_r[i]++, cnt_c[j]++;
      if (cnt_r[i] > 1 || cnt_c[j] > 1) return false;
    }
  return true;
}