#include "greedy.h"

namespace AsiaCrypt24 {

double h_sq(const Mat &m) {
  int SIZE = m.size();
  double ret = 0;
  for (int i = 0; i < SIZE; i++) ret += m[i].count() * m[i].count();
  return ret;
}

double h_sum(const Mat &m) {
  int SIZE = m.size();
  double ret = 0;
  for (int i = 0; i < SIZE; i++) ret += m[i].count();
  return ret;
}

double h_prod(const Mat &m) {
  int SIZE = m.size();
  double ret = 0;
  for (int i = 0; i < SIZE; i++) ret += log2(m[i].count());
  return ret;
}

int cost_type = 2;
void rand_cost_type() { cost_type = rand() % 3; }

double cost(const Mat &m) {
  switch (cost_type) {
    case 0:
      return h_sum(m);
    case 1:
      return h_prod(m);
    case 2:
      return h_sq(m);
    default:
      return h_sq(m);
  }
}

constexpr double epsilon = 1e-3;
double cost_min(int size) {
  switch (cost_type) {
    case 0:
      return size + epsilon;
    case 1:
      return epsilon;
    case 2:
      return size + epsilon;
    default:
      return size + epsilon;
  }
}

std::optional<std::vector<std::vector<Op>>> greedy(Mat mat, int depth_limit,
                                                   int size_limit) {
  rand_cost_type();

  int SIZE = mat.size();
  Mat origin = mat;
  Mat inv = mat.inv();
  Mat tmp_mat;
  Mat tmp_inv;

  int depth = 0;
  int size = 0;
  double minm_cst = DBL_MAX;  // large
  double tmp_cst;

  std::vector<Op> select_list;
  std::vector<Op> r_layer, c_layer;
  std::vector<std::vector<Op>> r_layers, c_layers;
  dynamic_bitset row_visi(SIZE), col_visi(SIZE);
  std::vector<Op> r_op, c_op;
  bool one = false;
  while (cost(mat) > cost_min(SIZE)) {
    select_list.clear();
    minm_cst =
        std::max(cost(mat) + cost(inv.trans()), cost(mat.trans()) + cost(inv));

    if (!one) {
      for (int i = 0; i < SIZE; i++) {
        if (row_visi[i] == 1) continue;
        for (int j = 0; j < SIZE; j++) {
          if (row_visi[j] == 1 || j == i) continue;
          tmp_mat = mat.row_i2j(i, j);
          tmp_inv = inv.col_i2j(j, i);
          tmp_cst = cost(tmp_mat) + cost(tmp_inv.trans());
          if (tmp_cst < minm_cst + epsilon) {
            if (tmp_cst < minm_cst - epsilon) {
              select_list.clear();
              Op op = {i, j, 0};
              select_list.push_back(op);
              minm_cst = tmp_cst;
            } else {
              Op op = {i, j, 0};
              select_list.push_back(op);
            }
          }
        }
      }
    }

    // col trans
    for (int i = 0; i < SIZE; i++) {
      if (col_visi[i] == 1) continue;
      for (int j = 0; j < SIZE; j++) {
        if (col_visi[j] == 1 || j == i) continue;
        tmp_mat = mat.col_i2j(i, j);
        tmp_inv = inv.row_i2j(j, i);
        tmp_cst = cost(tmp_mat.trans()) + cost(tmp_inv);
        if (tmp_cst < minm_cst + epsilon) {
          if (tmp_cst < minm_cst - epsilon) {
            select_list.clear();
            Op op = {i, j, 1};
            select_list.push_back(op);
            minm_cst = tmp_cst;
          } else {
            Op op = {i, j, 1};
            select_list.push_back(op);
          }
        }
      }
    }

    if (select_list.size() == 0) {
      if (r_layer.size()) {
        r_layers.push_back(r_layer);
        r_layer.clear();
        row_visi.reset();
      }
      if (c_layer.size()) {
        c_layers.push_back(c_layer);
        c_layer.clear();
        col_visi.reset();
      }
      if (mat.can_depthone()) {
        one = true;
      }
    } else {
      int rd = rand() % select_list.size();
      auto [c, t, t01] = select_list[rd];
      if (t01 == 0) {
        mat = mat.row_i2j(c, t);
        inv = inv.col_i2j(t, c);
        r_layer.push_back({c, t, 0});
        r_op.push_back({c, t, 0});
        if (row_visi.count() == 0) depth++;
        row_visi[c] = 1;
        row_visi[t] = 1;
      } else {
        mat = mat.col_i2j(c, t);
        inv = inv.row_i2j(t, c);
        c_layer.push_back({c, t, 1});
        c_op.push_back({c, t, 1});
        if (col_visi.count() == 0) depth++;
        col_visi[c] = 1;
        col_visi[t] = 1;
      }
    }

    if (depth > depth_limit) {
      return std::nullopt;
    }
    if (++size > size_limit) {
      return std::nullopt;
    }
  }

  if (r_layer.size()) {
    r_layers.push_back(r_layer);
    r_layer.clear();
    row_visi.reset();
  }
  if (c_layer.size()) {
    c_layers.push_back(c_layer);
    c_layer.clear();
    col_visi.reset();
  }

  Mat reduce = origin;
  int sz = 0;  // cnot count
  for (auto [c, t, _] : r_op) reduce = reduce.row_i2j(c, t), sz++;
  for (auto [c, t, _] : c_op) reduce = reduce.col_i2j(c, t), sz++;

  bool ok = true;
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (reduce.m[i][j] != mat.m[i][j]) {
        ok = false;
      }
    }
  }
  if (!ok) return std::nullopt;

  std::vector<int> per(SIZE, 0);
  for (int i = 0; i < SIZE; i++)
    for (int j = 0; j < SIZE; j++)
      if (mat.m[i][j] == 1) per[i] = j;

  std::vector<std::vector<Op>> layers;
  for (auto &l : c_layers) {
    std::vector<Op> nl;
    for (auto [t, c, _] : l) {
      nl.push_back({c, t, 1});
    }
    layers.push_back(nl);
  }
  reverse(r_layers.begin(), r_layers.end());  // 列1-d， 行d到1
  for (auto &l : r_layers) {
    std::vector<Op> nl;
    for (auto [c, t, _] : l) {
      nl.push_back({per[c], per[t], 0});
    }
    layers.push_back(nl);
  }
  return layers;
}

};  // namespace AsiaCrypt24