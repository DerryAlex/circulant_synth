#include "reduce.h"
#include "greedy.h"
#include <algorithm>

namespace FSE20 {

// --- strgy begin

std::mt19937 rand_generator(time(NULL));

std::vector<int> build_table(const Mat& m) {
  size_t SIZE = m.size();
  std::vector<int> tab(SIZE, -1);
  for (unsigned i = 0; i < SIZE; ++i) {
    int ind = -1;
    const auto& layer = m[i];
    for (unsigned j = 0; j < SIZE; ++j) {
      if (layer[SIZE - 1 - j]) {
        ind = j;
        break;
      }
    }
    tab.at(ind) = i;
  }
  for (auto v : tab) {
    if (v == -1) throw std::out_of_range("expected permutation matrix");
  }
  return tab;
}

std::vector<Op> update_seq_str(const std::vector<Op>& seq,
                               const std::vector<int>& tab) {
  std::vector<Op> tmp_seq;
  for (const auto& op : seq) {
    if (!op.type) {
      tmp_seq.push_back(op);
    }
  }

  for (const auto& op : seq | std::views::reverse) {
    if (op.type) {
      tmp_seq.push_back({tab[op.dst], tab[op.src], 0});
    }
  }
  return tmp_seq;
}

void select_oper(const Mat& m, std::vector<Op>& max_seq, int& no_reduced,
                 int opr_type) {
  size_t SIZE = m.size();
  int no_before = 0;
  int no_after = 0;

  for (unsigned i = 0; i < SIZE; i++) {
    for (unsigned j = 0; j < SIZE; j++) {
      if (i != j) {
        auto tmp = m[j];
        no_before = tmp.count();
        tmp ^= m[i];
        no_after = tmp.count();
        if ((no_before - no_after) > 0) {
          if (no_reduced < (no_before - no_after)) {
            no_reduced = no_before - no_after;
            max_seq.clear();
          }
          if (no_reduced == (no_before - no_after)) {
            max_seq.push_back({(int)i, (int)j, opr_type});
          }
        }
      }
    }
  }
}

Mat get_trans_matrix(const Mat& m) {
  size_t SIZE = m.size();
  Mat trans_m(SIZE);
  for (unsigned i = 0; i < SIZE; i++) {
    for (unsigned j = 0; j < SIZE; j++) {
      trans_m[i][SIZE - 1 - j] = m[j][SIZE - i - 1];
    }
  }
  return trans_m;
}

std::vector<Op> strgy1(Mat& m, int opr_type = 0) {
  /* compatible with non-square matrix.*/
  std::vector<Op> seq;
  size_t row_size = m.size();
  size_t col_size = m[0].size();

  dynamic_bitset mark(row_size);
  for (unsigned col = 1; col <= col_size; ++col) {
    unsigned r = 0;
    while (((r < row_size) && (!m[r][col_size - col])) || (mark[r] == 1)) {
      ++r;
    }
    if (r >= row_size) {
      continue;
    } else {
      mark[r] = 1;
    }

    for (unsigned i = 0; i < row_size; ++i) {
      if (m[i][col_size - col] && (i != r)) {
        m[i] ^= m[r];
        seq.push_back({(int)r, (int)i, opr_type});
      }
    }
  }
  return seq;
}

std::vector<Op> strgy2(Mat& m) {
  Mat trans_m = get_trans_matrix(m);
  std::vector<Op> seq = strgy1(trans_m, 1);
  m = get_trans_matrix(trans_m);

  std::vector<int> tab = build_table(m);  // get TABLE
  std::vector<Op> final_seq = update_seq_str(seq, tab);
  return final_seq;
}

std::vector<Op> strgy3(Mat& m) {
  std::vector<Op> tmp_seq;
  while (m.count() != m.size()) {
    std::vector<Op> base_oper;
    int no_reduced = 0;
    select_oper(m, base_oper, no_reduced, 0);
    Mat trans_m = get_trans_matrix(m);
    select_oper(trans_m, base_oper, no_reduced, 1);

    if (base_oper.size() >= 1) {
      int rand_num = rand_generator() % base_oper.size();
      const auto& op = base_oper[rand_num];
      tmp_seq.push_back(op);

      if (!op.type) {
        m[op.dst] ^= m[op.src];
      } else if (op.type) {
        auto trans_tmp_m = get_trans_matrix(m);
        trans_tmp_m[op.dst] ^= trans_tmp_m[op.src];
        m = get_trans_matrix(trans_tmp_m);
      }
    } else {
      break;
    }
  }

  if (m.count() != m.size()) {
    // int rnd = rand()%2;
    int rnd = rand_generator() % 2;
    if (rnd == 0) {
      std::vector<Op> seq_2(strgy1(m));
      tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
    } else if (rnd == 1) {
      std::vector<Op> seq_2(strgy2(m));
      tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
    }
  }

  std::vector<int> tab = build_table(m);
  std::vector<Op> final_seq = update_seq_str(tmp_seq, tab);
  return final_seq;
}

// --- strgy end

// --- reduce begin

bool exchange(Op a, Op b) {
  if (a.dst == b.dst) return true;
  if (a.src == b.src) return true;
  if ((a.src != b.dst) && (a.dst != b.src)) return true;
  return false;
}

bool exchange_set(const std::vector<Op>& seq, int start, int end, Op p) {
  if (start <= end) {
    for (int i = start; i <= end; i++) {
      if (!exchange(p, seq[i])) return false;
    }
  } else {
    for (int i = start; i >= end; i--) {
      if (!exchange(p, seq[i])) return false;
    }
  }
  return true;
}

Matrix<int> get_table(const std::vector<Op>& seq) {
  int nsize = seq.size();
  Matrix<int> table(nsize, std::vector<int>(nsize, 0));

  int s = seq.size();
  for (int i = 0; i < s; i++) {
    table[i][i] = 1;
    for (int j = i + 1; j < s; j++) {
      if (exchange(seq[i], seq[j]))
        table[i][j] = 1;
      else {
        for (int k = j; k < s; k++) table[i][k] = 0;
        break;
      }
    }
    for (int j = i - 1; j >= 0; j--) {
      if (exchange(seq[i], seq[j]))
        table[i][j] = 1;
      else {
        for (int k = j; k >= 0; k--) table[i][k] = 0;
        break;
      }
    }
  }
  return table;
}

bool reduce0(std::vector<Op>& seq, const Matrix<int>& table);
bool reduce1(std::vector<Op>& seq, const Matrix<int>& table);
bool reduce2(std::vector<Op>& seq, const Matrix<int>& table);
bool reduce3(std::vector<Op>& seq, const Matrix<int>& table);

int reduce_step(std::vector<Op>& seq) {
  auto tab = get_table(seq);
  int i = 0;
  int NUM = 4;
  int counter = 0;
  while (counter != NUM) {
    switch (i) {
      case 0: {
        if (reduce0(seq, tab)) {
          // cout << "activate case 0"<<endl;
          tab = get_table(seq);
          // i = 0;
          counter = 0;
        } else {
          i = (i + 1) % NUM;
          counter++;
        }
        break;
      }
      case 1: {
        if (reduce1(seq, tab)) {
          // cout << "activate case 1" << endl;
          tab = get_table(seq);
          // i = 1;
          counter = 0;
        } else {
          i = (i + 1) % NUM;
          counter++;
        }
        break;
      }
      case 2: {
        if (reduce2(seq, tab)) {
          // cout << "activate case 2" <<endl;
          tab = get_table(seq);
          counter = 0;
        } else {
          i = (i + 1) % NUM;
          counter++;
        }
        break;
      }
      case 3: {
        if (reduce3(seq, tab)) {
          // cout << "activate case 3" << endl;
          tab = get_table(seq);
          counter = 0;
        } else {
          // cout << "case 3" << endl;
          i = (i + 1) % NUM;
          counter++;
        }
        break;
      }
    }
  }
  return seq.size();
}

void update_seq(std::vector<Op>& seq, const std::vector<int>& tab, int start) {
  for (auto& op : std::views::drop(seq, start)) {
    op.src = tab[op.src];
    op.dst = tab[op.dst];
  }
}

void get_equivalent_seq(int SIZE, std::vector<Op>& seq, int gap, int start) {
  Mat m(SIZE);
  for (int i = 0; i < SIZE; i++) {
    m[i][SIZE - 1 - i] = 1;
  }
  for (int i = start + gap - 1; i >= start; i--) {
    m[seq[i].dst] ^= m[seq[i].src];
  }

  std::vector<Op> seq_here(strgy3(m));
  seq.erase(seq.begin() + start, seq.begin() + start + gap);
  seq.insert(seq.begin() + start, seq_here.begin(), seq_here.end());
  auto tab = build_table(m);
  update_seq(seq, tab, start + seq_here.size());
}

typedef struct {
  std::vector<Op> seq;
  int gap;
  int start;
  int len;
  int SIZE;
} thread_data;

void* reduce_thread(void* d) {
  thread_data* data = (thread_data*)d;
  get_equivalent_seq(data->SIZE, data->seq, data->gap, data->start);
  data->len = reduce_step(data->seq);
  pthread_exit(NULL);
}

Mat get_reduced_matrix(const std::vector<Op>& seq, const Mat& m) {
  Mat tmp_m(m);
  for (const auto& op : seq) tmp_m[op.dst] ^= tmp_m[op.src];
  return tmp_m;
}

std::vector<Op> strgy4(Mat &m) {
  auto result = AsiaCrypt24::greedy(m);
  if (!result.has_value()) return strgy3(m);
  std::vector<Op> seq;
  for (const auto &layer: *result) {
    for (const auto &op: layer) {
      seq.push_back(op);
    }
  }
  std::reverse(seq.begin(), seq.end());
  return seq;
}

std::vector<Op> do_reduce(Mat m) {
  int SIZE = m.size();
  std::vector<Op> seq = strgy3(m);
  int l = reduce_step(seq);

  int gap = l;
  int start = 0;

  constexpr int THREAD_NUM = 4;
  thread_data data[THREAD_NUM];
  bool flag = true;
  while (flag) {
    for (int j = 0; j < THREAD_NUM; j++) {
      data[j].SIZE = SIZE;
      data[j].seq = seq;
      data[j].gap = gap;
      data[j].start = start;
      if ((gap != 3) || (start != l - gap)) {
        if (start == l - gap) {
          gap--;
          start = 0;
        } else {
          start++;
        }
      } else {
        flag = false;
      }
    }
    int rc;
    pthread_t threads[THREAD_NUM];
    void* status;
    for (int j = 0; j < THREAD_NUM; j++) {
      rc = pthread_create(&threads[j], NULL, reduce_thread, (void*)(data + j));
      if (rc) {
        std::cout << "Error: unable to create thread!" << std::endl;
        exit(-1);
      }
    }
    // pthread_attr_destroy(&attr);
    for (int j = 0; j < THREAD_NUM; j++) {
      rc = pthread_join(threads[j], &status);
      if (rc) {
        std::cout << "Error: unable to join, " << rc << std::endl;
        exit(-1);
      }
    }
    for (int j = 0; j < THREAD_NUM; j++) {
      if (data[j].len < l) {
        seq = data[j].seq;
        l = data[j].len;
        gap = l;
        start = 0;
      }
    }
    // cout << "gap = " << gap << endl;
  }

  return seq;
}

std::vector<Op> reduce_main(const Mat& mat) {
  int SIZE = mat.size();
  auto tmp_m = mat;
  auto seq = do_reduce(mat);
  auto tmp = get_reduced_matrix(seq, mat);
  std::vector<int> tab(SIZE, -1);
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      if (tmp[i][j]) {
        tab[i] = j;
        break;
      }
    }
  }
  update_seq(seq, tab, 0);
  return seq;
}

// --- reduce end

std::vector<std::vector<Op>> reduce(Mat mat) {
  size_t SIZE = mat.size();
  auto seq = reduce_main(mat);
  std::reverse(seq.begin(), seq.end());

  dynamic_bitset vis(SIZE);
  std::vector<std::vector<Op>> result;
  std::vector<Op> layer;
  for (auto op : seq) {
    if (vis[op.src] || vis[op.dst]) {
      result.push_back(layer);
      layer.clear();
      vis.reset();
    }
    vis[op.src] = true;
    vis[op.dst] = true;
    layer.push_back(op);
  }
  if (!layer.empty()) result.push_back(layer);
  return result;
}

bool reduce0(std::vector<Op>& seq, const Matrix<int>& table) {
  int s = seq.size();
  int index;

  for (int i = 0; i < s; i++) {
    for (int j = i + 1; j < s; j++) {
      if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0)) break;
      if (seq[i].src != seq[j].src) continue;
      for (int k = j + 1; k < s; k++) {
        if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0)) break;
        if (seq[k].dst != seq[i].dst) {
          if (seq[k].dst != seq[j].dst)
            continue;
          else {
            if (seq[k].src != seq[i].dst) continue;
          }
        } else {
          if (seq[k].src != seq[j].dst) continue;
        }
        if ((table[i][j - 1]) && table[k][j + 1])
          index = j - 1;  // fixpoint = j;
        else if (table[j][i + 1] && table[k][j + 1] &&
                 exchange_set(seq, j - 1, i + 1, seq[k]))
          index = i;  // fixpoint = i;
        else if (table[i][j - 1] && table[j][k - 1] &&
                 exchange_set(seq, j + 1, k - 1, seq[i]))
          index = k - 2;  // fixpoint = k;
        else
          continue;
        Op tmp1, tmp2;
        tmp1 = seq[k];
        tmp2.src = seq[i].src;
        tmp2.dst = seq[k].src;
        seq.erase(seq.begin() + k);
        seq.erase(seq.begin() + j);
        seq.erase(seq.begin() + i);
        seq.insert(seq.begin() + index, tmp1);
        seq.insert(seq.begin() + index + 1, tmp2);
        return true;
      }
    }
  }
  return false;
}

bool reduce1(std::vector<Op>& seq, const Matrix<int>& table) {
  int s = seq.size();
  int index;

  for (int i = 0; i < s; i++) {
    for (int j = i + 1; j < s; j++) {
      if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0)) break;
      if ((seq[i].src != seq[j].dst) && (seq[i].dst != seq[j].dst)) continue;
      for (int k = j + 1; k < s; k++) {
        if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0)) break;
        if (seq[k].src != seq[j].src)
          continue;
        else {
          if (seq[j].dst == seq[i].dst) {
            if (seq[k].dst != seq[i].src) continue;
          }
          if (seq[j].dst == seq[i].src) {
            if (seq[k].dst != seq[i].dst) continue;
          }
        }
        if ((table[i][j - 1]) && table[k][j + 1])
          index = j - 1;  // fixpoint = j;
        else if (table[j][i + 1] && table[k][j + 1] &&
                 exchange_set(seq, j - 1, i + 1, seq[k]))
          index = i;  // fixpoint = i;
        else if (table[i][j - 1] && table[j][k - 1] &&
                 exchange_set(seq, j + 1, k - 1, seq[i]))
          index = k - 2;  // fixpoint = k;
        else
          continue;
        Op tmp1, tmp2;
        tmp1.src = seq[j].src;
        tmp1.dst = seq[i].src;
        tmp2 = seq[i];
        seq.erase(seq.begin() + k);
        seq.erase(seq.begin() + j);
        seq.erase(seq.begin() + i);
        seq.insert(seq.begin() + index, tmp1);
        seq.insert(seq.begin() + index + 1, tmp2);
        return true;
      }
    }
  }
  return false;
}

bool reduce2(std::vector<Op>& seq, const Matrix<int>& table) {
  int s = seq.size();
  int index;

  for (int i = 0; i < s; i++) {
    for (int j = i + 1; j < s; j++) {
      if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0)) break;
      if ((seq[j].src != seq[i].dst) && (seq[j].dst != seq[i].dst)) continue;
      for (int k = j + 1; k < s; k++) {
        if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0)) break;
        if (seq[k].src != seq[i].src)
          continue;
        else {
          if (seq[j].src == seq[i].dst) {
            if (seq[k].dst != seq[j].dst) continue;
          }
          if (seq[j].dst == seq[i].dst) {
            if (seq[k].dst != seq[j].src) continue;
          }
        }
        if ((table[i][j - 1]) && table[k][j + 1])
          index = j - 1;  // fixpoint = j;
        else if (table[j][i + 1] && table[k][j + 1] &&
                 exchange_set(seq, j - 1, i + 1, seq[k]))
          index = i;  // fixpoint = i;
        else if (table[i][j - 1] && table[j][k - 1] &&
                 exchange_set(seq, j + 1, k - 1, seq[i]))
          index = k - 2;  // fixpoint = k;
        else
          continue;
        Op tmp1, tmp2;
        tmp1.dst = seq[k].dst;
        tmp2.dst = seq[i].dst;
        if (seq[j].src == seq[i].dst) {
          tmp1.src = seq[j].src;
          tmp2.src = seq[i].src;
        } else {
          tmp1.src = seq[i].src;
          tmp2.src = seq[j].src;
        }
        seq.erase(seq.begin() + k);
        seq.erase(seq.begin() + j);
        seq.erase(seq.begin() + i);
        seq.insert(seq.begin() + index, tmp1);
        seq.insert(seq.begin() + index + 1, tmp2);
        return true;
      }
    }
  }
  return false;
}

bool reduce3(std::vector<Op>& seq, const Matrix<int>& table) {
  int s = seq.size();
  int index;

  for (int i = 0; i < s; i++) {
    for (int j = i + 1; j < s; j++) {
      if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0)) break;
      if (seq[j].dst != seq[i].src) continue;
      if (seq[j].src != seq[i].dst) continue;
      if (table[i][j - 1])
        index = j - 1;
      else
        index = i;
      Op tmp = seq[j];
      seq.erase(seq.begin() + j);
      seq.erase(seq.begin() + i);
      seq.insert(seq.begin() + index, tmp);
      int mask = seq[index].dst ^ seq[index].src;
      for (int k = index + 1; k < s; k++) {
        if ((seq[k].dst == seq[index].dst) || (seq[k].dst == seq[index].src))
          seq[k].dst ^= mask;
        if ((seq[k].src == seq[index].dst) || (seq[k].src == seq[index].src))
          seq[k].src ^= mask;
      }
      return true;
    }
  }
  return false;
}

};  // namespace FSE20
