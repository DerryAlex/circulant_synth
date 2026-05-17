#include <bits/stdc++.h>

#include "greedy.h"
#include "reduce.h"

enum class Mode {
  depth,
  size,
};
Mode mode = Mode::size;

std::map<std::string, int> depth_record = {
    {"AES.txt", 10},
    {"Anubis.txt", 10},
    {"C_BeiKraLea16_4x4_4.txt", 11},
    {"C_BeiKraLea16_4x4_8.txt", 18},
    {"Clefia_M0.txt", 10},
    {"Clefia_M1.txt", 11},
    {"FSE_LiWang16_4x4_4.txt", 10},
    {"FSE_LiWang16_4x4_8.txt", 16},
    {"FSE_LiWang16_i_4x4_4.txt", 9},
    {"FSE_LiWang16_i_4x4_8.txt", 8},
    {"FSE_LiWang16_i_4x4_8_2.txt", 11},
    {"FSE_LiuSim16_4x4_4.txt", 10},
    {"FSE_LiuSim16_4x4_8.txt", 17},
    {"FSE_SKOP15_4x4_4.txt", 8},
    {"FSE_SKOP15_4x4_8.txt", 9},
    {"FSE_SKOP15_8x8_4.txt", 20},
    {"FSE_SKOP15_i_4x4_4.txt", 8},
    {"FSE_SKOP15_i_4x4_8.txt", 8},
    {"FSE_SKOP15_i_8x8_4.txt", 23},
    {"Joltik.txt", 8},
    {"M_4_4.txt", 7},
    {"M_4_8.txt", 7},
    {"SM4_L=83.txt", 9},
    {"SmallScale_AES.txt", 8},
    {"Whirlwind_M0.txt", 17},
    {"Whirlwind_M1.txt", 19},
};

std::map<std::string, int> size_record = {
    {"AES.txt", 105},
    {"Anubis.txt", 116}, // 111
    {"C_BeiKraLea16_4x4_4.txt", 44},
    {"C_BeiKraLea16_4x4_8.txt", 147},
    {"Clefia_M0.txt", 115}, // 111
    {"Clefia_M1.txt", 117}, // 115
    {"FSE_LiWang16_4x4_4.txt", 52},
    {"FSE_LiWang16_4x4_8.txt", 130},
    {"FSE_LiWang16_i_4x4_4.txt", 48},
    {"FSE_LiWang16_i_4x4_8.txt", 100},
    {"FSE_LiWang16_i_4x4_8_2.txt", 100},
    {"FSE_LiuSim16_4x4_4.txt", 55},
    {"FSE_LiuSim16_4x4_8.txt", 144},
    {"FSE_SKOP15_4x4_4.txt", 50},
    {"FSE_SKOP15_4x4_8.txt", 110},
    {"FSE_SKOP15_8x8_4.txt", 176},
    {"FSE_SKOP15_i_4x4_4.txt", 48},
    {"FSE_SKOP15_i_4x4_8.txt", 99},
    {"FSE_SKOP15_i_8x8_4.txt", 186},
    {"Joltik.txt", 50}, // 48
    {"MIDORI.txt", 24},
    {"M_4_4.txt", 42},
    {"M_4_8.txt", 91},
    {"PRIDE_L_0.txt", 24},
    {"PRIDE_L_3.txt", 24},
    {"PRINCE_M_0.txt", 24},
    {"PRINCE_M_1.txt", 24},
    {"QARMA128.txt", 48},
    {"QARMA64.txt", 24},
    {"SM4_L=83.txt", 96},
    {"SmallScale_AES.txt", 50},
    {"Whirlwind_M0.txt", 159},
    {"Whirlwind_M1.txt", 169},
};

Matrix<int> read_matrix(std::istream &is) {
  int n, x, y;
  is >> n;
  assert(n == 1);
  is >> x >> y;
  Matrix<int> matrix(x, std::vector<int>(y, 0));
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      is >> matrix[i][j];
    }
  }
  return matrix;
}

std::tuple<int, int> calc_depth_and_size(
    const std::vector<std::vector<Op>> &ops) {
  int depth = 0, size = 0;
  for (auto layer : ops) {
    depth += 1;
    size += layer.size();
  }
  if (depth == 0) depth = 100;
  if (size == 0) size = 400;
  return std::make_pair(depth, size);
}

std::vector<std::vector<Op>> process(Mat mat, double UPDATE_TIME_LIMIT = 1.0,
                                     double TIME_LIMIT = 20.0) {
  std::vector<std::vector<Op>> best;
  int best_depth = INT_MAX, best_size = INT_MAX;
  clock_t start_time = clock();
  clock_t update_time = start_time;
  while (clock() < update_time + UPDATE_TIME_LIMIT * CLOCKS_PER_SEC &&
         clock() < start_time + TIME_LIMIT * CLOCKS_PER_SEC) {
    std::vector<std::vector<Op>> result;
    std::optional<std::vector<std::vector<Op>>> result_;
    switch (mode) {
      case Mode::depth:
        result_ = AsiaCrypt24::greedy(mat);
        if (result_.has_value()) {
          result = *result_;
        }
        break;
      case Mode::size:
        result = FSE20::reduce(mat);
        break;
    }
    if (!result.empty()) {
      auto [depth, size] = calc_depth_and_size(result);
      if ((mode == Mode::depth &&
           (depth < best_depth || (depth == best_depth && size < best_size))) ||
          (mode == Mode::size &&
           (size < best_size || (size == best_size && depth < best_depth)))) {
        best = result;
        best_depth = depth;
        best_size = size;
        update_time = clock();
      }
    }
  }
  if (best_depth == INT_MAX) {
    std::cerr << "greedy failed" << std::endl;
  }
  return best;
}

std::optional<int> check_circulant(const Matrix<int> &matrix);
std::vector<std::pair<Matrix<int>, int>> process_circulant(Matrix<int> matrix,
                                                           int block_size = 1);

int main(int argc, char *argv[]) {
  if (argc != 2 || (argc == 2 && (strcmp("-h", argv[1]) == 0 ||
                                  strcmp("--help", argv[1]) == 0))) {
    std::cerr << std::format("Usage: {} path_to_matrices", argv[0])
              << std::endl;
    return (argc == 2) ? 0 : 1;
  }
  srand(time(NULL));
  std::string path = argv[1];
  for (const auto &entry : std::filesystem::directory_iterator(path)) {
    std::ifstream file(entry.path());
    Matrix<int> matrix = read_matrix(file);
    int n = matrix.size();
    auto _block_size = check_circulant(matrix);
    if (!_block_size.has_value()) continue;
    std::string filename = entry.path().filename();
    if (!filename.starts_with("PRIDE")) continue;
    int baseline_depth = 100;
    if (depth_record.contains(filename)) {
      baseline_depth = depth_record[filename];
    }
    int baseline_size = 400;
    if (size_record.contains(filename)) {
      baseline_size = size_record[filename];
    }
    for (int block_size = *_block_size; block_size <= n; block_size *= 2) {
      int depth_overhead = 2 * (int)std::log2(n / block_size);
      int size_overhead = depth_overhead / 2 * n;
      auto _matrix = process_circulant(matrix, block_size);
      for (auto [m, id] : _matrix) {
        std::vector<std::vector<Op>> result =
            process(m, 20, 60);
        auto [depth, size] = calc_depth_and_size(result);
        std::cout << std::format("{} {} {} ({}+{} {}+{} id={})",
                                 entry.path().filename().string(),
                                 depth + depth_overhead, size + size_overhead,
                                 depth, depth_overhead, size, size_overhead, id)
                  << std::endl;
        depth += depth_overhead;
        size += size_overhead;
        if ((mode == Mode::depth && depth <= baseline_depth) ||
            (mode == Mode::size && size <= baseline_size)) {
          for (auto layer : result) {
            for (auto op : layer) {
              std::cout << std::format("{} {} ", op.src, op.dst);
            }
            std::cout << std::endl;
          }
        }
      }
    }
  }
  return 0;
}

std::optional<int> check_circulant(const Matrix<int> &matrix) {
  int n = matrix.size();
  if (matrix[0].size() != (size_t)n) return false;
  for (int step = 1; step < n; step *= 2) {
    if (n % step != 0) break;
    bool flag = true;
    for (int i = 0; i < n / step && flag; i++) {
      int i_prev = i - 1;
      if (i_prev < 0) i_prev += n / step;
      for (int j = 0; j < n / step && flag; j++) {
        int j_prev = j - 1;
        if (j_prev < 0) j_prev += n / step;
        for (int ii = 0; ii < step && flag; ii++) {
          for (int jj = 0; jj < step && flag; jj++) {
            if (matrix[i * step + ii][j * step + jj] !=
                matrix[i_prev * step + ii][j_prev * step + jj]) {
              flag = false;
            }
          }
        }
      }
    }
    if (flag) return step;
  }
  return std::nullopt;
}

std::vector<std::pair<Matrix<int>, int>> process_circulant_ex(
    Matrix<int> _matrix, int block_size, int step) {
  int n = _matrix.size();
  std::vector<std::pair<Matrix<int>, int>> result;
  for (int S = 0; S < 1 << (n / step); S++) {
    std::stringstream ss;
    ss << std::format("{:0{}b}", S, n / step);
    auto matrix = _matrix;
    for (int i = 0; i < n / 2 / step; i++) {
      char type;
      ss >> type;
      for (int j = 0; j < step; j++) {
        if (type == '0') {
          for (int k = 0; k < n; k++) {
            matrix[i * step * 2 + j][k] ^= matrix[i * step * 2 + j + step][k];
          }
        } else {
          for (int k = 0; k < n; k++) {
            matrix[i * step * 2 + j + step][k] ^= matrix[i * step * 2 + j][k];
          }
        }
      }
    }
    for (int i = 0; i < n / 2 / step; i++) {
      char type;
      ss >> type;
      for (int j = 0; j < step; j++) {
        if (type == '0') {
          for (int k = 0; k < n; k++) {
            matrix[k][i * step * 2 + j] ^= matrix[k][i * step * 2 + j + step];
          }
        } else {
          for (int k = 0; k < n; k++) {
            matrix[k][i * step * 2 + j + step] ^= matrix[k][i * step * 2 + j];
          }
        }
      }
    }
    if (step > block_size) {
      for (auto [m, id] : process_circulant_ex(matrix, block_size, step / 2)) {
        result.push_back({m, S * 1000 + id});
      }
    } else {
      result.push_back({matrix, S});
    }
  }
  return result;
}

std::vector<std::pair<Matrix<int>, int>> process_circulant(Matrix<int> matrix,
                                                           int block_size) {
  int n = matrix.size();
  if (__builtin_popcount(n) != 1) return {{matrix, 0}};
  if (block_size >= n) return {{matrix, 0}};
  if (n / block_size <= 4)
    return process_circulant_ex(matrix, block_size, n / 2);
  for (int step = n / 2; step >= block_size; step /= 2) {
    for (int i = 0; i < n / 2 / step; i++) {
      for (int j = 0; j < step; j++) {
        for (int k = 0; k < n; k++) {
          matrix[i * step * 2 + j + step][k] ^= matrix[i * step * 2 + j][k];
        }
      }
    }
    for (int i = 0; i < n / 2 / step; i++) {
      for (int j = 0; j < step; j++) {
        for (int k = 0; k < n; k++) {
          matrix[k][i * step * 2 + j] ^= matrix[k][i * step * 2 + j + step];
        }
      }
    }
  }
  return {{matrix, 0}};
}
