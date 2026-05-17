// The original program is written by Haotian Shi under GPLv3 License
//
// It is modified to support
// - use as library
// - dynamic size of matrix
// - H_sum cost function
#pragma once

#include <bits/stdc++.h>

#include "matrix.h"

namespace AsiaCrypt24 {

std::optional<std::vector<std::vector<Op>>> greedy(Mat mat,
                                                   int depth_limit = 100,
                                                   int size_limit = 400);

};  // namespace AsiaCrypt24
