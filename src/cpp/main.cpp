#include <cmath>
#include <iostream>
#include <vector>

std::pair<float, float> axialCoord2Euclidean(int q, int r) {
  return {q + r / 2.0, r * std::sqrt(3) / 2.0};
}

std::vector<std::pair<float, float>> hexagon(int n, float k) {
  std::vector<std::pair<float, float>> ret = {{0, 0}};
  if (n == 1) return ret;

  std::vector<std::tuple<int, int, int>> diffs = {{0, 1, -1}, {-1, 1, 0}, {-1, 0, 1},
                                                  {0, -1, 1}, {1, -1, 0}, {1, 0, -1}};

  int distFromCenter = 0;
  while (true) {
    int q = distFromCenter, r = -distFromCenter, s = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < distFromCenter; ++j) {
        ret.push_back(axialCoord2Euclidean(k * q, k * r));
        if (ret.size() == n) return ret;
        q += std::get<0>(diffs[i]);
        r += std::get<1>(diffs[i]);
        s += std::get<2>(diffs[i]);
      }
    }
    distFromCenter += 1;
  }
}

int main() {
  int n = 61;
  float k = 1.0;
  std::vector<std::pair<float, float>> hexagonList = hexagon(n, k);

  for (const auto& coord : hexagonList) {
    std::cout << "(" << coord.first << ", " << coord.second << ")\n";
  }

  return 0;
}