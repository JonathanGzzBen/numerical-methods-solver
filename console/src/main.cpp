#include <iostream>

#include "nms/newton_solver.h"

void newton_solve() {
  std::vector<int> x = {-3, -1, 1, 3};
  std::vector<int> y = {125, 127, 132, 133};
  int x0 = 0;
  std::cout << "newton solving\n";
  nms::NewtonForwardInput input({x, y}, x0);
  try {
    auto output = nms::NewtonForwardSolver::Solve(input);
    std::cout << "g(0) = " << output.gx(0) << "\n";
    std::cout << "g(0) = " << output.gx(0) << "\n";
    std::cout << "g(1) = " << output.gx(1) << "\n";
    std::cout << "g(-1) = " << output.gx(-1) << "\n";
  } catch (const std::domain_error &e) {
    std::cerr << e.what();
  }
}

int main() {
  newton_solve();
  return 0;
}
