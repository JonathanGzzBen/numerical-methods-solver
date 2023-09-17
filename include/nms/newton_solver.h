#ifndef LIB_NEWTON_SOLVER_H
#define LIB_NEWTON_SOLVER_H

#include <vector>

namespace nms {

class TableOfValues {
private:
  std::vector<int> x_;
  std::vector<int> y_;

public:
  TableOfValues(std::vector<int> x, std::vector<int> y);
  auto x() -> std::vector<int> &;
  auto y() -> std::vector<int> &;
};

class NewtonForwardInput {
  TableOfValues table_of_values_;
  int x0_; // Desired value

public:
  NewtonForwardInput(TableOfValues table_of_values, int x0);
  auto table_of_values() -> TableOfValues &;
  auto x0() const -> int;
};

class NewtonForwardOutput {
private:
  int h_;
  int x1_;
  std::vector<std::vector<int>> differences_;

public:
  NewtonForwardOutput(int h, int x1, std::vector<std::vector<int>> differences);
  auto get_delta_y(int i) const -> int;
  auto s(double x) const -> double;
  auto gx(double x) const -> double;
};

class NewtonForwardSolver {
public:
  static auto Solve(NewtonForwardInput input) -> NewtonForwardOutput;
};

} // namespace nms

#endif // LIB_NEWTON_SOLVER_H
