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
  auto x() const -> const std::vector<int> &;
  auto y() const -> const std::vector<int> &;
};

class NewtonInterpolationSolution {
private:
  int h_;
  int x1_;
  std::vector<std::vector<int>> differences_;

  auto s(double x) const -> double;

public:
  NewtonInterpolationSolution(int h, int x1,
                              std::vector<std::vector<int>> differences);
  auto get_delta_y(int i) const -> int;
  auto gx(double x) const -> double;
};

class NewtonInterpolationSolver {
public:
  static auto Solve(const TableOfValues &table_of_values)
      -> NewtonInterpolationSolution;
};

} // namespace nms

#endif // LIB_NEWTON_SOLVER_H
