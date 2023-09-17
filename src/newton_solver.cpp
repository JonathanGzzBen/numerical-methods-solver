#include "nms/newton_solver.h"
#include <stdexcept>

namespace nms {

TableOfValues::TableOfValues(std::vector<int> x, std::vector<int> y)
    : x_{x}, y_{y} {
  if (x.size() != y.size())
    throw std::domain_error{"x and y vectors do not have the same length"};
}

auto TableOfValues::x() const -> const std::vector<int> & { return x_; }
auto TableOfValues::y() const -> const std::vector<int> & { return y_; }

NewtonInterpolationSolution::NewtonInterpolationSolution(
    int h, int x1, std::vector<std::vector<int>> differences)
    : h_{h}, x1_{x1}, differences_{differences} {}

auto NewtonInterpolationSolution::s(double x) const -> double {
  return static_cast<double>(x - x1_) / h_;
}

auto NewtonInterpolationSolution::get_delta_y(int i) const -> int {
  return differences_[i][0];
}

auto NewtonInterpolationSolution::gx(double x) const -> double {
  const auto &factorial = [](int n) {
    if (n == 0)
      return 1;
    int res = n;
    for (n = n - 1; n > 1; n--)
      res *= n;
    return res;
  };

  double gx = differences_[0][0];
  double s_value = s(x);
  double total_s = s_value;
  for (int i = 1; i < differences_.size(); i++) {
    gx += (differences_[i][0] * total_s) / factorial(i);
    total_s *= (--s_value);
  }
  return gx;
};

auto NewtonInterpolationSolver::Solve(const TableOfValues &table_of_values)
    -> NewtonInterpolationSolution {
  if (table_of_values.x().size() < 2)
    throw std::domain_error{"x size less than 2"};
  const int h = std::abs(table_of_values.x()[1] - table_of_values.x()[0]);
  for (size_t i = 1; i < table_of_values.x().size() - 1; i++) {
    if (std::abs(table_of_values.x()[i + 1] - table_of_values.x()[i]) != h)
      throw std::domain_error{"x differences are not equally distributed"};
  }

  auto differences = std::vector<std::vector<int>>(table_of_values.y().size());
  differences[0] = std::vector<int>(table_of_values.y().size());
  for (size_t i = 0; i < table_of_values.y().size(); i++)
    differences[0][i] = table_of_values.y()[i];

  const auto &calculate_differences = [&](size_t difference_i) {
    differences[difference_i] =
        std::vector<int>(differences[difference_i - 1].size() - 1);
    for (unsigned int row_i = 0; row_i < differences[difference_i].size();
         row_i++) {
      differences[difference_i][row_i] =
          differences[difference_i - 1][row_i + 1] -
          differences[difference_i - 1][row_i];
    }
  };

  for (size_t i = 1; i < table_of_values.y().size(); i++) {
    calculate_differences(i);
  }
  return {h, table_of_values.x()[0], std::move(differences)};
}

} // namespace nms
