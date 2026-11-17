#include "shape_obb.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main()
{
  std::size_t n;
  if (!(std::cin >> n)) return 2;
  std::vector<LAMMPS_NS::ShapeOBB::Point> points(n);
  for (auto &point : points)
    for (double &value : point) {
      std::string token;
      if (!(std::cin >> token)) return 2;
      char *end;
      value = std::strtod(token.c_str(), &end);
      if (*end != '\0') return 2;
    }
  const auto result = LAMMPS_NS::ShapeOBB::min_volume(points);
  std::cout << result.valid << std::setprecision(17);
  for (double value : result.length) std::cout << ' ' << value;
  std::cout << '\n';
}
