
#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "nd_array/nd_array.hpp"

#include <functional>

using real = double;

constexpr real pi = 3.1415926535897932;

static constexpr real min_x = 0.0;
static constexpr real max_x = 1.0;

static constexpr real min_y = 0.0;
static constexpr real max_y = 1.0;

template <int _ctrl_vols_x, int _ctrl_vols_y>
class Mesh {
 public:
  static constexpr int ctrl_vols_x = _ctrl_vols_x;
  static constexpr int ctrl_vols_y = _ctrl_vols_y;

  // Our mesh of control volume has 1 ghost cell on each border
  using ControlVolumes = ND_Array<real, ctrl_vols_x + 2, ctrl_vols_y + 2>;

  static constexpr real dx = (max_x - min_x) / ctrl_vols_x;
  static constexpr real dy = (max_y - min_y) / ctrl_vols_y;

  static constexpr real x1(int cell_x) noexcept {
    return min_x + (cell_x - 1) * dx;
  }

  static constexpr real x2(int cell_x) noexcept { return min_x + cell_x * dx; }

  static constexpr real y1(int cell_y) noexcept {
    return min_y + (cell_y - 1) * dy;
  }

  static constexpr real y2(int cell_y) noexcept { return min_y + cell_y * dy; }

  static constexpr real median_x(int cell_x) noexcept {
    return min_x - dx / 2 + cell_x * dx;
  }

  static constexpr real median_y(int cell_y) noexcept {
    return min_y - dy / 2 + cell_y * dy;
  }

  constexpr Mesh() noexcept {
    for(int i = 0; i < cva.extent(0); i++) {
      for(int j = 0; j < cva.extent(0); j++) {
        cva(i, j) = 0.0;
      }
    }
  }

  ControlVolumes error() noexcept;

  void print() noexcept;
  void print_error() noexcept;

  real accumulate_error(std::function<real(real, real)> accumulator) noexcept;
  real linf_error() noexcept;
  real avg_error() noexcept;
  real rms_error() noexcept;

  real interpolate(real x, real y) noexcept;

 protected:
  static void print_cva(const ControlVolumes &cva) noexcept;
  ControlVolumes cva;
};

#endif // _MESH_HPP_
