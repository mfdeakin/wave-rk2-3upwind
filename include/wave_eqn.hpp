
#ifndef _WAVE_EQN_HPP_
#define _WAVE_EQN_HPP_

#include <assert.h>
#include <cmath>

#include "nd_array/nd_array.hpp"

using real = double;

constexpr real pi = 3.1415926535897932;

static constexpr real min_x = 0.0;
static constexpr real max_x = 1.0;

class Part1 {
 public:
  template <typename MeshT>
  static void initial_conds(MeshT &mesh) noexcept {
    for(int i = 0; i < mesh.extent(0); i++) {
      mesh(i) =
          (std::cos(2.0 * pi * mesh.x2(i)) - std::cos(2.0 * pi * mesh.x1(i))) /
          (2.0 * pi * mesh.dx);
    }
  }

  static constexpr real boundary_val(const real time) noexcept {
    return std::sin(4.0 * pi * time);
  }

  template <typename MeshT>
  static void fill_ghostcells(MeshT &mesh, real time) noexcept {
    mesh(0) = 2.0 * boundary_val(time) - mesh(1);
  }
};

template <int _ctrl_vols, typename _boundary_conds = Part1>
class Mesh : public ND_Array<real, _ctrl_vols + 2> {
 public:
  using BoundaryConds = _boundary_conds;

  static constexpr real dx = (max_x - min_x) / _ctrl_vols;

  static constexpr real x1(const int cell_x) noexcept {
    return min_x + (cell_x - 1) * dx;
  }

  static constexpr real x2(const int cell_x) noexcept {
    return min_x + cell_x * dx;
  }

  static constexpr real median_x(const int cell_x) noexcept {
    return min_x - dx / 2 + cell_x * dx;
  }

  static constexpr int cell_idx(const real x) noexcept {
    return static_cast<int>((x - min_x) / dx) + 1;
  }

  constexpr Mesh() noexcept {
    for(int i = 0; i < this->extent(0); i++) {
      (*this)(i) = 0.0;
    }
  }

  // Computes the CVA of dT/dx
  constexpr real flux_integral(const int i, const real time) const noexcept {
    assert(i > 0);
    assert(i < this->extent(0) - 1);
    if(i > 1) {
      return (3.0 * (*this)(i)-4.0 * (*this)(i - 1) + (*this)(i - 2)) /
             (2.0 * dx);
    } else {
      return (3.0 * (*this)(i)-1.0 * (*this)(i - 1) -
              2.0 * BoundaryConds::boundary_val(time)) /
             (2.0 * dx);
    }
  }
};

enum class TimeStage { stage_0 = 0, stage_1 = 1, stage_partial = 2 };

template <int _ctrl_vols, typename _boundary_conds>
class WaveEqnSolver {
 public:
  using BoundaryConds = _boundary_conds;
  using MeshT         = Mesh<_ctrl_vols, BoundaryConds>;

  static constexpr int time_stages = 3;

  constexpr WaveEqnSolver(const real dtdx = 1.0) noexcept
      : dt(dtdx * cva(0).dx), t(0.0), cur_ts(TimeStage::stage_0) {
    BoundaryConds::initial_conds(cur_mesh());
  }

  constexpr void fill_ghostcells() noexcept {
    BoundaryConds::fill_ghostcells(cur_mesh(), time());
  }

  void timestep_rk2() noexcept;
  void timestep_rk3() noexcept;

  real time() const noexcept { return t; }

  const real dt;

  // Performs a linear interpolation between the cells to compute x
  real operator()(const real x) const noexcept;

  real operator[](int i) const noexcept {
    assert(i >= 0);
    assert(i < cva(static_cast<int>(cur_ts)).extent(0));
    return cva(static_cast<int>(cur_ts))(i);
  }

  real &operator[](int i) noexcept {
    assert(i >= 0);
    assert(i < cva(static_cast<int>(cur_ts)).extent(0));
    return cva(static_cast<int>(cur_ts))(i);
  }

  void print() const noexcept;

  Mesh<_ctrl_vols> &cur_mesh() noexcept {
    return cva(static_cast<int>(cur_ts));
  }

 protected:
  void flux_integration(const MeshT &cur_ts, const MeshT &partial_ts,
                        MeshT &next_ts, const real stage_dt) noexcept;

  ND_Array<Mesh<_ctrl_vols>, time_stages> cva;

  // current time
  real t;

  TimeStage cur_ts;
  static constexpr TimeStage partial_ts = TimeStage::stage_partial;
};

#endif  // _WAVE_EQN_HPP_
