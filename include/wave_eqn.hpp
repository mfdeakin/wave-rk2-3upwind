
#ifndef _WAVE_EQN_HPP_
#define _WAVE_EQN_HPP_

#include <assert.h>
#include <cmath>
#include <limits>

#include "nd_array/nd_array.hpp"

using real = double;

constexpr real pi = 3.1415926535897932;

class Part1 {
 public:
  static constexpr real velocity() noexcept { return 2.0; }

  static real solution_val(const real x, const real time) noexcept {
    return std::sin(2.0 * pi * (velocity() * time - x));
  }

  static real avg_solution_val(const real x1, const real x2,
                               const real time) noexcept {
    const real dx = x2 - x1;
    return (std::cos(2.0 * pi * (velocity() * time - x2)) -
            std::cos(2.0 * pi * (velocity() * time - x1))) /
           (2.0 * pi * dx);
  }

  static real deriv_val(const real x, const real time) noexcept {
    return 2.0 * pi * std::cos(2.0 * pi * (velocity() * time - x));
  }

  static real avg_deriv_val(const real x1, const real x2,
                            const real time) noexcept {
    return (solution_val(x2, time) - solution_val(x1, time)) / (x2 - x1);
  }

  static real boundary_val(const real time) noexcept {
    return solution_val(0.0, time);
  }

  static real initial_val(const real x) noexcept {
    return solution_val(x, 0.0);
  }

  template <typename MeshT>
  static void fill_ghostcells(MeshT &mesh, real time) noexcept {
    mesh(0) = 2.0 * boundary_val(time) - mesh(1);
  }

  template <typename MeshT>
  static void initial_conds(MeshT &mesh) noexcept {
    for(int i = 0; i < mesh.extent(0); i++) {
      mesh(i) = avg_solution_val(mesh.x1(i), mesh.x2(i), 0.0);
    }
  }
};

class Part2 {
 public:
  static constexpr real velocity() noexcept { return 2.0; }

  constexpr static real solution_val(const real x, const real time) noexcept {
    return std::numeric_limits<real>::quiet_NaN();
  }

  constexpr static real avg_solution_val(const real x1, const real x2,
                                         const real time) noexcept {
    return std::numeric_limits<real>::quiet_NaN();
  }

  static real boundary_val(const real time) noexcept {
    return Part1::boundary_val(time);
  }

  template <typename MeshT>
  static void fill_ghostcells(MeshT &mesh, real time) noexcept {
    mesh(0) = 2.0 * boundary_val(time) - mesh(1);
  }

  template <typename MeshT>
  static void initial_conds(MeshT &mesh) noexcept {
    for(int i = 0; i < mesh.extent(0); i++) {
      if(mesh.median_x(i) < 1.0) {
        mesh(i) = -mesh.median_x(i);
      } else {
        mesh(i) = 0.0;
      }
    }
  }
};

template <int _ctrl_vols, typename _boundary_conds>
class Mesh : public ND_Array<real, _ctrl_vols + 1> {
 public:
  using BoundaryConds = _boundary_conds;

  constexpr real x1(const int cell_x) const noexcept {
    return min_x() + (cell_x - 1) * dx();
  }

  constexpr real x2(const int cell_x) const noexcept {
    return min_x() + cell_x * dx();
  }

  constexpr real median_x(const int cell_x) const noexcept {
    return min_x() - dx() / 2 + cell_x * dx();
  }

  constexpr int cell_idx(const real x) const noexcept {
    return static_cast<int>((x - min_x()) / dx()) + 1;
  }

  constexpr Mesh(const real _min_x = 0.0, const real _max_x = 1.0) noexcept
      : _min_x(_min_x), _max_x(_max_x), _dx((_max_x - _min_x) / _ctrl_vols) {}

  // Computes the CVA of dT/dx with a first order approximation
  constexpr real flux_integral_1st(const int i, const real time) const
      noexcept {
    assert(i > 0);
    assert(i < this->extent(0));
    return ((*this)(i) - (*this)(i - 1)) / dx();
  }

  // Computes the CVA of dT/dx with a second order approximation
  constexpr real flux_integral_2nd(const int i, const real time) const
      noexcept {
    assert(i > 0);
    assert(i < this->extent(0));
    if(i > 1) {
      return (3.0 * (*this)(i)-4.0 * (*this)(i - 1) + (*this)(i - 2)) /
             (2.0 * dx());
    } else {
      // Evaluate the flux at the right side of the cell, then subtracting the
      // given value at the left side of the cell and dividing by dx gives the
      // approximate derivative
      return ((3.0 * (*this)(1) - 1.0 * (*this)(0)) / 2.0 -
              BoundaryConds::boundary_val(time)) /
             dx();
    }
  }

  // Computes the CVA of dT/dx with a third order approximation
  constexpr real flux_integral_3rd(const int i, const real time) const
      noexcept {
    assert(i > 0);
    assert(i < this->extent(0));
    if(i > 3) {
      return (11.0 * (*this)(i)-18.0 * (*this)(i - 1) + 9.0 * (*this)(i - 2) -
              2.0 * (*this)(i - 3)) /
             (6.0 * dx());
    } else {
      // const real T0 = -6.0 * BoundaryConds::boundary_val(time) +
      //                 5.0 * (*this)(i - 1) - (*this)(i);
      const real T0  = BoundaryConds::avg_solution_val(x1(0), x2(0), time);
      // const real Tm1 = BoundaryConds::avg_solution_val(x1(-1), x2(-1), time);
      const real T_5half = (11.0 * (*this)(2) - 7.0 * (*this)(1) + T0) / 6.0;
      // const real T_5half = BoundaryConds::solution_val(x1(3), time);
      if(i == 3) {
        const real T_7half =
            (11.0 * (*this)(3) - 7.0 * (*this)(2) + 2.0 * (*this)(1)) / 6.0;
        return (T_7half - T_5half) / (dx());
      } else {
        // const real Tm1 = 2.0 * (*this)(i - 1) + 5.0 * T0 -
        //                  6.0 * BoundaryConds::boundary_val(time);
        // const real T_3half = (11.0 * (*this)(1) - 7.0 * T0 + 2.0 * Tm1)
        // / 6.0;
        const real T_3half = BoundaryConds::solution_val(x1(2), time);
        if(i == 2) {
          // Evaluate the flux at the right side of the cell, then subtracting
          // the given value at the left side of the cell and dividing by dx
          // gives the approximate derivative
          return (T_5half - T_3half) / (dx());
        } else {
          return (T_3half - BoundaryConds::boundary_val(time)) / dx();
        }
      }
    }
  }

  constexpr real min_x() const noexcept { return _min_x; }

  constexpr real max_x() const noexcept { return _max_x; }

  constexpr real dx() const noexcept { return _dx; }

 protected:
  real _min_x;
  real _max_x;
  real _dx;
};

enum class TimeStage { stage_0 = 0, stage_1 = 1, stage_2 = 2 };

template <int _ctrl_vols, typename _boundary_conds>
class WaveEqnSolver {
 public:
  using BoundaryConds = _boundary_conds;
  using MeshT         = Mesh<_ctrl_vols, BoundaryConds>;

  static constexpr int time_stages = 3;
  static constexpr int ctrl_vols   = _ctrl_vols;

  constexpr WaveEqnSolver(const real min_x = 0.0, const real max_x = 1.0,
                          const real dtdx = 0.2) noexcept
      : mesh_1(min_x, max_x),
        mesh_2(min_x, max_x),
        mesh_3(min_x, max_x),
        _time(0.0),
        _dt(dtdx * mesh_1.dx()),
        _cur_ts(TimeStage::stage_0) {
    BoundaryConds::initial_conds(cur_mesh());
  }

  constexpr void fill_ghostcells() noexcept {
    BoundaryConds::fill_ghostcells(cur_mesh(), time());
  }

  constexpr void fill_ghostcells(TimeStage ts, const real t) noexcept {
    BoundaryConds::fill_ghostcells(mesh(ts), t);
  }

  void timestep_rk1_1st() noexcept;
  void timestep_rk2_1st() noexcept;
  void timestep_rk4_1st() noexcept;
  void timestep_3rd_1st() noexcept;
  void timestep_4th_1st() noexcept;

  void timestep_rk1() noexcept;
  void timestep_rk2() noexcept;
  void timestep_rk4() noexcept;
  void timestep_3rd() noexcept;
  void timestep_4th() noexcept;

  real time() const noexcept { return _time; }
  real dx() const noexcept { return mesh_1.dx(); }
  real dt() const noexcept { return _dt; }

  // Performs a linear interpolation between the cells to compute x
  real operator()(const real x) const noexcept;

  real &operator[](int i) noexcept { return cur_mesh()(i); }

  real operator[](int i) const noexcept { return cur_mesh()(i); }

  const MeshT &cur_mesh() const noexcept { return mesh(_cur_ts); }

  MeshT &cur_mesh() noexcept { return mesh(_cur_ts); }

 protected:
  constexpr const MeshT &mesh(TimeStage ts) const noexcept {
    switch(ts) {
      case TimeStage::stage_0:
        return mesh_1;
      case TimeStage::stage_1:
        return mesh_2;
      default:
        return mesh_3;
    }
  }

  constexpr MeshT &mesh(TimeStage ts) noexcept {
    switch(ts) {
      case TimeStage::stage_0:
        return mesh_1;
      case TimeStage::stage_1:
        return mesh_2;
      default:
        return mesh_3;
    }
  }

  void flux_integration_1st(const MeshT &cur_ts, const MeshT &partial_ts,
                            MeshT &next_ts, const real bc_time,
                            const real stage_dt) noexcept;

  void flux_integration(const MeshT &cur_ts, const MeshT &partial_ts,
                        MeshT &next_ts, const real bc_time,
                        const real stage_dt) noexcept;

  void flux_integration_3rd(const MeshT &cur_ts, const MeshT &partial_ts,
                            MeshT &next_ts, const real bc_time,
                            const real stage_dt) noexcept;

  MeshT mesh_1, mesh_2, mesh_3;

  // current time information
  real _time;
  const real _dt;

  TimeStage _cur_ts;
};

#endif  // _WAVE_EQN_HPP_
