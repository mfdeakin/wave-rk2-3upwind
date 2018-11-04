
#ifndef _WAVE_EQN_INTERNAL_HPP_
#define _WAVE_EQN_INTERNAL_HPP_

#include "wave_eqn.hpp"

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk1_1st() noexcept {
  // Implements the RK2 timestep
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  fill_ghostcells();
  flux_integration_1st(mesh(_cur_ts), mesh(_cur_ts), mesh(next_ts), time(),
                       dt());

  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk2_1st() noexcept {
  // Implements the RK2 timestep
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts = TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 2.0);

  // Next compute the complete timestep
  fill_ghostcells(partial_ts, time() + dt() / 2.0);
  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(next_ts),
                   time() + dt() / 2.0, dt());
  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk4_1st() noexcept {
  // Implements the RK4 timestep
  // TODO: Understand the cause of the errors
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;

  Mesh<_ctrl_vols, _bc> stage_1;
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), stage_1, time(), dt() / 2.0);

  Mesh<_ctrl_vols, _bc> stage_2;
  real stage_t = time() + dt() / 2.0;
  BoundaryConds::fill_ghostcells(stage_1, stage_t);
  flux_integration(mesh(_cur_ts), stage_1, stage_2, stage_t, dt() / 2.0);

  Mesh<_ctrl_vols, _bc> stage_3;
  stage_t = time() + dt();
  BoundaryConds::fill_ghostcells(stage_2, stage_t);
  flux_integration(mesh(_cur_ts), stage_2, stage_3, stage_t, dt());

  for(int i = 1; i < stage_3.extent(0); i++) {
    stage_3(i) += stage_1(i) + 2.0 * stage_2(i) + stage_3(i);
  }

  stage_t = time() + dt() * 4.0 / 6.0;
  BoundaryConds::fill_ghostcells(stage_3, stage_t);
  flux_integration(mesh(_cur_ts), stage_3, mesh(next_ts), stage_t, dt() / 6.0);

  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_3rd_1st() noexcept {
  // Implements the RK3 timestep
  const TimeStage partial2_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts =
      (_cur_ts == TimeStage::stage_2) ? TimeStage::stage_1 : TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 3.0);

  // Second stage
  // compute w(2) based on w(1) and the current timestep
  fill_ghostcells(partial_ts, time() + dt() / 3.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 3.0, dt() / 2.0);

  // Third stage
  fill_ghostcells(partial2_ts, time() + dt() / 2.0);

  flux_integration(mesh(_cur_ts), mesh(partial2_ts), mesh(partial_ts),
                   time() + dt() / 2.0, dt());

  _cur_ts = partial_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_4th_1st() noexcept {
  // Implements the RK3 timestep
  const TimeStage partial2_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts =
      (_cur_ts == TimeStage::stage_2) ? TimeStage::stage_1 : TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 4.0);

  // Second stage
  // compute w(2) based on w(1) and the current timestep
  fill_ghostcells(partial_ts, time() + dt() / 4.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 4.0, dt() / 3.0);

  // Third stage
  fill_ghostcells(partial2_ts, time() + dt() / 3.0);

  flux_integration(mesh(_cur_ts), mesh(partial2_ts), mesh(partial_ts),
                   time() + dt() / 3.0, dt() / 2.0);

  // Fourth stage
  fill_ghostcells(partial_ts, time() + dt() / 2.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 2.0, dt());

  _cur_ts = partial2_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk1() noexcept {
  // Implements the RK2 timestep
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(next_ts), time(), dt());

  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk2() noexcept {
  // Implements the RK2 timestep
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts = TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 2.0);

  // Next compute the complete timestep
  fill_ghostcells(partial_ts, time() + dt() / 2.0);
  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(next_ts),
                   time() + dt() / 2.0, dt());
  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk4() noexcept {
  // Implements the RK4 timestep
  // TODO: Understand the cause of the errors
  const TimeStage next_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;

  Mesh<_ctrl_vols, _bc> stage_1;
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), stage_1, time(), dt() / 2.0);

  Mesh<_ctrl_vols, _bc> stage_2;
  real stage_t = time() + dt() / 2.0;
  BoundaryConds::fill_ghostcells(stage_1, stage_t);
  flux_integration(mesh(_cur_ts), stage_1, stage_2, stage_t, dt() / 2.0);

  Mesh<_ctrl_vols, _bc> stage_3;
  stage_t = time() + dt();
  BoundaryConds::fill_ghostcells(stage_2, stage_t);
  flux_integration(mesh(_cur_ts), stage_2, stage_3, stage_t, dt());

  for(int i = 1; i < stage_3.extent(0); i++) {
    stage_3(i) += stage_1(i) + 2.0 * stage_2(i) + stage_3(i);
  }

  stage_t = time() + dt() * 4.0 / 6.0;
  BoundaryConds::fill_ghostcells(stage_3, stage_t);
  flux_integration(mesh(_cur_ts), stage_3, mesh(next_ts), stage_t, dt() / 6.0);

  _cur_ts = next_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_3rd() noexcept {
  // Implements the RK3 timestep
  const TimeStage partial2_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts =
      (_cur_ts == TimeStage::stage_2) ? TimeStage::stage_1 : TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 3.0);

  // Second stage
  // compute w(2) based on w(1) and the current timestep
  fill_ghostcells(partial_ts, time() + dt() / 3.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 3.0, dt() / 2.0);

  // Third stage
  fill_ghostcells(partial2_ts, time() + dt() / 2.0);

  flux_integration(mesh(_cur_ts), mesh(partial2_ts), mesh(partial_ts),
                   time() + dt() / 2.0, dt());

  _cur_ts = partial_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_4th() noexcept {
  // Implements the RK3 timestep
  const TimeStage partial2_ts =
      (_cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  const TimeStage partial_ts =
      (_cur_ts == TimeStage::stage_2) ? TimeStage::stage_1 : TimeStage::stage_2;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(mesh(_cur_ts), mesh(_cur_ts), mesh(partial_ts), time(),
                   dt() / 4.0);

  // Second stage
  // compute w(2) based on w(1) and the current timestep
  fill_ghostcells(partial_ts, time() + dt() / 4.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 4.0, dt() / 3.0);

  // Third stage
  fill_ghostcells(partial2_ts, time() + dt() / 3.0);

  flux_integration(mesh(_cur_ts), mesh(partial2_ts), mesh(partial_ts),
                   time() + dt() / 3.0, dt() / 2.0);

  // Fourth stage
  fill_ghostcells(partial_ts, time() + dt() / 2.0);

  flux_integration(mesh(_cur_ts), mesh(partial_ts), mesh(partial2_ts),
                   time() + dt() / 2.0, dt());

  _cur_ts = partial2_ts;

  _time += dt();
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::flux_integration_1st(
    const MeshT &cur_ts, const MeshT &partial_ts, MeshT &next_ts,
    const real bc_time, const real stage_dt) noexcept {
  assert(&next_ts != &cur_ts);
  assert(&next_ts != &partial_ts);
  for(int i = 1; i < next_ts.extent(0); i++) {
    // Compute the dT/dx term with the second order upwinding scheme
    // Then add the previous time stage to this times the stage_dt
    next_ts(i) =
        -2.0 * partial_ts.flux_integral_1st(i, bc_time) * stage_dt + cur_ts(i);
  }
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::flux_integration(
    const MeshT &cur_ts, const MeshT &partial_ts, MeshT &next_ts,
    const real bc_time, const real stage_dt) noexcept {
  assert(&next_ts != &cur_ts);
  assert(&next_ts != &partial_ts);
  for(int i = 1; i < next_ts.extent(0); i++) {
    // Compute the dT/dx term with the second order upwinding scheme
    // Then add the previous time stage to this times the stage_dt
    next_ts(i) =
        -2.0 * partial_ts.flux_integral_2nd(i, bc_time) * stage_dt + cur_ts(i);
  }
}

template <int _ctrl_vols, typename _bc>
real WaveEqnSolver<_ctrl_vols, _bc>::operator()(const real x) const noexcept {
  assert(cur_mesh().max_x() >= x);
  assert(cur_mesh().min_x() <= x);
  const real right_cell = x + cur_mesh().dx() / 2.0;
  const int i           = cur_mesh().cell_idx(right_cell);
  const real weight     = (x - cur_mesh().x1(i)) / cur_mesh().dx() * 2.0;
  return weight * cur_mesh()(i) + (1.0 - weight) * cur_mesh()(i - 1);
}

#endif  // _WAVE_EQN_INTERNAL_HPP_
