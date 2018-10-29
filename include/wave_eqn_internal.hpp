
#ifndef _WAVE_EQN_INTERNAL_HPP_
#define _WAVE_EQN_INTERNAL_HPP_

#include "wave_eqn.hpp"

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk2() noexcept {
	// Implements the RK2 timestep
  const TimeStage next_ts =
      (cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(cva(static_cast<int>(cur_ts)), cva(static_cast<int>(cur_ts)),
                   cva(static_cast<int>(partial_ts)), dt / 2.0);

  // Next compute the completed timestep
  _bc::fill_ghostcells(cva(static_cast<int>(partial_ts)), time() + dt / 2.0);

  flux_integration(cva(static_cast<int>(cur_ts)),
                   cva(static_cast<int>(partial_ts)),
                   cva(static_cast<int>(next_ts)), dt);
  cur_ts = next_ts;

  t += dt;
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::timestep_rk3() noexcept {
	// Implements the RK3 timestep
  const TimeStage partial2_ts =
      (cur_ts == TimeStage::stage_0) ? TimeStage::stage_1 : TimeStage::stage_0;
  // First compute w(1) as the partial timestep term
  // w(1) only relies on the current timestep,
  // so treat it as the partial timestep as well as the current one
  fill_ghostcells();
  flux_integration(cva(static_cast<int>(cur_ts)), cva(static_cast<int>(cur_ts)),
                   cva(static_cast<int>(partial_ts)), dt / 3.0);

  // Second stage
  _bc::fill_ghostcells(cva(static_cast<int>(partial_ts)), time() + dt / 3.0);

  flux_integration(cva(static_cast<int>(cur_ts)),
                   cva(static_cast<int>(partial_ts)),
                   cva(static_cast<int>(partial2_ts)), dt / 2.0);

	// Third stage
  _bc::fill_ghostcells(cva(static_cast<int>(partial2_ts)), time() + 2.0 * dt / 3.0);

  flux_integration(cva(static_cast<int>(cur_ts)),
                   cva(static_cast<int>(partial2_ts)),
                   cva(static_cast<int>(partial_ts)), dt);

	cva(static_cast<int>(cur_ts)) = cva(static_cast<int>(partial_ts));

  t += dt;
}

template <int _ctrl_vols, typename _bc>
void WaveEqnSolver<_ctrl_vols, _bc>::flux_integration(
    const MeshT &cur_ts, const MeshT &partial_ts, MeshT &next_ts,
    const real stage_dt) noexcept {
  assert(&next_ts != &cur_ts);
  assert(&next_ts != &partial_ts);
  for(int i = 1; i < next_ts.extent(0) - 1; i++) {
    // Compute the dT/dx term with the second order upwinding scheme
    // Then add the previous time stage to this times the stage_dt
    next_ts(i) = -2.0 * partial_ts.flux_integral(i) * stage_dt + cur_ts(i);
  }
}

template <int _ctrl_vols, typename _bc>
real WaveEqnSolver<_ctrl_vols, _bc>::operator()(const real x) const noexcept {
  assert(max_x >= x);
  assert(min_x <= x);
  const real right_cell = x + cva(0).dx / 2.0;
  const int i           = cva(0).cell_idx(right_cell);
  const real weight     = (x - cva(0).x1(i)) / cva(0).dx * 2.0;
  return weight * cva(static_cast<int>(cur_ts))(i) +
         (1.0 - weight) * cva(static_cast<int>(cur_ts))(i - 1);
}

#endif // _WAVE_EQN_INTERNAL_HPP_
