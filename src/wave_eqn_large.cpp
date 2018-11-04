
#include "wave_eqn_internal.hpp"

// The following mesh sizes were used to determine when the l2 error goes
// beneath 1e-4
template class WaveEqnSolver<4096, Part1>;
template class WaveEqnSolver<5120, Part1>;
template class WaveEqnSolver<5152, Part1>;
template class WaveEqnSolver<5160, Part1>;
template class WaveEqnSolver<5162, Part1>;
template class WaveEqnSolver<5164, Part1>;
template class WaveEqnSolver<5168, Part1>;
template class WaveEqnSolver<5184, Part1>;
template class WaveEqnSolver<5248, Part1>;
template class WaveEqnSolver<5376, Part1>;
template class WaveEqnSolver<5632, Part1>;
template class WaveEqnSolver<6144, Part1>;
template class WaveEqnSolver<8192, Part1>;
