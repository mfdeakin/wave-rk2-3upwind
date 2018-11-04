
#include "wave_eqn_internal.hpp"

// The following mesh sizes were used to determine when the l2 error goes
// beneath 1e-3
template class WaveEqnSolver<1024, Part1>;
template class WaveEqnSolver<1088, Part1>;
template class WaveEqnSolver<1104, Part1>;
template class WaveEqnSolver<1108, Part1>;
template class WaveEqnSolver<1110, Part1>;
template class WaveEqnSolver<1112, Part1>;
template class WaveEqnSolver<1120, Part1>;
template class WaveEqnSolver<1152, Part1>;
template class WaveEqnSolver<1280, Part1>;
template class WaveEqnSolver<1536, Part1>;
template class WaveEqnSolver<2048, Part1>;

template class WaveEqnSolver<1024, Part2>;
template class WaveEqnSolver<2048, Part2>;
