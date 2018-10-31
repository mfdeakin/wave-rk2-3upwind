
#include "wave_eqn_internal.hpp"

// The following 8 mesh sizes were used to determine when the l2 error goes
// beneath 1e-3
template class Mesh<512>;
template class WaveEqnSolver<512, Part1>;
template class Mesh<640>;
template class WaveEqnSolver<640, Part1>;
template class Mesh<704>;
template class WaveEqnSolver<704, Part1>;
// Error is less than 1e-3 for 712 cells
template class Mesh<712>;
template class WaveEqnSolver<712, Part1>;
template class Mesh<720>;
template class WaveEqnSolver<720, Part1>;
template class Mesh<736>;
template class WaveEqnSolver<736, Part1>;
template class Mesh<768>;
template class WaveEqnSolver<768, Part1>;
template class Mesh<1024>;
template class WaveEqnSolver<1024, Part1>;
template class Mesh<2048>;
template class WaveEqnSolver<2048, Part1>;
template class Mesh<4096>;
template class WaveEqnSolver<4096, Part1>;
