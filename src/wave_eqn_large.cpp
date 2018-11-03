
#include "wave_eqn_internal.hpp"

// The following mesh sizes were used to determine when the l2 error goes
// beneath 1e-4
template class Mesh<4096>;
template class WaveEqnSolver<4096, Part1>;
template class Mesh<5120>;
template class WaveEqnSolver<5120, Part1>;
template class Mesh<5152>;
template class WaveEqnSolver<5152, Part1>;
template class Mesh<5160>;
template class WaveEqnSolver<5160, Part1>;
template class Mesh<5162>;
template class WaveEqnSolver<5162, Part1>;
template class Mesh<5164>;
template class WaveEqnSolver<5164, Part1>;
template class Mesh<5168>;
template class WaveEqnSolver<5168, Part1>;
template class Mesh<5184>;
template class WaveEqnSolver<5184, Part1>;
template class Mesh<5248>;
template class WaveEqnSolver<5248, Part1>;
template class Mesh<5376>;
template class WaveEqnSolver<5376, Part1>;
template class Mesh<5632>;
template class WaveEqnSolver<5632, Part1>;
template class Mesh<6144>;
template class WaveEqnSolver<6144, Part1>;
template class Mesh<8192>;
template class WaveEqnSolver<8192, Part1>;
