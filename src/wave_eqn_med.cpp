
#include "wave_eqn_internal.hpp"

// The following mesh sizes were used to determine when the l2 error goes
// beneath 1e-3
template class Mesh<1024>;
template class WaveEqnSolver<1024, Part1>;
template class Mesh<1088>;
template class WaveEqnSolver<1088, Part1>;
template class Mesh<1104>;
template class WaveEqnSolver<1104, Part1>;
template class Mesh<1108>;
template class WaveEqnSolver<1108, Part1>;
template class Mesh<1110>;
template class WaveEqnSolver<1110, Part1>;
template class Mesh<1112>;
template class WaveEqnSolver<1112, Part1>;
template class Mesh<1120>;
template class WaveEqnSolver<1120, Part1>;
template class Mesh<1152>;
template class WaveEqnSolver<1152, Part1>;
template class Mesh<1280>;
template class WaveEqnSolver<1280, Part1>;
template class Mesh<1536>;
template class WaveEqnSolver<1536, Part1>;
template class Mesh<2048>;
template class WaveEqnSolver<2048, Part1>;
