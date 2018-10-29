
#include "wave_eqn_internal.hpp"

template class Mesh<1024>;
template class WaveEqnSolver<1024, Part1>;
template class Mesh<2048>;
template class WaveEqnSolver<2048, Part1>;
template class Mesh<4096>;
template class WaveEqnSolver<4096, Part1>;
template class Mesh<8192>;
template class WaveEqnSolver<8192, Part1>;
template class Mesh<16384>;
template class WaveEqnSolver<16384, Part1>;
