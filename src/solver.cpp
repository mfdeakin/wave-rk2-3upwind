
#include "wave_eqn.hpp"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace py::literals;

template <int ctrl_vols, typename _boundaries>
void plot_err_p1(py::object Plot, const Mesh<ctrl_vols, _boundaries> &mesh,
                 const real time) {
  ND_Array<real, ctrl_vols> x_vals;
  for(int i = 1; i < mesh.extent(0); i++) {
    x_vals(i - 1) = mesh.median_x(i);
  }
  ND_Array<real, ctrl_vols> errs;
  for(int i = 1; i < mesh.extent(0); i++) {
    errs(i - 1) =
        mesh(i) - _boundaries::avg_solution_val(mesh.x1(i), mesh.x2(i), time);
  }

  std::stringstream label_arg;
  label_arg << "Error for " << ctrl_vols << " CVs";

  Plot(py::array(ctrl_vols, reinterpret_cast<real *>(&x_vals)),
       py::array(ctrl_vols, reinterpret_cast<real *>(&errs)),
       "label"_a = label_arg.str());
}

template <typename MeshT>
void plot_solution(py::object Plot, const MeshT &mesh, std::string label) {
  ND_Array<real, MeshT::extent(0) - 1> x_vals;
  for(int i = 1; i < mesh.extent(0); i++) {
    x_vals(i - 1) = mesh.median_x(i);
  }
  ND_Array<real, MeshT::extent(0) - 1> y_vals;
  for(int i = 1; i < mesh.extent(0); i++) {
    y_vals(i - 1) = mesh(i);
  }

  Plot(py::array(MeshT::extent(0) - 1, reinterpret_cast<real *>(&x_vals)),
       py::array(MeshT::extent(0) - 1, reinterpret_cast<real *>(&y_vals)),
       "label"_a = label);
}

template <int ctrl_vols>
void run_p1(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.2);
  while(solver.time() < 1.0) {
    solver.timestep_rk2();
  }
  plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

void part_1(py::object Plot, py::object Legend, py::object Show) {
  run_p1<20>(Plot);
  run_p1<40>(Plot);
  run_p1<80>(Plot);
  Legend();
  Show();
}

template <typename Mesh>
real integrate_abs(const real x1, const real x2, const Mesh &mesh) {
  const int start = mesh.cell_idx(x1);
  const int end   = mesh.cell_idx(x2);
  assert(end) < mesh.extent(0);
  real integral = 0.0;
  for(int i = start; i < end; i++) {
    integral += std::abs(mesh(i));
  }
  return integral;
}

bool check_stability(py::object Plot, const real cfl) {
  constexpr int ctrl_vols = 500;
  WaveEqnSolver<ctrl_vols, Part2> solver(0.0, 20.0, cfl);
  while(solver.time() < 8.0) {
    solver.timestep_rk2();
  }

	std::stringstream label;
	label << "Solution for CFL = " << cfl;
  plot_solution(Plot, solver.cur_mesh(), label.str());
  return false;
}

void part_2(py::object Plot, py::object Legend, py::object Show) {
  check_stability(Plot, 0.501);
  check_stability(Plot, 0.5);
  check_stability(Plot, 0.4);
	Legend();
	Show();
}

int main(int argc, char **argv) {
  // Our Python instance
  py::scoped_interpreter _{};
  py::object Plot   = py::module::import("matplotlib.pyplot").attr("plot");
  py::object Legend = py::module::import("matplotlib.pyplot").attr("legend");
  py::object Show   = py::module::import("matplotlib.pyplot").attr("show");
  part_1(Plot, Legend, Show);
  part_2(Plot, Legend, Show);
}
