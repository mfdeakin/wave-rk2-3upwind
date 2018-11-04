
#include "wave_eqn.hpp"

#include <string>
#include <vector>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace py::literals;

template <typename MeshT>
void plot_mesh(py::object Plot, const MeshT &mesh, std::string label,
               const real min_x = 0.0, const real max_x = 20.0) {
  std::vector<real> x_vals;
  std::vector<real> y_vals;
  for(int i = 1; i < mesh.extent(0); i++) {
    const real x = mesh.median_x(i);
    if(x > min_x) {
      x_vals.push_back(x);
      y_vals.push_back(mesh(i));
    }
  }

  Plot(py::array(x_vals.size(), x_vals.data()),
       py::array(y_vals.size(), y_vals.data()), "label"_a = label);
}

template <int ctrl_vols, typename _boundaries>
void plot_err_p1(py::object Plot, const Mesh<ctrl_vols, _boundaries> &src_mesh,
                 const real time) {
  Mesh<ctrl_vols, _boundaries> mesh;
  for(int i = 1; i < mesh.extent(0); i++) {
    mesh(i) = src_mesh(i) -
              _boundaries::avg_solution_val(mesh.x1(i), mesh.x2(i), time);
  }
  std::stringstream label_arg;
  label_arg << "Error for " << ctrl_vols << " CVs";
  plot_mesh(Plot, mesh, label_arg.str());
}

template <int ctrl_vols>
void run_p1_rk1(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.2);
  while(solver.time() < 1.0) {
    solver.timestep_rk1();
  }
  plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

template <int ctrl_vols>
void run_p1_rk2(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.2);
  while(solver.time() < 1.0) {
    solver.timestep_rk2();
  }
	plot_mesh(Plot, solver.cur_mesh(), "RK2 Solution");
  // plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

template <int ctrl_vols>
void run_p1_rk4(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.2);
  while(solver.time() < 1.0) {
    solver.timestep_rk4();
  }
	plot_mesh(Plot, solver.cur_mesh(), "RK4 Solution");
  // plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

void part_1(py::object Plot, py::object Title, py::object Legend,
            py::object Show) {
  run_p1_rk1<20>(Plot);
  run_p1_rk1<40>(Plot);
  run_p1_rk1<80>(Plot);
  Title("RK 1 Errors");
  Legend();
  Show();

  run_p1_rk2<20>(Plot);
  run_p1_rk2<40>(Plot);
  run_p1_rk2<80>(Plot);
  Title("RK 2 Errors");
  Legend();
  Show();

  run_p1_rk4<20>(Plot);
  run_p1_rk4<40>(Plot);
  run_p1_rk4<80>(Plot);
  Title("RK 4 Errors");
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
  WaveEqnSolver<ctrl_vols, Part2> solver(0.0, 20.0, cfl / Part2::velocity());
  while(solver.time() < 8.0) {
    solver.timestep_rk2();
  }

  std::stringstream label;
  label << "Solution for dt = " << solver.dt() << " (CFL = " << cfl << ")";
  plot_mesh(Plot, solver.cur_mesh(), label.str(), 12.0);
  return false;
}

void part_2(py::object Plot, py::object Legend, py::object Show) {
  check_stability(Plot, 0.505);
  check_stability(Plot, 0.5);
  check_stability(Plot, 0.495);
  Legend();
  Show();
}

int main(int argc, char **argv) {
  // Our Python instance
  py::scoped_interpreter _{};
  py::object Plot   = py::module::import("matplotlib.pyplot").attr("plot");
  py::object Title  = py::module::import("matplotlib.pyplot").attr("title");
  py::object Legend = py::module::import("matplotlib.pyplot").attr("legend");
  py::object Show   = py::module::import("matplotlib.pyplot").attr("show");
  part_1(Plot, Title, Legend, Show);
  part_2(Plot, Legend, Show);
}
