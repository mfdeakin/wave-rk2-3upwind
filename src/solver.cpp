
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
    solver.timestep_rk1_1st();
  }
  plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

template <int ctrl_vols>
void run_p1_rk2(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.2);
  while(solver.time() < 1.0) {
    solver.timestep_rk2_1st();
  }
  plot_err_p1(Plot, solver.cur_mesh(), solver.time());
}

template <int ctrl_vols>
void run_p1_4th(py::object Plot) {
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, 0.02);
  while(solver.time() < 1.0) {
    solver.timestep_4th();
  }
  plot_err_p1(Plot, solver.cur_mesh(), solver.time());
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

  run_p1_4th<20>(Plot);
  run_p1_4th<40>(Plot);
  run_p1_4th<80>(Plot);
  Title("4th Order Errors");
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

bool check_stability_rk1_1(py::object Plot, const real cfl) {
  constexpr int ctrl_vols = 500;
  WaveEqnSolver<ctrl_vols, Part2> solver(0.0, 5.0, cfl / Part2::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk1_1st();
  }

  std::stringstream label;
  label << "Solution for dt = " << solver.dt() << " (CFL = " << cfl << ")";
  plot_mesh(Plot, solver.cur_mesh(), label.str(), 0.0);
  return false;
}

bool check_stability_rk2_1(py::object Plot, const real cfl) {
  constexpr int ctrl_vols = 500;
  WaveEqnSolver<ctrl_vols, Part2> solver(0.0, 20.0, cfl / Part2::velocity());
  while(solver.time() < 8.0) {
    solver.timestep_rk2_1st();
  }

  std::stringstream label;
  label << "Solution for dt = " << solver.dt() << " (CFL = " << cfl << ")";
  plot_mesh(Plot, solver.cur_mesh(), label.str(), 12.0);
  return false;
}

bool check_stability_rk1_2(py::object Plot, const real cfl) {
  constexpr int ctrl_vols = 500;
  WaveEqnSolver<ctrl_vols, Part2> solver(0.0, 20.0, cfl / Part2::velocity());
  while(solver.time() < 8.0) {
    solver.timestep_rk1();
  }

  std::stringstream label;
  label << "Solution for dt = " << solver.dt() << " (CFL = " << cfl << ")";
  plot_mesh(Plot, solver.cur_mesh(), label.str(), 10.0);
  return false;
}

bool check_stability_rk2_2(py::object Plot, const real cfl) {
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

void part_2(py::object Plot, py::object Title, py::object Legend,
            py::object Show) {
  check_stability_rk1_1(Plot, 0.505);
  check_stability_rk1_1(Plot, 0.5);
  check_stability_rk1_1(Plot, 0.495);
  Legend();
  Title("1st Order Upwind with RK1");
  Show();

  check_stability_rk2_1(Plot, 0.505);
  check_stability_rk2_1(Plot, 0.5);
  check_stability_rk2_1(Plot, 0.495);
  Legend();
  Title("1st Order Upwind with RK2");
  Show();

  check_stability_rk1_2(Plot, 0.1625);
  check_stability_rk1_2(Plot, 0.15);
  check_stability_rk1_2(Plot, 0.125);
  check_stability_rk1_2(Plot, 0.1);
  Legend();
  Title("2nd Order Upwind with RK1");
  Show();

  check_stability_rk2_2(Plot, 0.505);
  check_stability_rk2_2(Plot, 0.5);
  check_stability_rk2_2(Plot, 0.495);
  Legend();
  Title("2nd Order Upwind with RK2");
  Show();
}

template <typename MeshT>
real l1_error(const MeshT &mesh, const real time) {
  real err = 0.0;
  for(int i = 1; i < mesh.extent(0); i++) {
    const real sol = Part1::avg_solution_val(mesh.x1(i), mesh.x2(i), time);
    err += std::abs(mesh(i) - sol);
  }
  return err / (mesh.extent(0) - 1);
}

template <typename MeshT>
real linf_error(const MeshT &mesh, const real time) {
  real err = 0.0;
  for(int i = 1; i < mesh.extent(0); i++) {
    const real sol = Part1::avg_solution_val(mesh.x1(i), mesh.x2(i), time);
    if(err < std::abs(mesh(i) - sol)) {
      err = std::abs(mesh(i) - sol);
    }
  }
  return err;
}

void plot_time_errs(py::object Plot, const real cfl, const std::string name,
                    std::function<void(WaveEqnSolver<80, Part1> &)> timestep) {
  constexpr int ctrl_vols = 80;
  std::vector<real> times;
  std::vector<real> l1_errs;
  std::vector<real> linf_errs;
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0, cfl / Part1::velocity());
  while(solver.time() < 1.0) {
    times.push_back(solver.time());
    l1_errs.push_back(l1_error(solver.cur_mesh(), solver.time()));
    linf_errs.push_back(linf_error(solver.cur_mesh(), solver.time()));
    timestep(solver);
  }
  Plot(py::array(times.size(), times.data()),
       py::array(l1_errs.size(), l1_errs.data()),
       "label"_a = name + " L1 Error");
  Plot(py::array(times.size(), times.data()),
       py::array(linf_errs.size(), linf_errs.data()),
       "label"_a = name + " Linf Error", "linestyle"_a = "dashed");
}

template <int ctrl_vols>
void compute_spatial_disc_rk1_1(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         1.0 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk1_1st();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk1_1<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk1_1<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

template <int ctrl_vols>
void compute_spatial_disc_rk1_2(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         0.5 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk1();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk1_2<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk1_2<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

template <int ctrl_vols>
void compute_spatial_disc_rk2_1(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         1.0 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk2_1st();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk2_1<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk2_1<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

template <int ctrl_vols>
void compute_spatial_disc_rk2_2(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         0.5 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk2();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk2_2<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk2_2<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

template <int ctrl_vols>
void compute_spatial_disc_rk4_1(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         1.0 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_rk4_1st();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk4_1<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk4_1<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

template <int ctrl_vols>
void compute_spatial_disc_rk4_2(std::vector<real> &vols,
                                std::vector<real> &l1_err,
                                std::vector<real> &linf_err) {
  vols.push_back(ctrl_vols);
  WaveEqnSolver<ctrl_vols, Part1> solver(0.0, 1.0,
                                         0.5 * 0.8 / Part1::velocity());
  while(solver.time() < 1.0) {
    solver.timestep_4th();
  }
  l1_err.push_back(l1_error(solver.cur_mesh(), solver.time()));
  linf_err.push_back(linf_error(solver.cur_mesh(), solver.time()));
  compute_spatial_disc_rk4_2<ctrl_vols * 2>(vols, l1_err, linf_err);
}

template <>
void compute_spatial_disc_rk4_2<640>(std::vector<real> &vols,
                                     std::vector<real> &l1_err,
                                     std::vector<real> &linf_err) {}

void plot_spatial_errs(
    py::object Plot, std::string label,
    std::function<void(std::vector<real> &, std::vector<real> &,
                       std::vector<real> &)>
        compute) {
  std::vector<real> vols;
  std::vector<real> l1_errs;
  std::vector<real> linf_errs;
  compute(vols, l1_errs, linf_errs);
  py::object LogLog = py::module::import("matplotlib.pyplot").attr("loglog");
  LogLog(py::array(vols.size(), vols.data()),
         py::array(l1_errs.size(), l1_errs.data()),
         "label"_a = label + " L1 Error");
  LogLog(py::array(vols.size(), vols.data()),
         py::array(linf_errs.size(), linf_errs.data()),
         "label"_a = label + " Linf Error", "linestyle"_a = "dashed");
}

void part_3(py::object Plot, py::object Title, py::object Legend,
            py::object Show) {
  plot_time_errs(
      Plot, 0.5 * 0.8, "1st order upwind RK1",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_rk1_1st(); });
  plot_time_errs(
      Plot, 0.5 * 0.8, "2nd order upwind RK1",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_rk1(); });
  plot_time_errs(
      Plot, 0.5 * 0.8, "1st order upwind RK2",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_rk2_1st(); });
  plot_time_errs(
      Plot, 0.5 * 0.8, "2nd order upwind RK2",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_rk2(); });
  plot_time_errs(
      Plot, 0.5 * 0.8, "1st order upwind 4th order",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_4th_1st(); });
  plot_time_errs(
      Plot, 0.5 * 0.8, "2nd order upwind 4th order",
      [](WaveEqnSolver<80, Part1> &solver) { solver.timestep_4th(); });
  Legend();
  Title("Errors of each scheme combination wrt. time");
  Show();

  plot_spatial_errs(Plot, "1st order upwind RK1",
                    compute_spatial_disc_rk1_1<10>);
  plot_spatial_errs(Plot, "2nd order upwind RK1",
                    compute_spatial_disc_rk1_2<10>);
  plot_spatial_errs(Plot, "1st order upwind RK2",
                    compute_spatial_disc_rk2_1<10>);
  plot_spatial_errs(Plot, "2nd order upwind RK2",
                    compute_spatial_disc_rk2_2<10>);
  // plot_spatial_errs(Plot, "1st order upwind RK4",
  //                   compute_spatial_disc_rk4_1<10>);
  plot_spatial_errs(Plot, "2nd order upwind RK4",
                    compute_spatial_disc_rk4_2<10>);
  Legend();
  Title("Errors of each scheme combination wrt. space convergence");
  Show();
}

int main(const int argc, const char **argv) {
  // Our Python instance
  py::scoped_interpreter _{};
  py::object Plot   = py::module::import("matplotlib.pyplot").attr("plot");
  py::object Title  = py::module::import("matplotlib.pyplot").attr("title");
  py::object Legend = py::module::import("matplotlib.pyplot").attr("legend");
  py::object Show   = py::module::import("matplotlib.pyplot").attr("show");
  part_1(Plot, Title, Legend, Show);
  part_2(Plot, Title, Legend, Show);
  part_3(Plot, Title, Legend, Show);
  return 0;
}
