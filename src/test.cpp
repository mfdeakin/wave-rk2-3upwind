
#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "wave_eqn.hpp"

using rngAlg = std::mt19937_64;

TEST(boundary_conds, t0) {
  constexpr int ctrl_vols = 160;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.0);
  for(int i = 0; i < ctrl_vols + 1; i++) {
    EXPECT_NEAR(solver[i], -std::sin(2 * pi * solver.cur_mesh().median_x(i)),
                1e-4);
  }
}

// TEST(boundary_conds, x0_tstep) {
//   constexpr int ctrl_vols = 160;

//   using Solver = WaveEqnSolver<ctrl_vols, Part1>;

//   Solver solver(0.125);

//   for(int i = 0; i < 50; i++) {
//     solver.fill_ghostcells();
//     printf("sin(4 pi t): % .6f; boundary_val: % .6f; delta: % .6e\n",
//            std::sin(4.0 * pi * solver.time()),
//            Part1::boundary_val(solver.time()),
//            std::sin(4.0 * pi * solver.time()) -
//                Part1::boundary_val(solver.time()));
//     printf("mesh(0): % .6e; mesh(1): % .6e; avg: % .6e, diff: % .6e\n",
//            solver[0], solver[1], (solver[0] + solver[1]) / 2.0,
//            (solver[0] + solver[1]) / 2.0 - std::sin(4.0 * pi *
//            solver.time()));
//     // EXPECT_NEAR((mesh(0) + mesh(1)) / 2.0,
//     //             std::sin(4.0 * pi * solver.time()), 2e-12);
//     solver.timestep_rk2();
//   }
// }

TEST(boundary_conds, x0_rand) {
  constexpr int ctrl_vols = 160;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.125);

  std::random_device rd;
  rngAlg engine(rd());
  std::uniform_real_distribution<real> pdf(-1.0, 1.0);

  for(int i = 0; i < 100; i++) {
    solver[1] = pdf(engine);
    solver.fill_ghostcells();
    EXPECT_NEAR((solver[0] + solver[1]) / 2.0,
                std::sin(4.0 * pi * solver.time()), 1e-15);
  }
}

TEST(flux_integral, sinusoid) {
  constexpr int ctrl_vols = 1024;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.0, 1.0, 0.5);

  for(int ts = 0; ts < 10; ts++) {
    auto &mesh = solver.cur_mesh();
    printf("Time: % .3e; dt: % .3e; cva.dx: % .3e\n", solver.time(),
           solver.dt(), mesh.dx());

    for(int i = 1; i <= 10; i++) {
      const real x = mesh.median_x(i);
      // EXPECT_NEAR(mesh.flux_integral(i, solver.time()),
      //             -2.0 * pi * std::cos(2.0 * pi * (2.0 * solver.time() - x)),
      //             5e-3);
      printf(
          "i: %2d, x: % .3e, prev: % .6e, val: % .6e, fi: % .6e, deriv: % .6e, "
          "e: % .3e\n",
          i, x, mesh(i - 1), mesh(i), mesh.flux_integral(i, solver.time()),
          Part1::avg_deriv_val(mesh.x1(i), mesh.x2(i), solver.time()),
          mesh.flux_integral(i, solver.time()) -
              Part1::avg_deriv_val(mesh.x1(i), mesh.x2(i), solver.time()));
    }
    printf("\n");

    // We need to increment the time; this is the only way to do so
    solver.timestep_rk1();

    // Then overwrite the values with the exact solution
    // for(int i = 1; i <= ctrl_vols; i++) {
    // 	solver[i] = Part1::avg_solution_val(mesh.x1(i), mesh.x2(i),
    // solver.time());
    // }
    solver.fill_ghostcells();
  }
}

TEST(flux_integral, poly) {
  constexpr int ctrl_vols = 1024;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver;
  auto mesh = solver.cur_mesh();

  for(int i = 0; i < mesh.extent(0); i++) {
    const real x = mesh.median_x(i);

    mesh(i) = x * x * x + 0.25 * std::sin(2.0 * pi * x) + 3.0;
  }

  for(int i = 2; i < mesh.extent(0) - 1; i++) {
    const real x = mesh.median_x(i);
    EXPECT_NEAR(mesh.flux_integral(i, std::numeric_limits<real>::quiet_NaN()),
                3.0 * x * x + 0.5 * pi * std::cos(2.0 * pi * x), 5e-3);
  }
}

TEST(flux_integral, exp) {
  constexpr int ctrl_vols = 4096;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver;
  auto mesh = solver.cur_mesh();

  // sqrt is nan for negative x, so skip the ghost cell
  for(int i = 1; i < mesh.extent(0); i++) {
    const real x = mesh.median_x(i);

    mesh(i) = std::exp(std::sqrt(x * x * x));
  }

  // Avoid the boundary conditions by starting i > 2
  for(int i = 3; i < mesh.extent(0) - 1; i++) {
    const real x = mesh.median_x(i);
    EXPECT_NEAR(mesh.flux_integral(i, std::numeric_limits<real>::quiet_NaN()),
                std::exp(std::sqrt(x * x * x)) * 3.0 / 2.0 * std::sqrt(x),
                1e-3);
  }
}

TEST(solver, rk1_avg_sol) {
  constexpr int ctrl_vols = 1024;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.0, 1.0, 0.25);

  for(int ts = 0; ts < 5; ts++) {
    solver.timestep_rk1();
    const real t = solver.time();
    auto &mesh   = solver.cur_mesh();
    for(int i = 1; i < mesh.extent(0); i++) {
      const real x1 = mesh.x1(i);
      const real x2 = mesh.x2(i);

      const real avg_val = Part1::avg_solution_val(x1, x2, t);

      EXPECT_NEAR(mesh(i), avg_val, 1e-4);
    }
  }
}

// TEST(solver, rk2_avg_sol) {
//   constexpr int ctrl_vols = 8192;

//   using Solver = WaveEqnSolver<ctrl_vols, Part1>;

//   Solver solver(0.2);

//   for(int ts = 0; ts < 25; ts++) {
//     solver.timestep_rk2();
//     const real t = solver.time();
//     auto &mesh   = solver.cur_mesh();
//     for(int i = 1; i < mesh.extent(0); i++) {
//       const real x1 = mesh.x1(i);
//       const real x2 = mesh.x2(i);

//       const real avg_val = Part1::avg_solution_val(x1, x2, t);

//       EXPECT_NEAR(mesh(i), avg_val, 1e-2);
//     }
//   }
// }

// template <typename MeshT>
// real l2_error(const MeshT &mesh, const real time) {
//   real l2_err = 0.0;
//   for(int i = 1; i < mesh.extent(0) - 1; i++) {
//     // const real x1 = mesh.x1(i);
//     // const real x2 = mesh.x2(i);

//     // const real avg_val =
//     //     (std::sin(pi * (x2 - x1)) * std::sin(pi * (4.0 * time - x1 - x2)))
//     /
//     //     (pi * (x2 - x1));
//     const real x   = mesh.median_x(i);
//     const real val = std::sin(2.0 * pi * (2.0 * time - x));

//     const real err = mesh(i) - val;
//     l2_err += err * err;
//   }
//   return std::sqrt(l2_err);
// }

// template <typename MeshT>
// void print_error(const MeshT &mesh, const real time) {
//   for(int i = 1; i < mesh.extent(0) - 1; i++) {
//     printf("% .6f,  ", mesh.median_x(i));
//   }
//   printf("\n");
//   for(int i = 1; i < mesh.extent(0) - 1; i++) {
//     const real x   = mesh.median_x(i);
//     const real val = std::sin(2.0 * pi * (2.0 * time - x));
//     printf("% .6f,  ", mesh(i) - val);
//   }
//   printf("\n");
// }

// TEST(solver, rk2_l2_err_1en3) {
//   constexpr int ctrl_vols_below = 20;
//   constexpr int ctrl_vols_above = 40;

//   WaveEqnSolver<ctrl_vols_below, Part1> solver_below(0.2);
//   WaveEqnSolver<ctrl_vols_above, Part1> solver_above(0.2);

//   while(solver_above.time() < 1.0) {
//     solver_below.timestep_rk2();
//     solver_above.timestep_rk2();
//   }
//   printf("20 element mesh\n");
//   print_error(solver_below.cur_mesh(), solver_below.time());
//   printf("40 element mesh\n");
//   print_error(solver_above.cur_mesh(), solver_above.time());
//   const auto &mesh_below = solver_below.cur_mesh();
//   EXPECT_GT(l2_error(mesh_below, solver_below.time()), 1e-3);

//   const auto &mesh_above = solver_above.cur_mesh();
//   EXPECT_LT(l2_error(mesh_above, solver_above.time()), 1e-3);
// }

// TEST(solver, rk2_l2_err_1en4) {
//   constexpr int ctrl_vols_below = 7816;
//   constexpr int ctrl_vols_above = 7824;

//   WaveEqnSolver<ctrl_vols_below, Part1> solver_below(0.2);
//   WaveEqnSolver<ctrl_vols_above, Part1> solver_above(0.2);

//   while(solver_above.time() < 1.0) {
//     solver_below.timestep_rk2();
//     solver_above.timestep_rk2();
//   }
//   const auto &mesh_below = solver_below.cur_mesh();
//   EXPECT_GT(l2_error(mesh_below, solver_below.time()), 1e-4);

//   const auto &mesh_above = solver_above.cur_mesh();
//   EXPECT_LT(l2_error(mesh_above, solver_above.time()), 1e-4);
// }

// TEST(solver, rk3_avg_sol) {
//   constexpr int ctrl_vols = 8192;

//   using Solver = WaveEqnSolver<ctrl_vols, Part1>;

//   Solver solver(0.2);

//   for(int ts = 0; ts < 25; ts++) {
//     solver.timestep_rk3();
//     const real t = solver.time();
//     auto &mesh   = solver.cur_mesh();
//     for(int i = 1; i < mesh.extent(0); i++) {
//       const real x1 = mesh.x1(i);
//       const real x2 = mesh.x2(i);

//       const real avg_val =
//           (std::sin(pi * (x2 - x1)) * std::sin(pi * (4.0 * t - x1 - x2))) /
//           (pi * (x2 - x1));

//       EXPECT_NEAR(mesh(i), avg_val, 1e-2);
//     }
//   }
// }

// std::pair<real, real> richardson(const real fine, const real medium,
//                                  const real coarse) {
//   real order;
//   real extrap;

//   const double order_r = (fine - medium) / (medium - coarse);
//   order                = -std::log2(order_r);

//   extrap = fine - (medium - fine) * (1 / (std::pow(2, order) - 1));
//   return {order, extrap};
// }

// // Run the simulation until the maximum specified time
// // Return the average value of the solution
// template <int ctrl_vols>
// real run_simulation_rk2_p1(const real cfl, const real max_t) {
//   WaveEqnSolver<ctrl_vols, Part1> solver(cfl);

//   while(solver.time() < max_t) {
//     solver.timestep_rk2();
//   }

//   real total_val = 0.0;
//   for(int i = 1; i <= ctrl_vols; i++) {
//     total_val += solver[i];
//   }
//   return total_val * solver.dx();
// }

// TEST(solver, mesh_convergence_order_rk2) {
//   constexpr int coarse_n     = 20;
//   constexpr int med_n        = coarse_n * 2;
//   constexpr int fine_n       = med_n * 2;
//   const real coarse_sol      = run_simulation_rk2_p1<coarse_n>(0.2, 1.0);
//   const real med_sol         = run_simulation_rk2_p1<med_n>(0.2, 1.0);
//   const real fine_sol        = run_simulation_rk2_p1<fine_n>(0.2, 1.0);
//   const auto [order, extrap] = richardson(fine_sol, med_sol, coarse_sol);

//   EXPECT_GT(order, 2.0);
// }

// TEST(solver, time_convergence_order_rk2) {
//   constexpr int ctrl_vols  = 4096;
//   constexpr real coarse_dt = 0.125;
//   const real coarse_sol    =
//   run_simulation_rk2_p1<ctrl_vols>(coarse_dt, 1.0); const real med_sol  =
//   run_simulation_rk2_p1<ctrl_vols>(coarse_dt / 2.0, 1.0); const real fine_sol
//   = run_simulation_rk2_p1<ctrl_vols>(coarse_dt / 4.0, 1.0); const auto
//   [order, extrap] = richardson(fine_sol, med_sol, coarse_sol);

//   EXPECT_GT(order, 2.0);
// }

// // Run the simulation until the maximum specified time
// // Return the average value of the solution
// template <int ctrl_vols>
// real run_simulation_rk3_p1(const real cfl, const real max_t) {
//   WaveEqnSolver<ctrl_vols, Part1> solver(cfl);

//   while(solver.time() < max_t) {
//     solver.timestep_rk3();
//   }

//   real total_val = 0.0;
//   for(int i = 1; i <= ctrl_vols; i++) {
//     total_val += solver[i];
//   }
//   return total_val * solver.dx();
// }

// TEST(solver, mesh_convergence_order_rk3) {
//   constexpr int coarse_n     = 1024;
//   constexpr int med_n        = coarse_n * 2;
//   constexpr int fine_n       = med_n * 2;
//   const real coarse_sol      = run_simulation_rk3_p1<coarse_n>(0.2, 1.0);
//   const real med_sol         = run_simulation_rk3_p1<med_n>(0.2, 1.0);
//   const real fine_sol        = run_simulation_rk3_p1<fine_n>(0.2, 1.0);
//   const auto [order, extrap] = richardson(fine_sol, med_sol, coarse_sol);

//   EXPECT_NEAR(order, 2.0, 1e-1);
// }

// TEST(solver, time_convergence_order_rk3) {
//   constexpr int ctrl_vols  = 4096;
//   constexpr real coarse_dt = 0.125;
//   const real coarse_sol    =
//   run_simulation_rk3_p1<ctrl_vols>(coarse_dt, 1.0); const real med_sol  =
//   run_simulation_rk3_p1<ctrl_vols>(coarse_dt / 2.0, 1.0); const real fine_sol
//   = run_simulation_rk3_p1<ctrl_vols>(coarse_dt / 4.0, 1.0); const auto
//   [order, extrap] = richardson(fine_sol, med_sol, coarse_sol);

//   EXPECT_GT(order, 2.0);
// }

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
