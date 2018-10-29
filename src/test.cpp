
#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "wave_eqn.hpp"

using rngAlg = std::mt19937_64;

TEST(boundary_conds, t0) {
  constexpr int ctrl_vols = 160;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.0);
  for(int i = 0; i < ctrl_vols + 2; i++) {
    EXPECT_NEAR(solver[i], -std::sin(2 * pi * Solver::MeshT::median_x(i)),
                1e-4);
  }
}

TEST(boundary_conds, x0) {
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

TEST(flux_integral, poly) {
  constexpr int ctrl_vols = 4096;

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

TEST(solver, known_sol) {
  constexpr int ctrl_vols = 8192;

  using Solver = WaveEqnSolver<ctrl_vols, Part1>;

  Solver solver(0.25);

  for(int ts = 0; ts < 6; ts++) {
    EXPECT_EQ(solver.time(), solver.dt * ts);
    solver.timestep_rk2();
    const real t = solver.time();
    auto &mesh   = solver.cur_mesh();
    for(int i = 1; i < mesh.extent(0); i++) {
      const real x1 = mesh.x1(i);
      const real x2 = mesh.x2(i);
      // It's not clear why this isn't a better approximation of the
      // const real avg_val = (std::cos(2.0 * pi * (2.0 * t - x2)) -
      //                       std::cos(2.0 * pi * (2.0 * t - x1))) /
      //                      (2.0 * pi * mesh.dx);
      const real avg_val =
          (std::sin(pi * (x2 - x1)) * std::sin(pi * (4.0 * t - x1 - x2))) /
          (pi * (x2 - x1));

      // const real x = mesh.median_x(i);

      // const real median_val = std::sin(2.0 * pi * (2.0 * t - x));
      // EXPECT_NEAR(mesh(i), median_val, 2e-2);
      EXPECT_NEAR(mesh(i), avg_val, 1e-2);
    }
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
