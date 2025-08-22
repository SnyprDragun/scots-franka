
/*
 * dcdc_n_link_arm.cc
 *
 *  created: Aug 2025 (adapted)
 *   author: adapted by ChatGPT
 *
 *  This file generalizes the original 2-link SCOTS example to an n-link
 *  articulated arm. It provides a framework that supports an arbitrary
 *  number of joints. Dynamics are implemented as a simple decoupled
 *  joint double-integrator (theta_dot = omega, omega_dot = (u - d*omega)/I).
 *
 *  NOTE: For a true coupled n-link rigid-body dynamics (Coriolis, gravity,
 *  mass matrix), replace the `joint_dynamics` implementation in `system_post`
 *  with your preferred model (e.g. recursive Newton-Euler or using a
 *  robotics library). The rest of the code (grid, abstraction, radius
 *  propagation) is written generally for arbitrary state/input dimensions.
 */

#include <limits>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

// === USER CONFIGURATION ===
// number of joints in the manipulator
static const int n_joints = 2; // change this to any positive integer

// derived dimensions
const int state_dim = 2 * n_joints; // [theta_0..theta_{n-1}, omega_0..omega_{n-1}]
const int input_dim = n_joints;     // torques at each joint

// integration time step
const double tau = 1.0;

using state_type = std::vector<double>;
using input_type = std::vector<double>;
using abs_type = scots::abs_type;

// Simple per-joint inertias and damping (users may modify)
static const double default_inertia = 1.0;
static const double default_damping = 0.1;

/* system_post integrates the continuous dynamics for time tau; prepared for
   arbitrary number of joints. The dynamics used here are intentionally simple
   (decoupled double integrators) so the example compiles and runs. Replace
   the body of `joint_dynamics` with the real rigid-body equations for true
   coupling between joints. */
auto system_post = [](state_type &x, input_type &u) -> void {

  // rhs: fills xx (time derivative) given state x and input u
  auto rhs = [](state_type &xx, const state_type &x, const input_type &u) -> void {
    // state layout: [theta_0..theta_{n-1}, omega_0..omega_{n-1}]
    int n = n_joints;
    // ensure sizes
    xx.assign(state_dim, 0.0);

    // theta_dot = omega
    for (int i = 0; i < n; ++i) {
      xx[i] = x[n + i];
    }

    // Simple decoupled dynamics: omega_dot = (u - d*omega)/I
    for (int i = 0; i < n; ++i) {
      double I = default_inertia;   // inertia of joint i (replace with vector if needed)
      double d = default_damping;   // damping
      double ui = (i < (int)u.size()) ? u[i] : 0.0;
      double omega = x[n + i];
      xx[n + i] = (ui - d * omega) / I;
    }
  };

  // call SCOTS fixed-step RK4 integrator (works with std::vector state)
  scots::runge_kutta_fixed4(rhs, x, u, state_dim, tau, 10);
};

/* radius_post computes a bound on the propagation of uncertainties
   (r) when applying input u at state x. For general n we compute a
   linearization (Jacobian) numerically with finite differences and then
   propagate r <- | r + tau * A * r | (conservative). */

auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
    const double eps = 1e-6; // finite difference epsilon
    int n = state_dim;

    // Safety: ensure sizes
    if ((int)r.size() != n) r.assign(n, 0.0);

    // Evaluate rhs at nominal (we re-use the dynamics from system_post but
    // extract it locally for numeric Jacobian). We'll implement the same
    // simple dynamics here so that radius_post is self-contained.
    auto dynamics = [&](const state_type &s, const input_type &in, state_type &out) {
        out.assign(n, 0.0);
        // theta_dot = omega
        for (int i = 0; i < n_joints; ++i) out[i] = s[n_joints + i];
        // omega_dot
        for (int i = 0; i < n_joints; ++i) {
            double I = default_inertia;
            double d = default_damping;
            double ui = (i < (int)in.size()) ? in[i] : 0.0;
            double omega = s[n_joints + i];
            out[n_joints + i] = (ui - d * omega) / I;
        }
    };

    // compute Jacobian A (n x n) with finite differences d f / d x
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    state_type f0;
    dynamics(x, u, f0);

    state_type x_pert = x;
    for (int j = 0; j < n; ++j) {
        x_pert = x;
        double h = eps * std::max(1.0, std::abs(x[j]));
        x_pert[j] += h;
        state_type f1;
        dynamics(x_pert, u, f1);
        for (int i = 0; i < n; ++i) {
            A(i, j) = (f1[i] - f0[i]) / h;
        }
    }

    // propagate radius conservatively: r_new = | r + tau * A * r |
    Eigen::VectorXd r_eig(n);
    for (int i = 0; i < n; ++i) r_eig(i) = r[i];

    Eigen::VectorXd result = r_eig + tau * A.cwiseAbs() * r_eig.cwiseAbs();

    for (int i = 0; i < n; ++i) r[i] = std::abs(result(i));
};

int main() {
  TicToc tt;

  // prepare box bounds and grid spacing for the state space
  state_type s_lb(state_dim), s_ub(state_dim), s_eta(state_dim);

  // Example bounds (angles between -pi and pi, velocities bounded)
  for (int i = 0; i < n_joints; ++i) {
    s_lb[i] = -M_PI;         // theta lower
    s_ub[i] =  M_PI;         // theta upper
    s_eta[i] = 0.1;          // angle grid resolution
  }
  for (int i = 0; i < n_joints; ++i) {
    s_lb[n_joints + i] = -5.0;  // omega lower
    s_ub[n_joints + i] =  5.0;  // omega upper
    s_eta[n_joints + i] = 0.2;  // velocity grid resolution
  }

  scots::UniformGrid ss(state_dim, s_lb, s_ub, s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();

  // input grid (torques)
  input_type i_lb(input_dim), i_ub(input_dim), i_eta(input_dim);
  for (int i = 0; i < input_dim; ++i) {
    i_lb[i] = -2.0;
    i_ub[i] =  2.0;
    i_eta[i] = 0.5;
  }

  scots::UniformGrid is(input_dim, i_lb, i_ub, i_eta);
  is.print_info();

  std::cout << "Computing the transition function: " << std::endl;
  scots::TransitionFunction tf;
  scots::Abstraction<state_type, input_type> abs(ss, is);

  tt.tic();
  abs.compute_gb(tf, system_post, radius_post);
  tt.toc();

  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

  // define a target set example: first joint angle near 0.6 (user can change)
  auto target = [&ss](const abs_type& idx) {
    state_type x(state_dim);
    ss.itox(idx, x);
    // first joint theta in [0.5, 0.7]
    int j = 0; // joint index to define target
    double eta = ss.get_eta()[j];
    double low = 0.5, high = 0.7;
    if ((low <= (x[j] - eta / 2.0) && (x[j] + eta / 2.0) <= high))
      return true;
    return false;
  };

  // write target to file
  write_to_file(ss, target, "target");

  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win = scots::solve_reachability_game(tf, target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if (write_to_file(scots::StaticController(ss, is, std::move(win)), "scots"))
    std::cout << "Done. \n";

  return 0;
}
