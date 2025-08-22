// nlink_proper.cc
//
// SCOTS-ready C++ program for an n-link planar revolute manipulator
// Uses CRBA (mass matrix via RNEA columns) + RNEA for bias term (h = C*qd + g).
// Supplies system_post and radius_post for SCOTS `compute_gb`.
// Produces controller table CSV and (optionally) a closed-loop trajectory CSV.
//
// Compile with:
//   g++ -O3 -std=c++17 nlink_proper.cc -I/path/to/eigen -o nlink_proper
//
// Run with (example):
//   N_LINKS=3 ./nlink_proper
//

#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sys/resource.h>

#include <Eigen/Dense>

#include "scots.hh"
#include "TicToc.hh"
#include "RungeKutta4.hh"

#ifndef N_LINKS
#define N_LINKS 2
#endif

using namespace std;
using namespace scots;
using namespace Eigen;

constexpr size_t NL = static_cast<size_t>(N_LINKS);
constexpr size_t STATE_DIM = 2 * NL;
constexpr size_t INPUT_DIM = NL;
constexpr double tau = 1.0; // SCOTS sampling time for integration in abstraction

using state_type = array<double, STATE_DIM>;
using input_type = array<double, INPUT_DIM>;

/* ----------------- Model parameters ----------------- */
struct ModelParams {
  VectorXd l;   // link lengths (n)
  VectorXd lc;  // COM offsets (n)
  VectorXd m;   // masses (n)
  VectorXd I;   // planar inertias about COM (n)
  double g = 9.81;
  double viscous = 0.05;
  ModelParams() {
    l  = VectorXd::Constant(NL, 1.0);
    lc = 0.5 * l;
    m  = VectorXd::Constant(NL, 1.0);
    I  = VectorXd::Constant(NL, 0.05);
  }
};

/* ----------------- Small helpers ----------------- */
static inline Matrix2d R(double a) {
  double c = cos(a), s = sin(a);
  Matrix2d M; M << c, -s, s, c; return M;
}
static inline Vector2d R90(const Vector2d &v) {
  return Vector2d(-v.y(), v.x());
}

/* ----------------- Kinematics cache ----------------- */
struct KinCache {
  vector<Vector2d> P; // joint positions 0..n
  vector<Vector2d> C; // COM positions 0..n-1
  VectorXd a;         // cumulative angles
};

static KinCache forward_positions(const VectorXd &q, const ModelParams &mp) {
  KinCache kc;
  kc.P.resize(NL + 1);
  kc.C.resize(NL);
  kc.a.setZero(NL);

  kc.P[0] = Vector2d::Zero();
  double a_sum = 0.0;
  for (size_t i = 0; i < NL; ++i) {
    a_sum += q(static_cast<int>(i));
    kc.a(static_cast<int>(i)) = a_sum;
    Matrix2d Ri = R(a_sum);
    Vector2d r_link = Ri * Vector2d(mp.l(static_cast<int>(i)), 0.0);
    Vector2d r_com  = Ri * Vector2d(mp.lc(static_cast<int>(i)), 0.0);
    kc.C[i]   = kc.P[i] + r_com;
    kc.P[i+1] = kc.P[i] + r_link;
  }
  return kc;
}

/* ----------------- RNEA (planar specialization) -----------------0
   Given q, qdot, qdd compute tau = M(q)*qdd + h(q,qdot)
   h includes gravity, Coriolis/centrifugal, and optional viscous damping.
------------------------------------------------------------------ */
static VectorXd RNEA(const VectorXd &q,
                            const VectorXd &qdot,
                            const VectorXd &qdd,
                            const ModelParams &mp) {
  const size_t n = NL;
  VectorXd tau = VectorXd::Zero(n);

  KinCache kc = forward_positions(q, mp);

  // angular velocities and accelerations (scalars about z)
  VectorXd omega(n), alpha(n);
  omega.setZero(); alpha.setZero();
  double sum_w = 0.0, sum_a = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum_w += qdot(static_cast<int>(i));
    sum_a += qdd(static_cast<int>(i));
    omega(static_cast<int>(i)) = sum_w;
    alpha(static_cast<int>(i)) = sum_a;
  }

  // link vectors in world frame
  vector<Vector2d> r_link(n), r_com(n);
  for (size_t i = 0; i < n; ++i) {
    Matrix2d Ri = R(kc.a(static_cast<int>(i)));
    r_link[i] = Ri * Vector2d(mp.l(static_cast<int>(i)), 0.0);
    r_com[i]  = Ri * Vector2d(mp.lc(static_cast<int>(i)), 0.0);
  }

  // forward pass: linear accelerations of COM and joint origins
  vector<Vector2d> aP(n+1), aC(n);
  aP[0] = Vector2d(0.0, -mp.g); // base acceleration (gravity)
  for (size_t i = 0; i < n; ++i) {
    Vector2d term1 = R90(r_com[i]) * alpha(static_cast<int>(i));
    Vector2d term2 = R90(R90(r_com[i])) * (omega(static_cast<int>(i)) * omega(static_cast<int>(i)));
    aC[i] = aP[i] + term1 + term2;

    Vector2d t1 = R90(r_link[i]) * alpha(static_cast<int>(i));
    Vector2d t2 = R90(R90(r_link[i])) * (omega(static_cast<int>(i)) * omega(static_cast<int>(i)));
    aP[i+1] = aP[i] + t1 + t2;
  }

  // backward pass: accumulate forces/moments
  vector<Vector2d> F(n, Vector2d::Zero());
  vector<double> N(n, 0.0);

  for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
    Vector2d Fi = mp.m(i) * aC[i];
    double Ni = mp.I(i) * alpha(i) + (r_com[i].x() * Fi.y() - r_com[i].y() * Fi.x());
    if (i < static_cast<int>(n) - 1) {
      Vector2d Fchild = F[i+1];
      Fi += Fchild;
      Ni += N[i+1] + (r_link[i].x() * Fchild.y() - r_link[i].y() * Fchild.x());
    }
    F[i] = Fi;
    N[i] = Ni;
    tau(static_cast<int>(i)) = Ni;
  }

  // viscous damping
  tau.noalias() += mp.viscous * qdot;

  return tau;
}

/* ----------------- Mass matrix (CRBA-style via RNEA columns) ------------- */
static MatrixXd mass_matrix(const VectorXd &q, const ModelParams &mp) {
  const size_t n = NL;
  MatrixXd M = MatrixXd::Zero(n, n);
  VectorXd qdot = VectorXd::Zero(n);
  VectorXd qdd  = VectorXd::Zero(n);

  ModelParams mp_nod = mp;
  mp_nod.viscous = 0.0; // avoid viscous term when building inertia

  for (size_t k = 0; k < n; ++k) {
    qdd.setZero();
    qdd(static_cast<int>(k)) = 1.0;
    VectorXd col = RNEA(q, qdot, qdd, mp_nod); // returns M(:,k)
    M.col(static_cast<int>(k)) = col;
  }
  return M;
}

/* ----------------- Bias term h(q,qd) = C*qd + g(q) --------------------- */
static VectorXd bias_term(const VectorXd &q, const VectorXd &qd, const ModelParams &mp) {
  VectorXd qdd_zero = VectorXd::Zero(NL);
  return RNEA(q, qd, qdd_zero, mp);
}

/* ----------------- System post used by SCOTS ---------------------------
   Integrate one sampling interval tau using RK4 with the accurate RHS.
   RHS: xdot = [ qdot ; qdd ], where qdd = M(q)^{-1} ( u - h(q,qd) ).
------------------------------------------------------------------------- */
static auto system_post = [](state_type &x, input_type &u) {
  // local static model parameters
  static ModelParams mp;
  static bool init = false;
  if (!init) {
    // default params (tweak here if desired)
    mp.l  = VectorXd::Constant(NL, 1.0);
    mp.lc = 0.5 * mp.l;
    mp.m  = VectorXd::Constant(NL, 1.0);
    mp.I  = VectorXd::Constant(NL, 0.05);
    mp.g = 9.81;
    mp.viscous = 0.05;
    init = true;
  }

  // wrapper RHS to match runge_kutta_fixed4 signature
  auto rhs = [&](state_type &xdot, const state_type &xcur, input_type &uu) {
    // convert to Eigen
    VectorXd q(NL), qd(NL), tau_u(NL);
    for (size_t i = 0; i < NL; ++i) { q(static_cast<int>(i)) = xcur[i]; qd(static_cast<int>(i)) = xcur[NL + i]; }
    for (size_t i = 0; i < NL; ++i) { tau_u(static_cast<int>(i)) = uu[i]; }

    // compute M and h
    MatrixXd M = mass_matrix(q, mp);
    VectorXd h = bias_term(q, qd, mp);

    // solve for accelerations
    VectorXd qdd = M.ldlt().solve(tau_u - h);

    // fill xdot
    for (size_t i = 0; i < NL; ++i) {
      xdot[i] = qd(static_cast<int>(i));
      xdot[NL + i] = qdd(static_cast<int>(i));
    }
  };

  // number of RK substeps can be tuned; using 5 substeps for accuracy
  runge_kutta_fixed4(rhs, x, u, STATE_DIM, tau, 5);
};

/* ----------------- Radius propagation (conservative Jacobian-based) -----
   We compute a conservative linear bound on the local Jacobian A(x,u) of f(x,u)
   using numerical directional derivatives (cheap FD) at the center point.
   Then propagate radii via r_next = | r + tau * A * r | (entrywise abs).
   This is tighter than a fixed ad-hoc bound but still conservative.
------------------------------------------------------------------------- */
static auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
  // convert center state to Eigen
  VectorXd xc(STATE_DIM), uu(INPUT_DIM);
  for (size_t i = 0; i < STATE_DIM; ++i) xc(static_cast<int>(i)) = x[i];
  for (size_t j = 0; j < INPUT_DIM; ++j) uu(static_cast<int>(j)) = u[j];

  // small FD step
  const double eps = 1e-6;

  // Evaluate f(x,u) via same RHS as system_post's inner rhs
  auto eval_f = [&](const VectorXd &xx, const VectorXd &uu_in, VectorXd &fx_out) {
    // split
    VectorXd q(NL), qd(NL), tau_u(NL);
    for (size_t i = 0; i < NL; ++i) { q(static_cast<int>(i)) = xx(static_cast<int>(i)); qd(static_cast<int>(i)) = xx(static_cast<int>(NL + i)); }
    for (size_t i = 0; i < NL; ++i) { tau_u(static_cast<int>(i)) = uu_in(static_cast<int>(i)); }

    // local params (match those in system_post)
    ModelParams mp_local;
    mp_local.l  = VectorXd::Constant(NL, 1.0);
    mp_local.lc = 0.5 * mp_local.l;
    mp_local.m  = VectorXd::Constant(NL, 1.0);
    mp_local.I  = VectorXd::Constant(NL, 0.05);
    mp_local.g = 9.81; mp_local.viscous = 0.05;

    MatrixXd M = mass_matrix(q, mp_local);
    VectorXd h = bias_term(q, qd, mp_local);
    VectorXd qdd = M.ldlt().solve(tau_u - h);

    fx_out.resize(STATE_DIM);
    for (size_t i = 0; i < NL; ++i) {
      fx_out(static_cast<int>(i)) = qd(static_cast<int>(i));
      fx_out(static_cast<int>(NL + i)) = qdd(static_cast<int>(i));
    }
  };

  // base f0
  VectorXd f0(STATE_DIM);
  eval_f(xc, uu, f0);

  // Form A matrix approx by columnwise finite differences: A(:,j) ≈ (f(x+eps e_j) - f0)/eps
  MatrixXd A(STATE_DIM, STATE_DIM);
  A.setZero();
  VectorXd xx = xc;
  for (int j = 0; j < static_cast<int>(STATE_DIM); ++j) {
    VectorXd xp = xc;
    xp(j) += eps;
    VectorXd fp(STATE_DIM);
    eval_f(xp, uu, fp);
    A.col(j) = (fp - f0) / eps;
  }

  // Now compute r_next = | r + tau * A * r |
  VectorXd rvec(STATE_DIM);
  for (size_t i = 0; i < STATE_DIM; ++i) rvec(static_cast<int>(i)) = r[i];

  VectorXd res = rvec + tau * A.cwiseAbs() * rvec.cwiseAbs(); // use abs(A) * abs(r) to be conservative
  for (size_t i = 0; i < STATE_DIM; ++i) r[i] = abs(res(static_cast<int>(i)));
};

/* ----------------- Main: SCOTS grids, abstraction, synthesis, export ---- */
int main(int argc, char** argv) {
  TicToc tt;
  cout << "nlink_proper: N_LINKS=" << NL << " STATE_DIM=" << STATE_DIM << " INPUT_DIM=" << INPUT_DIM << endl;

  // State bounds & eta (per-dimension)
  state_type s_lb, s_ub, s_eta;
  for (size_t i = 0; i < NL; ++i) {
    s_lb[i] = -0.8; s_ub[i] = 0.8; s_eta[i] = 0.06;            // angles
  }
  for (size_t i = 0; i < NL; ++i) {
    s_lb[NL + i] = -1.0; s_ub[NL + i] = 1.0; s_eta[NL + i] = 0.06; // angular velocities
  }

  UniformGrid ss(STATE_DIM, s_lb, s_ub, s_eta);
  cout << "State grid:" << endl; ss.print_info();

  // Input grid
  input_type i_lb, i_ub, i_eta;
  for (size_t i = 0; i < NL; ++i) { i_lb[i] = -2.0; i_ub[i] = 2.0; i_eta[i] = 0.1; }
  UniformGrid is(INPUT_DIM, i_lb, i_ub, i_eta);
  cout << "Input grid:" << endl; is.print_info();

  cout << "Computing transition function (may take time depending on grid)...\n";
  TransitionFunction tf;
  Abstraction<state_type, input_type> abs(ss, is);
  tt.tic();
  abs.compute_gb(tf, system_post, radius_post);
  tt.toc();

  cout << "Transitions: " << tf.get_no_transitions() << endl;

  // Define a simple target set: first min(2,NL) joints in [0.5,0.7]
  auto target = [&ss, &s_eta](const abs_type &idx) {
    state_type x; ss.itox(idx, x);
    const double a_lo = 0.5, a_hi = 0.7;
    size_t dims = min<size_t>(2, NL);
    for (size_t d = 0; d < dims; ++d) {
      double v = x[d]; double half = s_eta[d] / 2.0;
      if (!(a_lo <= (v - half) && (v + half) <= a_hi)) return false;
    }
    return true;
  };
  write_to_file(ss, target, "target");

  // Synthesis
  cout << "Solving reachability (synthesis)...\n";
  tt.tic();
  WinningDomain win = solve_reachability_game(tf, target);
  tt.toc();
  cout << "Winning domain size: " << win.get_size() << endl;

  // Write controller
  cout << "Writing controller to nlink_proper.scs\n";
  StaticController controller(ss, is, move(win));
  if (write_to_file(controller, "nlink_proper")) cout << "Controller written.\n";

  // Export controller table CSV (state-input pairs)
  {
    ofstream csv("nlink_proper_table.csv");
    for (size_t i = 0; i < STATE_DIM; ++i) { csv << "x" << i << (i+1<STATE_DIM ? "," : ""); }
    for (size_t i = 0; i < INPUT_DIM; ++i) { csv << ",u" << i; }
    csv << "\n";

    state_type x;
    size_t rows = 0;
    for (abs_type si = 0; si < ss.size(); ++si) {
      ss.itox(si, x);
      vector<input_type> ctrls;
      try { ctrls = controller.get_control<state_type, input_type>(x); }
      catch (const runtime_error&) { continue; }
      if (!ctrls.empty()) {
        for (auto &uc : ctrls) {
          for (size_t d = 0; d < STATE_DIM; ++d) csv << x[d] << ",";
          for (size_t d = 0; d < INPUT_DIM; ++d) csv << uc[d] << (d+1==INPUT_DIM ? "\n" : ",");
          ++rows;
        }
      }
    }
    csv.close();
    cout << "Controller table exported (" << rows << " rows) -> nlink_proper_table.csv\n";
  }

  // Optional: closed-loop simulation to produce a trajectory CSV for plotting
  {
    const double dt_sim = 0.01;
    const int steps = 20000; // produces 200s of data at dt_sim=0.01 -> 20000 rows
    state_type xcur;
    // initial state: small negative angles, zero velocities (ensure it's inside grid)
    for (size_t i = 0; i < NL; ++i) xcur[i] = -0.5 + 0.2 * static_cast<double>(i);
    for (size_t i = 0; i < NL; ++i) xcur[NL + i] = 0.0;

    // reuse same ModelParams as in system_post
    ModelParams mpSim;
    mpSim.l  = VectorXd::Constant(NL, 1.0);
    mpSim.lc = 0.5 * mpSim.l;
    mpSim.m  = VectorXd::Constant(NL, 1.0);
    mpSim.I  = VectorXd::Constant(NL, 0.05);
    mpSim.g = 9.81; mpSim.viscous = 0.05;

    auto sim_rhs = [&](state_type &xdot, const state_type &xloc, input_type &uloc) {
      // identical to rhs in system_post but using mpSim
      VectorXd q(NL), qd(NL), tau_u(NL);
      for (size_t i = 0; i < NL; ++i) { q(static_cast<int>(i)) = xloc[i]; qd(static_cast<int>(i)) = xloc[NL + i]; }
      for (size_t i = 0; i < NL; ++i) { tau_u(static_cast<int>(i)) = uloc[i]; }

      MatrixXd M = mass_matrix(q, mpSim);
      VectorXd h = bias_term(q, qd, mpSim);
      VectorXd qdd = M.ldlt().solve(tau_u - h);

      for (size_t i = 0; i < NL; ++i) {
        xdot[i] = qd(static_cast<int>(i));
        xdot[NL + i] = qdd(static_cast<int>(i));
      }
    };

    ofstream traj("nlink_proper.csv");
    double t = 0.0;
    input_type ucur;
    for (int k = 0; k < steps; ++k) {
      // get controller input for current quantized state
      vector<input_type> ctrls;
      try { ctrls = controller.get_control<state_type, input_type>(xcur); }
      catch (const runtime_error&) { ctrls.clear(); }
      if (!ctrls.empty()) ucur = ctrls[0]; else { for (size_t i=0;i<INPUT_DIM;++i) ucur[i]=0.0; }

      // write t, state, input (no header)
      traj << t;
      for (size_t d = 0; d < STATE_DIM; ++d) traj << "," << xcur[d];
      for (size_t d = 0; d < INPUT_DIM; ++d) traj << "," << ucur[d];
      traj << "\n";

      // integrate one dt_sim (RK4 substeps)
      runge_kutta_fixed4(sim_rhs, xcur, ucur, STATE_DIM, dt_sim, 5);
      // wrap angles into [-pi,pi]
      for (size_t i = 0; i < NL; ++i) {
        double a = xcur[i];
        a = fmod(a + M_PI, 2.0*M_PI);
        if (a < 0) a += 2.0*M_PI;
        xcur[i] = a - M_PI;
      }
      t += dt_sim;
    }
    traj.close();
    cout << "Closed-loop trajectory written to nlink_proper.csv (" << steps << " rows)\n";
  }

  return 0;
}
