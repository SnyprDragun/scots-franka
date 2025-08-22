// franka_nlink_slow.cc
// Generic n-link planar rigid-body dynamics + SCOTS abstraction/controller export.

#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>

#include <Eigen/Dense>

#include "scots.hh"
#include "TicToc.hh"
#include "RungeKutta4.hh"

#ifndef N_LINKS
#define N_LINKS 3
#endif

using namespace std;
using namespace scots;

// ---------- Typedefs / compile-time sizes ----------
constexpr size_t NL         = static_cast<size_t>(N_LINKS);
constexpr size_t STATE_DIM  = 2 * NL;     // [q(0..n-1), qdot(0..n-1)]
constexpr size_t INPUT_DIM  = NL;         // one torque per joint
constexpr double tau        = 1.0;        // sampling time used by SCOTS integrator (sec)

using state_type = array<double, STATE_DIM>;
using input_type = array<double, INPUT_DIM>;

// ---------- Planar model parameters ----------
struct ModelParams {
  // link lengths (joint i-1 to joint i)
  Eigen::VectorXd l;      // size n
  // COM distance from joint i-1 along link i
  Eigen::VectorXd lc;     // size n
  // masses and planar (about z through COM) inertias
  Eigen::VectorXd m;      // size n
  Eigen::VectorXd I;      // size n
  double g = 9.81;        // gravity
  double viscous = 0.05;  // viscous damping coefficient (Nm·s/rad), uniform (can be vectorized)
};

// Rotation by angle a (2D)
static inline Eigen::Matrix2d R(double a) {
  double c = std::cos(a), s = std::sin(a);
  Eigen::Matrix2d M; M << c, -s, s, c;
  return M;
}

// R90 * v (perp operator)
static inline Eigen::Vector2d R90(const Eigen::Vector2d &v) {
  // Rotates v by +90 degrees: [0 -1; 1 0] * v
  return Eigen::Vector2d(-v.y(), v.x());
}

// ---------- Kinematics helpers (world frame) ----------
struct KinCache {
  // Joint positions P_i in world (P_0 = base at origin)
  vector<Eigen::Vector2d> P;  // size n+1 (0..n); P_i = joint i position, P_0 = base, P_i is end of link i
  // Link COM positions C_i in world (1..n)
  vector<Eigen::Vector2d> C;  // size n (1..n mapped to [0..n-1])
  // Cumulative angles a_i = sum_{k=1..i} q_k
  Eigen::VectorXd a;          // size n
};

// Build joint and COM positions for given q
static KinCache forward_positions(const Eigen::VectorXd &q, const ModelParams &mp) {
  KinCache kc;
  kc.P.resize(NL + 1);
  kc.C.resize(NL);
  kc.a.setZero(NL);

  Eigen::Vector2d base = Eigen::Vector2d::Zero();
  kc.P[0] = base;

  double a_sum = 0.0;
  for (size_t i = 0; i < NL; ++i) {
    a_sum += q(static_cast<int>(i));
    kc.a(static_cast<int>(i)) = a_sum;

    Eigen::Matrix2d Ri = R(a_sum);
    // joint i position to next joint i+1
    Eigen::Vector2d r_link = Ri * Eigen::Vector2d(mp.l(static_cast<int>(i)), 0.0);
    // COM offset from joint i position
    Eigen::Vector2d r_com  = Ri * Eigen::Vector2d(mp.lc(static_cast<int>(i)), 0.0);

    kc.C[i]   = kc.P[i] + r_com;
    kc.P[i+1] = kc.P[i] + r_link;
  }
  return kc;
}

// ---------- Recursive Newton–Euler (planar) ----------
// Computes joint torques tau for given (q, qdot, qddot).
// World-frame formulation using positions/angles above.
// Gravity acts at the base as a_P0 = [0, -g].
static Eigen::VectorXd RNEA(const Eigen::VectorXd &q,
                            const Eigen::VectorXd &qdot,
                            const Eigen::VectorXd &qddot,
                            const ModelParams &mp) {
  const size_t n = NL;
  Eigen::VectorXd tau(n); tau.setZero();

  // Kinematic pass to get positions (also reuse for shifting wrenches)
  KinCache kc = forward_positions(q, mp);

  // For each link i, compute:
  // angular vel/acc: omega_i = sum_{k<=i} qdot_k ; alpha_i = sum_{k<=i} qddot_k
  Eigen::VectorXd omega(n), alpha(n);
  omega.setZero(); alpha.setZero();
  double sum_w = 0.0, sum_a = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum_w += qdot(static_cast<int>(i));
    sum_a += qddot(static_cast<int>(i));
    omega(static_cast<int>(i)) = sum_w;
    alpha(static_cast<int>(i)) = sum_a;
  }

  // Base linear acceleration (gravity)
  Eigen::Vector2d aP_prev(0.0, -mp.g); // acceleration at joint 0 (world origin)
  // We also need the *joint* positions P_i to propagate to next joint
  // Already in kc.P. For convenience, we will also need link vectors (joint->next joint)
  vector<Eigen::Vector2d> r_link(n), r_com(n);
  for (size_t i = 0; i < n; ++i) {
    Eigen::Matrix2d Ri = R(kc.a(static_cast<int>(i)));
    r_link[i] = Ri * Eigen::Vector2d(mp.l(static_cast<int>(i)), 0.0);
    r_com[i]  = Ri * Eigen::Vector2d(mp.lc(static_cast<int>(i)), 0.0);
  }

  // Forward pass: compute linear accelerations of each COM and next joint origins
  vector<Eigen::Vector2d> aP(n+1); // joint accelerations
  vector<Eigen::Vector2d> aC(n);   // COM accelerations
  aP[0] = aP_prev;

  for (size_t i = 0; i < n; ++i) {
    // Acceleration from joint i to COM of link i:
    // a_Ci = a_Pi + alpha_i x r_com + omega_i x (omega_i x r_com)
    Eigen::Vector2d term1 = R90(r_com[i]) * alpha(static_cast<int>(i));
    Eigen::Vector2d term2 = R90(R90(r_com[i])) * (omega(static_cast<int>(i)) * omega(static_cast<int>(i)));
    aC[i] = aP[i] + term1 + term2;

    // Acceleration to next joint origin:
    Eigen::Vector2d t1 = R90(r_link[i]) * alpha(static_cast<int>(i));
    Eigen::Vector2d t2 = R90(R90(r_link[i])) * (omega(static_cast<int>(i)) * omega(static_cast<int>(i)));
    aP[i+1] = aP[i] + t1 + t2;
  }

  // Backward pass: accumulate forces and moments to get joint torques
  vector<Eigen::Vector2d> F(n, Eigen::Vector2d::Zero()); // net force at joint i from link i..n
  vector<double>          N(n, 0.0);                     // net moment about joint i (z-scalar)

  for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
    // Link i's own inertial force and moment about joint i:
    Eigen::Vector2d Fi = mp.m(i) * aC[i];                      // linear
    double Ni = mp.I(i) * alpha(i) + (r_com[i].x()*Fi.y() - r_com[i].y()*Fi.x()); // torque about joint i

    // Add child (i+1) contributions shifted to joint i
    if (i < static_cast<int>(n) - 1) {
      // Shift child net force F[i+1] from joint (i+1) to joint i
      Eigen::Vector2d Fchild = F[i+1];
      Fi += Fchild;
      // Add moment contribution: Nchild + (r_link[i] x Fchild)
      Ni += N[i+1] + (r_link[i].x()*Fchild.y() - r_link[i].y()*Fchild.x());
    }

    // Store totals at joint i
    F[i] = Fi;
    N[i] = Ni;

    // Joint torque is z-component of moment at joint i
    tau(i) = Ni;
  }

  // Add simple viscous damping at joints (Nm = d * qdot)
  tau.noalias() += mp.viscous * qdot;

  return tau;
}

// Build mass matrix M(q) using RNEA columns: tau = M(q) * qddot + h(q,qdot)
// With qdot = 0 and qddot = e_k, tau = M.col(k)
static Eigen::MatrixXd mass_matrix(const Eigen::VectorXd &q, const ModelParams &mp) {
  const size_t n = NL;
  Eigen::MatrixXd M(n, n);
  M.setZero();

  Eigen::VectorXd qdot = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd qdd  = Eigen::VectorXd::Zero(n);

  for (size_t k = 0; k < n; ++k) {
    qdd.setZero();
    qdd(static_cast<int>(k)) = 1.0;
    // No damping when building M (set viscous=0 temporarily)
    ModelParams mp_nod = mp;
    mp_nod.viscous = 0.0;
    Eigen::VectorXd tau_col = RNEA(q, qdot, qdd, mp_nod); // = M(:,k)
    M.col(static_cast<int>(k)) = tau_col;
  }
  return M;
}

// Bias term h(q,qdot) = RNEA(q, qdot, 0)  (includes gravity + Coriolis/centrifugal + viscous)
static Eigen::VectorXd bias_term(const Eigen::VectorXd &q,
                                 const Eigen::VectorXd &qdot,
                                 const ModelParams &mp) {
  Eigen::VectorXd qdd = Eigen::VectorXd::Zero(NL);
  return RNEA(q, qdot, qdd, mp);
}

// Convert between SCOTS arrays and Eigen vectors
static inline void to_eigen(const state_type &x, Eigen::VectorXd &q, Eigen::VectorXd &qd) {
  q.resize(NL); qd.resize(NL);
  for (size_t i = 0; i < NL; ++i) {
    q(static_cast<int>(i))  = x[i];
    qd(static_cast<int>(i)) = x[NL + i];
  }
}
static inline void from_eigen(state_type &x, const Eigen::VectorXd &q, const Eigen::VectorXd &qd) {
  for (size_t i = 0; i < NL; ++i) {
    x[i]       = q(static_cast<int>(i));
    x[NL + i]  = qd(static_cast<int>(i));
  }
}
static inline void u_to_eigen(const input_type &u, Eigen::VectorXd &tau_u) {
  tau_u.resize(NL);
  for (size_t i = 0; i < NL; ++i) tau_u(static_cast<int>(i)) = u[i];
}

// ---------- Continuous-time RHS: xdot = f(x,u) ----------
static void rhs(state_type &xdot, const state_type &x, input_type &u, const ModelParams &mp) {
  Eigen::VectorXd q, qd;
  to_eigen(x, q, qd);
  Eigen::VectorXd tau_u; u_to_eigen(u, tau_u);

  // M(q) qdd + h(q,qd) = tau_u
  Eigen::MatrixXd M = mass_matrix(q, mp);
  Eigen::VectorXd h = bias_term(q, qd, mp);

  Eigen::VectorXd qdd = M.ldlt().solve(tau_u - h);

  // Assemble xdot
  for (size_t i = 0; i < NL; ++i) {
    xdot[i]       = qd(static_cast<int>(i));      // qdot
    xdot[NL + i]  = qdd(static_cast<int>(i));     // qddot
  }
}

// ---------- System post (integrate one tau with RK4) ----------
static auto system_post = [](state_type &x, input_type &u) {
  // Wrap RHS for runge_kutta_fixed4
  auto wrapped = [&](state_type &xdot, const state_type &xcur, input_type &uu) {
    // Model parameters (could be made global, but kept here for clarity)
    static ModelParams mp;
    static bool init = false;
    if (!init) {
      mp.l     = Eigen::VectorXd::Constant(NL, 1.0);     // link lengths
      mp.lc    = 0.5 * mp.l;                              // COM at mid-link
      mp.m     = Eigen::VectorXd::Constant(NL, 1.0);     // 1 kg per link
      mp.I     = Eigen::VectorXd::Constant(NL, 0.05);    // modest planar inertia
      mp.g     = 9.81;
      mp.viscous = 0.05;
      init = true;
    }
    rhs(xdot, xcur, uu, mp);
  };
  runge_kutta_fixed4(wrapped, x, u, STATE_DIM, tau, 5); // 5 substeps per sample (tune)
};

// ---------- Radius propagation (conservative generic bound) ----------
// For a nonlinear manipulator, a tight bound requires local Lipschitz constants.
// As a simple conservative bound, we use:
//   r_q_next     = | r_q + tau * r_qdot |
//   r_qdot_next  = | r_qdot | + tau * ( a0 + a1 * |r_q| + a2 * |r_qdot| )
// with small (a0,a1,a2) to avoid zero growth; tune as needed for soundness.
static auto radius_post = [](state_type &r, const state_type &, const input_type &) {
  const double a0 = 1e-4, a1 = 0.5, a2 = 0.1; // conservative tunables
  state_type rnext;
  for (size_t i = 0; i < NL; ++i) {
    double rq   = std::abs(r[i]);
    double rqd  = std::abs(r[NL + i]);
    rnext[i]       = std::abs(r[i] + tau * r[NL + i]);
    rnext[NL + i]  = rqd + tau * (a0 + a1 * rq + a2 * rqd);
  }
  r = rnext;
};

// ---------- Main (SCOTS grids, abstraction, synthesis, export) ----------
int main(int, char**) {
  struct rusage usage;
  TicToc tt;

  cout << "N_LINKS = " << NL
       << "  STATE_DIM = " << STATE_DIM
       << "  INPUT_DIM = " << INPUT_DIM << "\n";

  // State grid bounds/etas (uniform per joint; tune per task)
  state_type s_lb, s_ub, s_eta;
  for (size_t i = 0; i < NL; ++i) {
    s_lb[i]      = -0.8;   s_ub[i]      = 0.8;    s_eta[i]      = 0.06; // angles
    s_lb[NL + i] = -0.8;   s_ub[NL + i] = 0.8;    s_eta[NL + i] = 0.06; // velocities
  }
  UniformGrid ss(STATE_DIM, s_lb, s_ub, s_eta);
  ss.print_info();

  // Input grid (torques)
  input_type i_lb, i_ub, i_eta;
  for (size_t i = 0; i < NL; ++i) {
    i_lb[i] = -2.0; i_ub[i] = 2.0; i_eta[i] = 0.1;
  }
  UniformGrid is(INPUT_DIM, i_lb, i_ub, i_eta);
  is.print_info();

  cout << "Computing transition function...\n";
  TransitionFunction tf;
  Abstraction<state_type, input_type> abs(ss, is);

  tt.tic();
  abs.compute_gb(tf, system_post, radius_post);
  tt.toc();

  if (!getrusage(RUSAGE_SELF, &usage) && tf.get_no_transitions() > 0) {
    cout << "Memory per transition (approx): "
         << usage.ru_maxrss / static_cast<double>(tf.get_no_transitions()) << " [kB/trans]\n";
  }
  cout << "No. of transitions: " << tf.get_no_transitions() << "\n";

  // Example target: require first two angles (or fewer if n==1) in [0.5, 0.7]
  auto target = [&ss, &s_eta](const abs_type& idx) {
    state_type x; ss.itox(idx, x);
    const double a_lo = 0.5, a_hi = 0.7;
    size_t dims_to_check = std::min<size_t>(2, NL);
    for (size_t d = 0; d < dims_to_check; ++d) {
      double v = x[d], half = s_eta[d] / 2.0;
      if (!(a_lo <= (v - half) && (v + half) <= a_hi))
        return false;
    }
    return true;
  };
  write_to_file(ss, target, "target");

  cout << "Synthesis (reachability)...\n";
  tt.tic();
  WinningDomain win = solve_reachability_game(tf, target);
  tt.toc();
  cout << "Winning domain size: " << win.get_size() << "\n";

  cout << "Write controller to franka_nlink_slow.scs ...\n";
  StaticController controller(ss, is, std::move(win));
  if (write_to_file(controller, "franka_nlink_slow")) {
    cout << "Done writing controller.\n";
  }

  // -------- CSV export (state-input map) --------
  ofstream csv("franka_nlink_slow.csv");
  for (size_t i = 0; i < STATE_DIM; ++i) csv << "x" << i << (i+1<STATE_DIM ? "," : "");
  for (size_t i = 0; i < INPUT_DIM; ++i) csv << ",u" << i;
  csv << "\n";

  state_type x;
  size_t rows_written = 0;
  for (abs_type si = 0; si < ss.size(); ++si) {
    ss.itox(si, x);
    vector<input_type> ctrls;
    try {
      ctrls = controller.get_control<state_type, input_type>(x);
    } catch (const runtime_error&) {
      continue;
    }
    if (!ctrls.empty()) {
      for (auto &uc : ctrls) {
        for (size_t d = 0; d < STATE_DIM; ++d) csv << x[d] << ",";
        for (size_t d = 0; d < INPUT_DIM; ++d)
          csv << uc[d] << (d+1==INPUT_DIM ? "\n" : ",");
        ++rows_written;
      }
    }
  }
  csv.close();
  cout << "Exported " << rows_written << " rows to franka_nlink_slow.csv\n";

  return 0;
}
