// // franka_nlink.cc
// #include <array>
// #include <vector>
// #include <cmath>
// #include <limits>
// #include <iostream>
// #include <fstream>
// #include <sys/time.h>
// #include <Eigen/Dense>
// #include <sys/resource.h>

// #include "scots.hh"
// #include "TicToc.hh"
// #include "RungeKutta4.hh"

// struct rusage usage;

// #ifndef N_LINKS
// #define N_LINKS 3
// #endif

// using namespace std;
// using namespace scots;

// constexpr size_t NL = (size_t)N_LINKS;
// constexpr size_t STATE_DIM = 2 * NL; // angles + angular velocities
// constexpr size_t INPUT_DIM = NL;     // one torque per joint
// constexpr double tau = 1.0;


// using state_type = array<double, STATE_DIM>;
// using input_type = array<double, INPUT_DIM>;
// using abs_type = abs_type;

// /* ----------------- Utilities / small helpers ----------------- */

// template <size_t N>
// void runge_kutta_wrap(state_type &x, const input_type &u,
//                       function<void(state_type&, const state_type&, input_type&)> rhs) {
//     // scots provided RK helper used previously: runge_kutta_fixed4
//     // But that helper expects rhs signature (state_type&, const state_type&, input_type&)
//     runge_kutta_fixed4(rhs, x, (input_type&)u, STATE_DIM, tau, 10);
// }

// /* ----------------- Dynamics for N_LINKS -----------------
//    We keep original 2-link exact dynamics when NL==2.
//    For NL>2 we use a simple chain of double integrators:
//      theta_i'     = omega_i
//      omega_i'     = u_i
//    This is a reasonable modular baseline that preserves functionality
//    and keeps the SCOTS interface unchanged.
// ----------------------------------------------------------- */

// auto system_post = [](state_type &x, input_type &u) -> void {
//   // Generic chain: for i in [0..NL-1]
//   // state ordering: [theta0, theta1, ..., theta{N-1}, omega0, omega1, ..., omega{N-1}]
//   auto rhs = [](state_type &xx, const state_type &x, input_type &u) -> void {
//     // first NL entries are theta derivatives = omega
//     for (size_t i = 0; i < NL; ++i) {
//       xx[i] = x[NL + i]; // theta_i' = omega_i
//     }
//     // angular accelerations = control torques (simple double integrator)
//     for (size_t i = 0; i < NL; ++i) {
//       xx[NL + i] = u[i];
//     }
//   };
//   runge_kutta_fixed4(rhs, x, u, STATE_DIM, tau, 10);
// };

// /* ----------------- Radius (linearization / bounding) -----------------
//    For NL==2 keep original radius_post. For general N use simple
//    Jacobian-based Euler step for chain model: r <- |r + tau * A * r|
// ------------------------------------------------------------------- */

// auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
//   // Generic linearization for chain double-integrator:
//   // state vector ordering as before [theta..., omega...]
//   // A = [ 0  I ; 0  0 ] so r_next = r + tau * A * r
//   // which reduces to:
//   // theta_rad_next_i = | r_theta_i + tau * r_omega_i |
//   // omega_rad_next_i = | r_omega_i |
//   state_type r_next;
//   for (size_t i = 0; i < NL; ++i) {
//     // theta part
//     double theta_r = r[i];
//     double omega_r = r[NL + i];
//     r_next[i] = abs( theta_r + tau * omega_r );
//     // omega part remains same magnitude (since acceleration linearization is identity to input)
//     r_next[NL + i] = abs( omega_r );
//   }
//   r = r_next;
// };

// /* ----------------- Grid / abstraction / synthesis helpers ----------------- */

// template <size_t N>
// UniformGrid make_state_grid(const state_type &lb, const state_type &ub, const state_type &eta) {
//   return UniformGrid(STATE_DIM, lb, ub, eta);
// }

// template <size_t M>
// UniformGrid make_input_grid(const input_type &lb, const input_type &ub, const input_type &eta) {
//   return UniformGrid(INPUT_DIM, lb, ub, eta);
// }

// /* ----------------- Main program ----------------- */

// int main(int argc, char** argv) {
//     TicToc tt;
//     cout << "N_LINKS = " << NL << ", STATE_DIM = " << STATE_DIM << ", INPUT_DIM = " << INPUT_DIM << endl;

//     // Example bounds/etas -- for general N we set same bounds per joint.
//     state_type s_lb;
//     state_type s_ub;
//     state_type s_eta;
//     for (size_t i = 0; i < NL; ++i) {
//         // angle bounds (radians)
//         s_lb[i] = -0.8;
//         s_ub[i] =  0.8;
//         s_eta[i] = 0.06;
//     }
//     for (size_t i = 0; i < NL; ++i) {
//         // angular velocity bounds
//         s_lb[NL + i] = -0.5;
//         s_ub[NL + i] =  0.5;
//         s_eta[NL + i] = 0.06;
//     }

//     UniformGrid ss(STATE_DIM, s_lb, s_ub, s_eta);
//     cout << "Uniform grid details:" << endl;
//     ss.print_info();

//     input_type i_lb, i_ub, i_eta;
//     for (size_t i = 0; i < NL; ++i) {
//         i_lb[i] = -1.0;
//         i_ub[i] =  1.0;
//         i_eta[i] = 0.057;
//     }

//     UniformGrid is(INPUT_DIM, i_lb, i_ub, i_eta);
//     cout << "Input grid details:" << endl;
//     is.print_info();

//     cout << "Computing the transition function: " << endl;
//     TransitionFunction tf;
//     Abstraction<state_type, input_type> abs(ss, is);

//     tt.tic();
//     abs.compute_gb(tf, system_post, radius_post);
//     tt.toc();

//     if(!getrusage(RUSAGE_SELF, &usage)) {
//         if (tf.get_no_transitions() > 0) {
//             cout << "Memory per transition: " << usage.ru_maxrss / (double)tf.get_no_transitions() << endl;
//         }
//     }
//     cout << "Number of transitions: " << tf.get_no_transitions() << endl;

//     // Example target set: we will set a small hyper-rectangle in the first two
//     // angle dimensions for demonstration. You may change to any goal spec.
//     auto target = [&ss, &s_eta](const abs_type& idx) {
//         state_type x;
//         ss.itox(idx, x);
//         // For demonstration, require first two angles to be in [0.5,0.7]
//         // if NL < 2 this condition will be adapted to just the first joint.
//         double a_lo = 0.5;
//         double a_hi = 0.7;
//         bool cond = true;
//         size_t dims_to_check = min<size_t>(2, NL);
//         for (size_t d = 0; d < dims_to_check; ++d) {
//             double v = x[d];
//             double half = s_eta[d] / 2.0;
//             if (!(a_lo <= (v - half) && (v + half) <= a_hi)) {
//                 cond = false;
//                 break;
//             }
//         }
//         return cond;
//     };

//     write_to_file(ss, target, "target");

//     cout << "\nSynthesis: " << endl;
//     tt.tic();
//     WinningDomain win = solve_reachability_game(tf, target);
//     tt.toc();
//     cout << "Winning domain size: " << win.get_size() << endl;

//     cout << "\nWrite controller to controller.scs \n";
//     StaticController controller(ss, is, move(win));
//     if (write_to_file(controller, "franka_nlink")) {
//         cout << "Done writing controller file. \n";
//     }

//     /* ================= CSV Export ================= */
//     ofstream csvfile("franka_nlink.csv");
//     csvfile << "x0";
//     for (size_t i = 1; i < STATE_DIM; ++i) csvfile << ",x" << i;
//     for (size_t i = 0; i < INPUT_DIM; ++i) csvfile << ",u" << i;
//     csvfile << "\n";

//     state_type x;
//     vector<input_type> controls;
//     size_t rows_written = 0;

//     for (abs_type si = 0; si < ss.size(); ++si) {
//         ss.itox(si, x);
//         try {
//             controls = controller.get_control<state_type, input_type>(x);
//         } catch (const runtime_error &e) {
//             continue;
//         }

//         if (!controls.empty()) {
//             for (size_t ui = 0; ui < controls.size(); ++ui) {
//                 input_type uc = controls[ui];
//                 for (size_t d = 0; d < STATE_DIM; ++d) {
//                     csvfile << x[d] << ",";
//                 }
//                 for (size_t d = 0; d < INPUT_DIM; ++d) {
//                     csvfile << uc[d] << (d + 1 == INPUT_DIM ? "\n" : ",");
//                 }
//                 ++rows_written;
//             }
//         }
//     }
//     csvfile.close();
//     cout << "State–input pairs written to franka_nlink.csv (rows: " << rows_written << ")\n";

//     return 0;
// }




// franka_nlink_single_integrator.cc
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <Eigen/Dense>
#include <sys/resource.h>

#include "scots.hh"
#include "TicToc.hh"
#include "RungeKutta4.hh"

struct rusage usage;

#ifndef N_LINKS
#define N_LINKS 7
#endif

using namespace std;
using namespace scots;

constexpr size_t NL = (size_t)N_LINKS;
constexpr size_t STATE_DIM = NL;     // ONLY joint angles
constexpr size_t INPUT_DIM = NL;     // one velocity command per joint
constexpr double tau = 1.0;

using state_type = array<double, STATE_DIM>;
using input_type = array<double, INPUT_DIM>;

/* ----------------- Dynamics: single integrator x' = u ----------------- */
auto system_post = [](state_type &x, input_type &u) -> void {
  // rhs: xdot = u  (independent of x)
  auto rhs = [](state_type &xdot, const state_type &x, input_type &uu) -> void {
    for (size_t i = 0; i < NL; ++i) xdot[i] = uu[i];
  };
  runge_kutta_fixed4(rhs, x, u, STATE_DIM, tau, 10);
};

/* ----------------- Radius propagation for x' = u -----------------
   Linearization: f_x = 0  => r_next = |r| (no growth from dynamics)
------------------------------------------------------------------- */
auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
  state_type r_next;
  for (size_t i = 0; i < NL; ++i) r_next[i] = std::abs(r[i]);
  r = r_next;
};

/* ----------------- Main program ----------------- */
int main(int argc, char** argv) {
    TicToc tt;
    cout << "N_LINKS = " << NL << ", STATE_DIM = " << STATE_DIM << ", INPUT_DIM = " << INPUT_DIM << endl;

    // State grid: angles only
    state_type s_lb, s_ub, s_eta;
    for (size_t i = 0; i < NL; ++i) {
        s_lb[i]  = -0.8;   // angle lower bound (rad)
        s_ub[i]  =  0.8;   // angle upper bound (rad)
        s_eta[i] =  0.27;  // cell size for angle i //0.06
    }
    UniformGrid ss(STATE_DIM, s_lb, s_ub, s_eta);
    cout << "Uniform grid details:" << endl;
    ss.print_info();

    // Input grid: u is angular velocity command
    input_type i_lb, i_ub, i_eta;
    for (size_t i = 0; i < NL; ++i) {
        i_lb[i]  = -0.5;   // min angular rate (rad/s)
        i_ub[i]  =  0.5;   // max angular rate (rad/s)
        i_eta[i] =  0.27; // quantization for control //0.057
    }
    UniformGrid is(INPUT_DIM, i_lb, i_ub, i_eta);
    cout << "Input grid details:" << endl;
    is.print_info();

    cout << "Computing the transition function: " << endl;
    TransitionFunction tf;
    Abstraction<state_type, input_type> abs(ss, is);

    tt.tic();
    abs.compute_gb(tf, system_post, radius_post);
    tt.toc();

    if(!getrusage(RUSAGE_SELF, &usage)) {
        if (tf.get_no_transitions() > 0) {
            cout << "Memory per transition: "
                 << usage.ru_maxrss / (double)tf.get_no_transitions() << endl;
        }
    }
    cout << "Number of transitions: " << tf.get_no_transitions() << endl;

    // Example target set: first up to 2 angles in [0.5, 0.7]
    auto target = [&ss, &s_eta](const scots::abs_type& idx) {
        state_type x;
        ss.itox(idx, x);
        const double a_lo = 0.5;
        const double a_hi = 0.7;
        bool cond = true;
        size_t dims_to_check = std::min<size_t>(2, NL);
        for (size_t d = 0; d < dims_to_check; ++d) {
            double v = x[d];
            double half = s_eta[d] / 2.0;
            if (!(a_lo <= (v - half) && (v + half) <= a_hi)) {
                cond = false;
                break;
            }
        }
        return cond;
    };

    write_to_file(ss, target, "target");

    cout << "\nSynthesis: " << endl;
    tt.tic();
    WinningDomain win = solve_reachability_game(tf, target);
    tt.toc();
    cout << "Winning domain size: " << win.get_size() << endl;

    cout << "\nWrite controller to controller.scs \n";
    StaticController controller(ss, is, std::move(win));
    if (write_to_file(controller, "franka_nlink_single_integrator")) {
        cout << "Done writing controller file. \n";
    }

    /* ================= CSV Export ================= */
    ofstream csvfile("franka_nlink_single_integrator.csv");
    // header
    csvfile << "x0";
    for (size_t i = 1; i < STATE_DIM; ++i) csvfile << ",x" << i;
    for (size_t i = 0; i < INPUT_DIM; ++i) csvfile << ",u" << i;
    csvfile << "\n";

    state_type x;
    vector<input_type> controls;
    size_t rows_written = 0;

    for (scots::abs_type si = 0; si < ss.size(); ++si) {
        ss.itox(si, x);
        try {
            controls = controller.get_control<state_type, input_type>(x);
        } catch (const std::runtime_error &) {
            continue;
        }
        for (const auto& uc : controls) {
            for (size_t d = 0; d < STATE_DIM; ++d) csvfile << x[d] << ",";
            for (size_t d = 0; d < INPUT_DIM; ++d)
                csvfile << uc[d] << (d + 1 == INPUT_DIM ? "\n" : ",");
            ++rows_written;
        }
    }
    csvfile.close();
    cout << "State–input pairs written to franka_nlink_single_integrator.csv (rows: "
         << rows_written << ")\n";

    return 0;
}
