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
  runge_kutta_fixed4(rhs, x, u, STATE_DIM, tau, 100);
};

/* ----------------- Radius propagation for x' = u -----------------
   Linearization: f_x = 0  => r_next = |r| (no growth from dynamics)
------------------------------------------------------------------- */
auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
  state_type r_next;
  for (size_t i = 0; i < NL; ++i) r_next[i] = abs(r[i]);
  r = r_next;
};

/* ----------------- Main program ----------------- */
int main(int argc, char** argv) {
    TicToc tt;
    cout << "N_LINKS = " << NL << ", STATE_DIM = " << STATE_DIM << ", INPUT_DIM = " << INPUT_DIM << endl;

    // State grid: angles only
    state_type s_lb, s_ub, s_eta;
    for (size_t i = 0; i < NL; ++i) {
        s_lb[i]  =  0.0;   // angle lower bound (rad) //-0.8
        s_ub[i]  =  4.0;   // angle upper bound (rad) //0.8
        s_eta[i] =  0.8;  // cell size for angle i //0.06
    }
    UniformGrid ss(STATE_DIM, s_lb, s_ub, s_eta);
    cout << "Uniform grid details:" << endl;
    ss.print_info();

    // Input grid: u is angular velocity command
    input_type i_lb, i_ub, i_eta;
    for (size_t i = 0; i < NL; ++i) {
        i_lb[i]  = -1.6;   // min angular rate (rad/s) //-0.5
        i_ub[i]  =  1.6;   // max angular rate (rad/s) //0.5
        i_eta[i] =  0.9; // quantization for control //0.057
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
    auto target = [&ss, &s_eta](const abs_type& idx) {
        state_type x;
        ss.itox(idx, x);
        const double a_lo = 1; // 6:0-3, 7:1-3
        const double a_hi = 2.5; 
        bool cond = true;
        size_t dims_to_check = min<size_t>(2, NL);
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
    StaticController controller(ss, is, move(win));
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

    for (abs_type si = 0; si < ss.size(); ++si) {
        ss.itox(si, x);
        try {
            controls = controller.get_control<state_type, input_type>(x);
        } catch (const runtime_error &) {
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
