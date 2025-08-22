#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <sys/time.h>
#include <Eigen/Dense>
#include <sys/resource.h>

#include "scots.hh"
#include "TicToc.hh"
#include "RungeKutta4.hh"

struct rusage usage;
const int state_dim=4; // state space dim
const int input_dim=2; // input space dim
const double tau=1;

using namespace std;
using namespace scots;
using state_type = array<double,state_dim>;
using input_type = array<double,input_dim>;
using abs_type = abs_type;

auto system_post = [](state_type &x, input_type &u) -> void {
  auto rhs =[](state_type& xx, const state_type &x, input_type &u) -> void {
    double x1,x2,x3,x4,u1,u2;
    x1=x[0];
    x2=x[1];
    x3=x[2];
    x4=x[3];
    u1=u[0];
    u2=u[1];
    xx[0] = x[2];
    xx[1] = x[3];
    xx[2]=-(3*(4*u1 - 4*u2 + 2*pow(x3,2)*sin(x2) + 2*pow(x4,2)*sin(x2) - 6*u2*cos(x2) + 4*x3*x4*sin(x2) + 3*pow(x3,2)*cos(x2)*sin(x2)))/(9*pow(cos(x2),2) - 16);
    xx[3]=(6*u1*(3*cos(x2) + 2))/(9*pow(cos(x2),2) - 16) - (12*u2*(3*cos(x2) + 5))/(9*pow(cos(x2),2) - 16) + (6*pow(x3,2)*sin(x2)*(3*cos(x2) + 5))/(9*pow(cos(x2),2) - 16) + (3*x4*sin(x2)*(3*cos(x2) + 2)*(2*x3 + x4))/(9*pow(cos(x2),2) - 16);
  };
  runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};

auto radius_post = [](state_type &r, const state_type &x, const input_type &u) {
	double x1=x[0];
	double x2=x[1];
	double x3=x[2];
	double x4=x[3];
	double u1=u[0];
	double u2=u[1];
  double A_3_3= -(6*sin(x2)*(2*x3 + 2*x4 + 3*x3*cos(x2)))/(9*pow(cos(x2),2) - 16);
  double A_3_4=(12*abs(sin(x2)*(x3 + x4)))/abs(9*pow(sin(x2),2) + 7);
  double A_3_2=abs(108*u2*sin(2*x2) - 12*pow(x4,2)*cos(x2) - 108*u1*sin(2*x2) - 12*pow(x3,2)*cos(x2) - 162*u2*pow(sin(x2),3) + 207*pow(x3,2)*pow(cos(x2),2) + 54*pow(x3,2)*pow(cos(x2),3) + 54*pow(x4,2)*pow(cos(x2),3) + 450*u2*sin(x2) - 144*pow(x3,2) - 24*x3*x4*cos(x2) + 108*x3*x4*pow(cos(x2),3))/pow(abs(9*pow(cos(x2),2) - 16),2);
  double A_4_4=(6*sin(x2)*(x3 + x4)*(3*cos(x2) + 2))/(9*pow(cos(x2),2) - 16);
  double A_4_3=(6*abs(sin(x2)*(10*x3 + 2*x4 + 6*x3*cos(x2) + 3*x4*cos(x2))))/abs(9*pow(cos(x2),2) - 16);
  double A_4_2=abs(1140*pow(x3,2)*cos(x2) + 228*pow(x4,2)*cos(x2) - 864*u1*sin(2*x2) - 324*u1*sin(3*x2) + 4320*u2*sin(2*x2) + 648*u2*sin(3*x2) - 648*x3*x4 + 1656*pow(x3,2)*cos(2*x2) + 540*pow(x3,2)*cos(3*x2) + 828*pow(x4,2)*cos(2*x2) + 108*pow(x4,2)*cos(3*x2) - 2628*u1*sin(x2) + 5256*u2*sin(x2) - 648*pow(x3,2)- 324*pow(x4,2) + 456*x3*x4*cos(x2) + 1656*x3*x4*cos(2*x2) + 216*x3*x4*cos(3*x2))/(2*pow(abs(9*cos(2*x2) - 23),2));

  Eigen::Matrix4d A;
  A << 0, 0, 1, 0,
      0, 0, 0, 1,
      0, A_3_2, A_3_3, A_3_4,
      0, A_4_2,A_4_3, A_4_4;

  Eigen::Vector4d r_eigen;
  for (int i = 0; i < 4; i++) {
      r_eigen(i) = r[i];
  }
  Eigen::Vector4d result = r_eigen + tau*A * r_eigen;
  for (int i = 0; i < 4; i++) {
      r[i] = abs(result(i));
  }
};

int main() {
  TicToc tt;
  state_type s_lb={{-0.8,-0.8,-0.5,-0.5}};
  state_type s_ub={{0.8,0.8,0.5,0.5}};
  state_type s_eta={{0.06,0.06,0.06,0.06}};

  UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  cout << "Uniform grid details:" << endl;
  ss.print_info();

  input_type i_lb={{-1,-1}};
  input_type i_ub={{ 1, 1}};
  input_type i_eta={{0.057,0.057}};

  UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  cout << "Computing the transition function: " << endl;
  TransitionFunction tf;
  Abstraction<state_type,input_type> abs(ss,is);

  tt.tic();
  abs.compute_gb(tf,system_post, radius_post);
  tt.toc();

  if(!getrusage(RUSAGE_SELF, &usage)) {
    cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << endl;
  }
  cout << "Number of transitions: " << tf.get_no_transitions() << endl;

  auto target = [&ss, &s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx, x);
    if((0.5 <= (x[0] - s_eta[0] / 2.0) && (x[0] + s_eta[0] / 2.0) <= 0.7) && (0.5 <= (x[1] - s_eta[1] / 2.0) && (x[1] + s_eta[1] / 2.0) <=0.7)) {
      return true;
    }
    return false;
  };

  write_to_file(ss,target,"target");
  cout << "\nSynthesis: " << endl;
  tt.tic();
  WinningDomain win=solve_reachability_game(tf,target);
  tt.toc();
  cout << "Winning domain size: " << win.get_size() << endl;

  cout << "\nWrite controller to controller.scs \n";
  StaticController controller(ss, is, move(win));
  if (write_to_file(controller, "franka")) {
    cout << "Done. \n";
  }

  /* ================= CSV Export ================= */
  ofstream csvfile("franka.csv");
  csvfile << "x0,x1,x2,x3,u0,u1\n";

  state_type x;
  vector<input_type> controls;
  size_t rows_written = 0;

  for (abs_type si = 0; si < ss.size(); ++si) {
    ss.itox(si, x);
    try {
      controls = controller.get_control<state_type, input_type>(x);
    } 
    catch (const runtime_error &e) {
      continue;
    }

    if (!controls.empty()) {
      for (size_t ui = 0; ui < controls.size(); ++ui) {
        input_type uc = controls[ui];
        csvfile << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << ","
                << uc[0] << "," << uc[1] << "\n";
        ++rows_written;
      }
    }
  }
  csvfile.close();
  cout << "State–input pairs written to state_input_pairs.csv (rows: "
       << rows_written << ")\n";

  return 0;
}
