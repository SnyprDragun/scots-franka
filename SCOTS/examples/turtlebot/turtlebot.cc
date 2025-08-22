/*
 * dcdc.cc
 *
 *  created: Oct 2015
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in the readme file
 */
#include <limits>
#include <iostream>
#include <array>
#include <cmath>

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


/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;

const double tau=0.5;
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

using abs_type = scots::abs_type;

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto  system_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };
  // ode_solver(rhs,x,u);
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};



auto radius_post = [](state_type &r, const state_type &, const input_type &u) {
				const state_type w = {{0.01, 0.01}};
				r[0] = r[0] + r[2] * std::abs(u[0]) * tau + w[0];
				r[1] = r[1] + r[2] * std::abs(u[0]) * tau + w[1];
			};


int main() {
  TicToc tt;
  state_type s_lb={{0,0,-3.2}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{3,3,3.2}};
  /* grid node distance diameter */
  state_type s_eta={{.05,.05,.2}};

  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();

  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-0.22,-2}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{ 0.22, 2}};
  /* grid node distance diameter */
  input_type i_eta={{0.02,0.2}};

  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* set up constraint functions with obtacles */
  double H[1][4] = {

    { 1.3, 1.6  ,1.3  ,   1.6 },

  };

  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&H,ss,s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    double c1= s_eta[0]/2.0+1e-10;
    double c2= s_eta[1]/2.0+1e-10;
    for(size_t i=0; i<3; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0] <= (H[i][1]+c1) &&
          (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };
  /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");

  std::cout << "Computing the transition function: " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);

  tt.tic();
  abs.compute_gb(tf,system_post, radius_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (2 <= (x[0]-s_eta[0]/2.0) && (x[0]+s_eta[0]/2.0) <= 2.6 &&
        2 <= (x[1]-s_eta[1]/2.0) && (x[1]+s_eta[1]/2.0) <= 2.6)
      return true;
    return false;
  };
   /* write target to file */
  write_to_file(ss,target,"target");


  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::solve_reachability_game(tf,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"scots"))
    std::cout << "Done. \n";

}
