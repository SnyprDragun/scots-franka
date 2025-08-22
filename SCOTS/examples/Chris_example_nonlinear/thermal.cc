/*
 * dcdc.cc
 *
 *  created: Oct 2015
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in the readme file
 */

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
const int state_dim=1;
/* input space dim */
const int input_dim=1;
/* sampling time */
const double tau = 0.01;

/*
 * data types for the elements of the state space 
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* parameters for system dynamics */
const double b=0.1;

/* parameters for radius calculation */
const double mu=std::sqrt(2);
/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
   /* the ode describing the dcdc converter */
   auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) noexcept {
       xx[0]=-b*sin(2*x[0])-sin(x[0])*sin(x[0])+cos(x[0])*cos(x[0])+u[0]*cos(x[0]);
     };
   scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);
};
/* we integrate the growth bound by 0.5 sec (the result is stored in r)  */

auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
	/* the ode for the growth bound */
    auto rhs =[](state_type& rr,  const state_type &r, const input_type &u) noexcept {
      rr[0] = std::abs(-2*b-4-u[0])*r[0];
      };
    scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{0.01}};
 /* lower bounds of the hyper-rectangle */
  state_type lb={{-1.1071}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{1.1071}};
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-1}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{1}};
  /* grid node distance diameter */
  input_type i_eta={{0.05}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* construct grid for the input alphabet */
  /* hyper-rectangle [1,2] with grid node distance 1 */
  /*scots::UniformGrid is(input_dim,input_type{{1}},input_type{{2}},input_type{{1}});
  is.print_info();*/

  /* compute transition function of symbolic model */
  std::cout << "Computing the transition function:\n";
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  abs.verbose_off();

  tt.tic();
  abs.compute_gb(tf,system_post, radius_post);
  tt.toc();
  std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";

  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n";

  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0])
      return true;
    return false;
  };
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf,safeset);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << "\n";

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";

  return 1;
}

