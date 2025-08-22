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
const int state_dim=10;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 1;

/*
 * data types for the elements of the state space 
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* parameters for system dynamics */
const double a=0.05;
const double ae2=0.005;
const double ae5=0.005;
const double ae=0.0033;
const double ah=0.0036;
const double te=12;
const double th=100;

/* parameters for radius calculation */
const double mu=std::sqrt(2);
/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
  /* the ode describing the dcdc converter */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) noexcept {
    xx[0]=(-a-ae)*x[0]+a*x[1]+ae*te;
    xx[1]=(-4*a-ae2-ah*u[0])*x[1]+a*x[0]+a*x[6]+a*x[8]+a*x[2]+ae2*te+ah*th*u[0];
    xx[2]=(-2*a-ae)*x[2]+a*x[1]+a*x[3]+ae*te;
    xx[3]=(-2*a-ae)*x[3]+a*x[2]+a*x[4]+ae*te;
    xx[4]=(-4*a-ae5-ah*u[1])*x[4]+a*x[3]+a*x[7]+a*x[5]+a*x[9]+ae5*te+ah*th*u[1];
    xx[5]=(-a-ae)*x[5]+a*x[4]+ae*te;
    xx[6]=(-a-ae)*x[6]+a*x[1]+ae*te;
    xx[7]=(-a-ae)*x[7]+a*x[4]+ae*te;
    xx[8]=(-a-ae)*x[8]+a*x[1]+ae*te;
    xx[9]=(-a-ae)*x[9]+a*x[4]+ae*te;
	};
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);
};
/* we integrate the growth bound by 0.5 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
  /* the ode for the growth bound */
  auto rhs =[](state_type& rr,  const state_type &r, const input_type &u) noexcept {
      rr[0]=0*r[0];
      rr[1]=0*r[1];
      rr[2]=0*r[2];
      rr[3]=0*r[3];
      rr[4]=0*r[4];
      rr[5]=0*r[5];
      rr[6]=0*r[6];
      rr[7]=0*r[7];
      rr[8]=0*r[8];
      rr[9]=0*r[9];
	};
  scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25}};
 /* lower bounds of the hyper-rectangle */
  state_type lb={{18,18,18,18,18,18,18,18,18,18}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{22,22,22,22,22,22,22,22,22,22}};
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{0,0}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{1,1}};
  /* grid node distance diameter */
  input_type i_eta={{0.5,0.5}};
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
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0] && 
        lb[1] <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= ub[1])
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

