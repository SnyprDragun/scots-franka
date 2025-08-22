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
const int state_dim=8;
/* input space dim */
const int input_dim=4;
/* sampling time */
//const double tau = 10;

/*
 * data types for the elements of the state space 
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


/* parameters for radius calculation */
const double mu=std::sqrt(2);
/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
   /* the ode describing the dcdc converter */
   state_type z = x;
   /* the ode describing the system */
    const double v = 120 ; 
    const double l = 1 ;
    const double T = 0.002777;
    //const double q = 0.25;
    x[0] = (0.5-T*v/l)*z[0] +T*v/l*x[7]+ 25*u[0];
    x[1] = T*v/l*z[0] + (0.65-T*v/l) *z[1];
    x[2] = (0.5-T*v/l)*z[2] +T*v/l*x[1]+ 25*u[1];
    x[3] = T*v/l*z[2] + (0.65-T*v/l) *z[3];
    x[4] = (0.5-T*v/l)*z[4] +T*v/l*x[3]+ 25*u[2];
    x[5] = T*v/l*z[4] + (0.65-T*v/l) *z[5];
    x[6] = (0.5-T*v/l)*z[6] +T*v/l*x[5]+ 25*u[3];
    x[7] = T*v/l*z[6] + (0.65-T*v/l) *z[7];
};
/* we integrate the growth bound by 0.5 sec (the result is stored in r)  */

auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
	r[0]=0*r[0];
      	r[1]=0*r[1];
      	r[2]=0*r[2];
     	r[3]=0*r[3];
     	r[4]=0*r[4];
	r[5]=0*r[5];
     	r[6]=0*r[6];
	r[7]=0*r[7];
  /* the ode for the growth bound */
//  auto rhs =[](state_type& rr,  const state_type &r, const input_type &u) noexcept {
//      rr[0]=0*r[0];
//      rr[1]=0*r[1];
//      rr[2]=0*r[2];
//      rr[3]=0*r[3];
//      rr[4]=0*r[4];
//	};
//  scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{3,3,3,3,3,3,3,3}};
 /* lower bounds of the hyper-rectangle */
  state_type lb={{2,2,2,2,2,2,2,2}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{30,30,30,30,30,30,30,30}};
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{0,0,0,0}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{1,1,1,1}};
  /* grid node distance diameter */
  input_type i_eta={{1,1,1,1}};
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

