/*
 * simulate.cc
 *
 *  created: Oct 2016
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
 */

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* stop the while loop for some milliseconds */
#include <chrono> 
#include <thread> 


/* sampling time */
const double tau = 0.3;

/* data types for the state space elements and input space
 * elements used in uniform grid and ode solvers */
using state_type = std::array<double,3>;
using input_type = std::vector<double>;

/* we integrate the vehicle ode by tau sec (the result is stored in x)  */
auto  vehicle_post = [](state_type& x, const input_type& u) {
  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type& x, const input_type& u) {
    double alpha=std::atan(std::tan(u[1])/2.0);
    xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
    xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
    xx[2] = u[0]*std::tan(u[1]);
  };
  /* simulate (use 10 intermediate steps in the ode solver) */
  scots::runge_kutta_fixed4(rhs,x,u,3,tau,10);
  if(x[2]>M_PI)
    x[2]=x[2]-2*M_PI;
  if(x[2]<-M_PI)
    x[2]=x[2]+2*M_PI;

};

int main() {

  /* Cudd manager */
  Cudd manager;

  /* define function to check if we are in target 1 */
  auto T1 = [](const state_type& x) {
    if (9 <= x[0] && x[0] <= 9.5 && 0 <= x[1] && x[1] <= 0.5)
      return true;
    return false;
  };

  /* define function to check if we are in target 1 */
  auto T2 = [](const state_type& x) {
    if (0 <= x[0] && x[0] <= 0.5 && 0 <= x[1] && x[1] <= 0.5)
      return true;
    return false;
  };

  /* read controller one from file */
  BDD C1;
  scots::SymbolicSet con1;
  if(!read_from_file(manager,con1,C1,"controller1")) {
    std::cout << "Could not read controller from controller1\n";
    return 0;
  }

  /* read controller two from file */
  BDD C2;
  scots::SymbolicSet con2;
  if(!read_from_file(manager,con2,C2,"controller2")) {
    std::cout << "Could not read controller from controller2\n";
    return 0;
  }  

  std::cout << "\nSimulation:\n " << std::endl;
  
  /* current target */
  int target=1;

  input_type u;

  state_type x={{0.6, 0.6, 0}};
  while(1) {
    /* returns a std vector with the valid control inputs */
    if(target==1)
      u = con1.restriction(manager,C1,x);
    else
      u = con2.restriction(manager,C2,x);

    std::cout << x[0] <<  " "  << x[1] << " " << x[2] << "\n";
    //std::cout << u[0] <<  " "  << u[1] << "\n";
    vehicle_post(x,u);
    if(T1(x) && (target==1)) {
      std::cout << "Arrived at target one: " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      std::this_thread::sleep_for(std::chrono::seconds{1});
      target=2;
    }
    if(T2(x) && (target==2)) {
      std::cout << "Arrived at target two: " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      std::this_thread::sleep_for(std::chrono::seconds{1});
      target=1;
    }
  }

  return 1;
}
