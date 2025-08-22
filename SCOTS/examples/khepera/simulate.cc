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

/* sampling time */
const double tau = 0.3;

using state_type = std::array<double,5>;
using input_type = std::vector<double>;

/* we integrate the unicycle ode by tau sec (the result is stored in x)  */
auto  system_post = [](state_type &x, const input_type &u) {
  /* the ode describing the unicycle */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) -> void {
    xx[0] = u[0]*std::cos(x[2]);
    xx[1] = u[0]*std::sin(x[2]);
    xx[2] = u[1];
  };
  scots::runge_kutta_fixed4(rhs,x,u,3,tau);
};

/* define first target set */
auto target1 = [](const state_type& x) {
  /* function returns 1 if cell associated with x is in target set  */
  if (4.5 <= (x[0]) && (x[0]) <= 5.001 && 4.5 <= (x[1]) && (x[1]) <= 5.001)
    return true;
  return false;
};

/* define second target set */
auto target2 = [](const state_type& x) {
  /* function returns 1 if cell associated with x is in target set  */
  if (0 <= (x[0]) && (x[0]) <= .5 && 0 <= (x[1]) && (x[1]) <= .5)
    return true;
  return false;
};


/* define atomic proposition for the door region */
auto door_region = [](const state_type& x) {
  /* function returns 1 if cell associated with x is in target set  */
  if (2.2 <= (x[0]) && (x[0]) <= 3 && 2.25 <= (x[1]) && (x[1]) <= 3.7)
    return true;
  return false;
};

int main() {

  /* Cudd manager */
  Cudd manager;

  /* read controller from file */
  BDD C;
  scots::SymbolicSet con;
  if(!read_from_file(manager,con,C,"controller",'A')) {
    std::cout << "Could not read controller from controller.scs\n";
    return 0;
  }


  std::cout << "\nSimulation:\n " << std::endl;

  //std::array<double,7> xx;
  //con.itox(0,xx);

  //BDD zero = con.id_to_bdd(0);

  //
  //BDD cube =  manager.bddVar(39) 
  //          & manager.bddVar(37) 
  //          & manager.bddVar(35) 
  //          & manager.bddVar(45) 
  //          & manager.bddVar(43) 
  //          & manager.bddVar(41) ;
  //                    
  //cube.PrintMinterm();

  //BDD tmp = zero.ExistAbstract(cube); 

  //tmp = tmp & C;


  //tmp.PrintMinterm();

  con.print_info(1);

 // for(size_t i=0; i<xx.size(); i++) 
 //   std::cout << xx[i] <<  " ";
 //   std::cout << "\n";

 state_type x={{1.2,1.2, -3.4,1,1}};
 std::vector<int> dom {0,1,2,5,6};

  std::srand(std::time(0));

  size_t t=200;
  while(--t) {
    /* door behavior */
    int random_variable = std::rand();
    x[3]  = random_variable%2+1;
    if(door_region(x))
      x[3]=1;
    std::cout << x[3] <<  " ";

    /* returns a std vector with the valid control inputs */
    if(x[2]<-M_PI)
      x[2]=x[2]+2*M_PI;
    if(x[2]>M_PI)
      x[2]=x[2]-2*M_PI;
    std::cout << x[0] <<  " "  << x[1] << " " << x[2] << " ";
    auto u = con.restriction(manager,C,x,dom);
    if(u.size()) {
      std::cout << u[0] <<  " "  << u[1] << "\n";
      system_post(x,u);
    } else
      break;

    if(target1(x)) {
      //std::cout << "Arrived at one " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      x[4]=2;
    }
    if(target2(x)) {
      //std::cout << "Arrived at two " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      x[4]=1;
    }
  }

  if(t)
    std::cout << "something wrong " << t << "\n";

  return 1;
}
