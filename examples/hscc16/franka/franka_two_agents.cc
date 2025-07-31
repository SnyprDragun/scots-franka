#include <array>
#include <iostream>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif


/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* sampling time */
const double tau = 0.3;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver ode_solver(sDIM,nint,tau);

/* we integrate the franka ode by 0.3 sec (the result is stored in x)  */
auto  franka_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the franka */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) {
      double alpha=std::atan(std::tan(u[1])/2.0);
      xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      xx[2] = u[0]*std::tan(u[1]);
  };
  ode_solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) {
    double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
    r[0] = r[0]+c*r[2]*0.3;
    r[1] = r[1]+c*r[2]*0.3;
};

/* forward declaration of the functions to setup the state space
 * input space and obstacles of the franka example */
scots::SymbolicSet frankaCreateStateSpace(Cudd &mgr);
scots::SymbolicSet frankaCreateInputSpace(Cudd &mgr);

void frankaCreateObstacles(scots::SymbolicSet &obs);


int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  scots::SymbolicSet ss=frankaCreateStateSpace(mgr);
  ss.writeToFile("franka_ss.bdd");

  /****************************************************************************/
  /* construct SymbolicSet for the obstacles */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet obs(ss);
  frankaCreateObstacles(obs);
  obs.writeToFile("franka_obst.bdd");

  /****************************************************************************/
  /* we define the target set */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet ts1(ss);
  /* define the target set as a symbolic set */
  double H[4*sDIM]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
  /* compute inner approximation of P={ x | H x<= h1 }  */
  double h1[4] = {-6,6.51,-0, .51};
  ts1.addPolytope(4,H,h1, scots::INNER);
  ts1.writeToFile("franka_target1.bdd");

  scots::SymbolicSet ts2(ss);
  double h2[4] = {-6,6.51,-4, 4.51};
  ts2.addPolytope(4,H,h2, scots::INNER);
  ts2.writeToFile("franka_target2.bdd");

  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  scots::SymbolicSet is=frankaCreateInputSpace(mgr);

  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  /* first create SymbolicSet of post variables 
   * by copying the SymbolicSet of the state space and assigning new BDD IDs */
  scots::SymbolicSet sspost(ss,1);
  /* instantiate the SymbolicModel */
  scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(franka_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

  /****************************************************************************/
  /* we continue with the controller synthesis */
  /****************************************************************************/
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp1(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T1 = ts1.getSymbolicSet();
  BDD O = obs.getSymbolicSet();
  tt.tic();
  /* compute controller */
  BDD C1=fp1.reachAvoid(T1,O,1);
  tt.toc();

  scots::FixedPoint fp2(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T2 = ts2.getSymbolicSet();
  tt.tic();
  /* compute controller */
  BDD C2=fp2.reachAvoid(T2,O,1);
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller1(ss,is);
  controller1.setSymbolicSet(C1);
  controller1.writeToFile("franka_controller1.bdd");

  scots::SymbolicSet controller2(ss,is);
  controller2.setSymbolicSet(C2);
  controller2.writeToFile("franka_controller2.bdd");

  return 1;
}

scots::SymbolicSet frankaCreateStateSpace(Cudd &mgr) {

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0,0,-M_PI-0.4};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={10,10,M_PI+0.4}; 
  /* grid node distance diameter */
  double eta[sDIM]={.2,.2,.2};   


  scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);

  /* add the grid points to the SymbolicSet ss */
  ss.addGridPoints();

  return ss;
}

void frankaCreateObstacles(scots::SymbolicSet &obs) {

  /* add the obstacles to the symbolic set */
  /* the obstacles are defined as polytopes */
  /* define H* x <= h */
  double H[4*sDIM]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

  double h2[4] = {-2.2,2.4,-0,5};
  obs.addPolytope(4,H,h2, scots::OUTER);
  double h5[4] = {-4.6 ,4.8,-1,10};
  obs.addPolytope(4,H,h5, scots::OUTER);

}

scots::SymbolicSet frankaCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={-1,-1};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={1,1}; 
  /* grid node distance diameter */
  double eta[sDIM]={.3,.3};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}