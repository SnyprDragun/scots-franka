/*
 * synthesis.cc
 *
 *  created: Jun 2017
 *   author: Matthias Rungger
 */

/* 
 * this example demonstrates the usage of SCOTS together with slugs 
 * 
 * for the unicycle dynamics that we use to model the motion of the Khepera4
 * robot
 *
 * 1. compute the symbolic model 
 * 2. compute BDDs that represent atomic propositions
 * 3. use slugs to synthesize a controller
 * 4. simulate the closed loop in MATLAB 
 *
 */

#include <iostream>
#include <array>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* customized computation of symbolic model to account for periodicity in third
 * coordiante of the unicycle */
#include "SymbolicModel.hh"

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
/* sampling time */
const double tau = 0.3;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* we integrate the unicycle ode by tau sec (the result is stored in x)  */
auto  system_post = [](state_type &x, const input_type &u) {
  /* the ode describing the unicycle */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) -> void {
    xx[0] = u[0]*std::cos(x[2]);
    xx[1] = u[0]*std::sin(x[2]);
    xx[2] = u[1];
  };
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau);
};

auto radius_post = [](state_type &r, const state_type &, const input_type &u) {
  r[0] = r[0]+r[2]*std::abs(u[0])*tau;
  r[1] = r[1]+r[2]*std::abs(u[0])*tau;
};

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable();
  //mgr.AutodynDisable();

  /* lower bounds of the hyper rectangle */
  state_type s_lb={{0,0,-M_PI-0.4}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{5,5,M_PI+0.4}};
  /* grid node distance diameter */
  state_type s_eta={{.2,.2,.2}};
  /* construct SymbolicSet with the UniformGrid information for the state space
   * and BDD variable IDs for the pre */
  scots::SymbolicSet ss_pre = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
  /* construct SymbolicSet with the UniformGrid information for the state space
   * and BDD variable IDs for the post */
  scots::SymbolicSet ss_post = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
  /* assign names to BDD variables, used to interface with slugs */
  ss_pre.set_slugs_var_names(std::vector<std::string>{"x","y","z"});
  /* post variables are marked with ' at the end of each name */
  ss_post.set_slugs_var_names(std::vector<std::string>{"x'","y'","z'"});
 
  std::cout << "Unfiorm grid details:" << std::endl;
  ss_pre.print_info(1);
  ss_pre.print_info(1);

  /* create SymbolicSet for door bit (indicating whether the door is open or closed) */
  using M = std::array<double,1>;
  scots::SymbolicSet ss_door(mgr,1,M{{1}},M{{2}},M{{1}});
  ss_door.set_slugs_var_names(std::vector<std::string>{"D"});
  ss_door.print_info(1);



  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-2,-2}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{ 2, 2}};
  /* grid node distance diameter */
  input_type i_eta={{.2,.2}};
  /* create symbolic set */
  scots::SymbolicSet ss_input = scots::SymbolicSet(mgr,  input_dim, i_lb, i_ub, i_eta);
  /* names of the variables used to interface with slugs */
  ss_input.set_slugs_var_names(std::vector<std::string>{"u'","v'"});
  ss_input.print_info(1);

  /* names of the variables used to interface with slugs */
  scots::SymbolicSet ss_input_pre = scots::SymbolicSet(mgr,  input_dim, i_lb, i_ub, i_eta);
  ss_input_pre.set_slugs_var_names(std::vector<std::string>{"u'","v'"});
 
  scots::SymbolicSet  ss_controller(ss_pre,ss_input);
  write_to_file(mgr,ss_controller,mgr.bddOne(),"controller1");

  BDD TF;
  scots::SymbolicModel<state_type,input_type> sym_model(ss_pre,ss_input,ss_post);

  /* set up constraint functions with obtacles */
  double H[2][4] = {
    { 2.5, 2.75,   0, 2.5 },
    { 2.5, 2.75, 3.5,   5 },
  };
  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&H,&ss_pre](const abs_type& idx) {
    state_type x;
    ss_pre.itox(idx,x);
    double c1= ss_pre.get_eta()[0]/2.0+1e-10;
    double c2= ss_pre.get_eta()[1]/2.0+1e-10;
    for(size_t i=0; i<2; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0] <= (H[i][1]+c1) && 
          (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };
  /* compute BDD for the avoid set (returns the number of elements) */ 
  BDD A = ss_pre.ap_to_bdd(mgr,avoid);
  /* write ap to files avoid.scs/avoid.bdd */
  scots::write_to_file(mgr,ss_pre,A,"obstacles");

  scots::SymbolicSet set = scots::SymbolicSet(scots::SymbolicSet(ss_pre,ss_input),ss_post);

  std::cout << "Computing the transition function: " << std::endl;
  tt.tic();
  size_t no_trans;
  TF = sym_model.compute_gb(mgr,system_post,radius_post,avoid,no_trans);
  tt.toc();
  std::cout << "Number of transitions: " << no_trans << std::endl;
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;


  if(scots::write_to_file(mgr,set,TF,"tf"))
    std::cout << "success in writing\n";
  set.print_info(1);

  /* define first target set */
  auto target1 = [&ss_pre](const abs_type& idx) {
    state_type x;
    ss_pre.itox(idx,x);
    double r0 = ss_pre.get_eta()[0]/2.0;
    double r1 = ss_pre.get_eta()[1]/2.0;
    /* function returns 1 if cell associated with x is in target set  */
    if (4.5 <= (x[0]-r0) && (x[0]+r0) <= 5.001 && 
        4.5 <= (x[1]-r1) && (x[1]+r1) <= 5.001)
      return true;
    return false;
  };
  BDD T1 = ss_pre.ap_to_bdd(mgr,target1);
  /* write target to file */
  write_to_file(mgr,ss_pre,T1,"target1");

  /* define second target set */
  auto target2 = [&ss_pre](const abs_type& idx) {
    state_type x;
    ss_pre.itox(idx,x);
    double r0 = ss_pre.get_eta()[0]/2.0;
    double r1 = ss_pre.get_eta()[1]/2.0;
    /* function returns 1 if cell associated with x is in target set  */
    if (0 <= (x[0]-r0) && (x[0]+r0) <= .5 && 
        0 <= (x[1]-r1) && (x[1]+r1) <= .5)
      return true;
    return false;
  };
  BDD T2 = ss_pre.ap_to_bdd(mgr,target2);
  /* write target to file */
  write_to_file(mgr,ss_pre,T2,"target2");


  /* define atomic proposition for the door region */
  auto door = [&ss_pre](const abs_type& idx) {
    state_type x;
    ss_pre.itox(idx,x);
    double r0 = ss_pre.get_eta()[0]/2.0;
    double r1 = ss_pre.get_eta()[1]/2.0;
    /* function returns 1 if cell associated with x is in target set  */
    if (2.5 <= (x[0]+r0) && (x[0]-r0) <= 2.75 && 
        2.5 <= (x[1]+r1) && (x[1]-r1) <= 3.5)
      return true;
    return false;
  };
  BDD D = ss_pre.ap_to_bdd(mgr,door);
  /* write target to file */
  write_to_file(mgr,ss_pre,D,"door");




  /* 
   * we implement the fixed point algorithm  to see if the spec is realizable
   *
   * nu Z. ( mu X1.  pre(X1) | ( T1 & Z) ) & ( mu X2.  pre(X2) | ( T2 & Z) )
   *
   */

  /* init controller as empty */
  std::vector<BDD> C;
  C.push_back(mgr.bddZero());
  C.push_back(mgr.bddZero());



  std::vector<BDD> T;
  T.push_back(T1);
  T.push_back(T2);

  /* setup enforcable predecessor */
  scots::EnfPre enf_pre(mgr,TF,sym_model);
  tt.tic();

  /* outer fp*/
  BDD Z   = mgr.bddZero();
  BDD ZZ  = mgr.bddOne();

  /* inner fps*/
  std::vector<BDD> Y;
  std::vector<BDD> YY;
  for(int i=0; i<2; i++) {
    Y.push_back(mgr.bddOne());
    YY.push_back(mgr.bddZero());
  }

  /* helper */
  BDD U=ss_input.get_cube(mgr);

  size_t t,j=0;

  tt.tic();
  /* outer fixed point iteration */
  while(Z != ZZ) {
    std::cout << "\nOuter loop " << j++;
    Z = ZZ;
    BDD preZ = enf_pre(Z);
    ZZ = mgr.bddOne();
    for(int i=0; i<2; i++) {
      std::cout << "\nInner loop " << i;
      YY[i] = mgr.bddZero();
      C[i] = mgr.bddZero();
      /* inner fixed point iteration */
      t=0;
      while(Y[i] != YY[i]) {
        t++;
        Y[i]   = YY[i];
        YY[i]  = (T[i] & preZ) | (enf_pre(Y[i]));
        BDD N = YY[i] & (!(C[i].ExistAbstract(U)));
        C[i] = C[i] | N;
        /* print progress */
        scots::print_progress(t);
      }

      ZZ = ZZ & YY[i];
    }
  }
  tt.toc();

  /*checking on the sizes of the winning domains*/
  std::cout << "Winning domain size of rechability controller 1: " << ss_pre.get_size(mgr,C[0]) << "\n";
  std::cout << "Winning domain size of rechability controller 2: " << ss_pre.get_size(mgr,C[1]) << "\n";



  return 1;
}
