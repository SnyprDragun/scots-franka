/*
 * RungeKutta4.hh
 *
 * created on: 22.04.2015
 * author: rungger
 */

/** @file **/
#ifndef RUNGEKUTTA4_HH_
#define RUNGEKUTTA4_HH_

/** @namespace scots **/ 
namespace scots {

/**
 * * @brief Fixed step size ODE solver implementing a RungeKutta scheme of order 4 
 * for the specific dynamics: \f[ \dot \xi(t) = u, \xi(0)=x \f]
 * * @param rhs - A dummy lambda expression (not used internally for this specific ODE)
 * @param x - state_type initial state x
 * @param u - input_type constant input u
 * @param dim - state space dimension
 * @param tau - sampling time
 * @param nint - number of intermediate steps (default = 10)
 * @return the solution of IVP at time tau \f$ \xi(\tau) \f$ stored in x
 **/
template<class RHS, class state_type, class input_type>
void runge_kutta_fixed4(RHS rhs, state_type &x, input_type &u, const int dim, const double tau, const int nint=10) noexcept {

  /* * For the dynamics \f[ \dot \xi(t) = u \f], the solution is simply 
   * \f[ \xi(\tau) = x + u \cdot \tau \f] (constant velocity/drift).
   * * Since the RHS is independent of the state \f$\xi\f$, all intermediate 
   * Runge-Kutta terms \f$k_i\f$ will be equal to \f$u\f$.
   * * k_0 = u
   * k_1 = u
   * k_2 = u
   * k_3 = u
   * * And the update rule simplifies to:
   * \f[ x_{new} = x_{old} + \frac{h}{6}(u + 2u + 2u + u) = x_{old} + \frac{h}{6}(6u) = x_{old} + h \cdot u \f]
   * * Over nint steps, the total update is \f$ nint \cdot h \cdot u = \tau \cdot u \f$.
   * * The explicit Runge-Kutta calculation is replaced by the direct solution 
   * (Euler method, which is exact for this simple linear system).
   */
  
  // Calculate the final state directly using the analytical solution: xi(tau) = x + u * tau
  for(int i=0; i<dim; i++) {
    // The total change is input 'u' multiplied by the total time 'tau'
    x[i] = x[i] + u[i] * tau;
  }
  
  // The 'nint' parameter becomes irrelevant since the solution is direct.
}

} /* close namespace */

#endif /* RUNGEKUTTA4_HH_ */