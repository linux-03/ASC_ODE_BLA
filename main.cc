#include "equation_rbs.h"

int main()  {
  size_t steps = 1000;
  Vector<double> q_(12);
  q_.setConstant(0);
  q_(3) = 1;
  q_(11) = 1;
  q_(7) = 1;
  Vector<double> p_(6);
  p_.setConstant(0);
  q_(0) = 1;

  RigidBodySystem rbs;

  RigidBody rb(q_, p_);

  rbs.add(rb);

  Connector connect_a({0, 0, 0});
  Connector connect_b({0, 0, 0}, rb.Index());

  Beam bm(connect_a, connect_b);

  rbs.add(bm);

  simulate_rbs(rbs, 0.15/4, steps);


  return 0;
}