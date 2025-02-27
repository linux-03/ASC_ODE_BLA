#include "equation_rbs.h"

int main()  {
  size_t steps = 10000;
  Vector<double> q_(12);
  q_.setConstant(0);
  q_(3) = 1;
  q_(11) = 1;
  q_(7) = 1;
  Vector<double> p_(6);
  p_.setConstant(0);
  q_(1) = 0.1;
  q_(2) = 1;

  RigidBodySystem rbs;

  RigidBody rb(q_, p_);

  rbs.add(rb);

  Connector connect_b({0, 0, 0});
  Connector connect_a({0, -0.5, 0}, rb.Index());

  Beam bm(connect_a, connect_b);

  rbs.add(bm);
  std::cout << bm.RelPosAToB() << " " << bm.RelPosBToA() << std::endl;
  std::cout << bm.Length();
  //for (size_t i = 0; i < steps; i++) {
  
  simulate_rbs(rbs, 0.004, steps, [](int i, double err, VectorView<double> q) {
                    std::cout << std::scientific << "Body1 newton-iteration: " << i << " newton-error: " << err << std::endl
                      <<"\t"<< "Translation = " << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                      }
                  );
  //}
  
  return 0;
}

