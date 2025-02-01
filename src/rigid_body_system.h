#ifndef RIGID_BODY_SYSTEM_H
#define RIGID_BODY_SYSTEM_H

#include "rigid_body.h"
#include "Beam.h"
#include "Spring.h"
#include "Connector.h"

using namespace std;

class RigidBodySystem {
    size_t num_beams_;
    size_t num_springs_;
    std::vector<std::reference_wrapper<RigidBody>> bodies_;
    std::vector<std::reference_wrapper<Beam>> beams_;
    std::vector<std::reference_wrapper<Spring>> springs_;
    Vector<double> gravity_ = {0, 9.81, 0};

  public:
    RigidBodySystem();

    template<typename T>
    T Potential(const VectorView<T> x);
    template<typename T>
    Vector<T> PotentialGradient(const VectorView<T> x, size_t body_index)
    {
      Vector<T> pot_grad(12);
      pot_grad.setConstant(0);

      pot_grad.segment(0, 3) = gravity_;

      return pot_grad;
    }
    template<typename T>
    Matrix<T> JacobianConstraint(const VectorView<T> x, size_t body_index)
    {
      Matrix<T> G(12, this->Bodies(body_index).Beams().size());
      G.setConstant(0);  
      
      for (size_t bm_index: this->Bodies(body_index).Beams())
      {
        G.Col(bm_index).segment(0, 3) += 2*x.segment(body_index*dim_per_body(), 3);
      }
      
      
      return G;
    }
    template<typename T>
    T Constraint(const VectorView<T> x, size_t beam_index)
    {
      Beam& bm = this->beams_[beam_index];
      size_t bd_index_a = bm.BodyIndexA();
      size_t bd_index_b = bm.BodyIndexB();
      Vector<T> q = bm.PositionA( x.segment(bd_index_a*dim_per_body(), 12) ) - bm.PositionB( x.segment(bd_index_b*dim_per_body(), 12) );
      T g; // One constraint for a simple pendulum
      g = q.squaredNorm() - bm.Length(); // Length of beam constraint
      return g;
    }

    template<typename T>
    T Velocity_Constraint(const VectorView<T> x, size_t beam_index)
    {
    Vector<T> G_a(12); // Jacobian for one constraint and three translational coordinates
    Vector<T> G_b(12);

    G_a.setConstant(0);
    G_b.setConstant(0);

    Beam& bm = this->beams_[beam_index];

    size_t bd_index_a = bm.BodyIndexA();
    size_t bd_index_b = bm.BodyIndexB();

    if (!bm.ConnectorA().Fix()) {
      G_a.segment(0, 3) = 2*x.segment(bd_index_a * dim_per_body(), 3);
    }
    if (!bm.ConnectorB().Fix()) {
      G_b.segment(0, 3) = 2*x.segment(bd_index_b * dim_per_body(), 3);
    }
    else {
      G_a.segment(0, 3) = 2*(x.segment(bd_index_a * dim_per_body(), 3) - x.segment(bd_index_b * dim_per_body(), 3));
    }

    T c_a = G_a.segment(0, 3)*x.segment(bd_index_a * dim_per_body() + 18, 3) + 
                  G_a.segment(3, 9)*ToVector(vectorToSkewSymmetric( x.segment(bd_index_a * dim_per_body() + 18 + 3, 3) ) );

    T c_b = G_b.segment(0, 3)*x.segment(bd_index_b * dim_per_body() + 18, 3) + 
                  G_b.segment(3, 9)*ToVector(vectorToSkewSymmetric( x.segment(bd_index_b * dim_per_body() + 18 + 3, 3) ) );

    return c_a + c_b;
  }

    size_t& NumBeams();
    size_t& NumSprings();
    size_t NumBodies();
    
    RigidBody& Bodies(size_t i);
    Beam& Beams(size_t i);
    Spring& Springs(size_t i);

    void SaveState(const VectorView<double> x);
    void ManageConstraints(const VectorView<double> x);
    Vector<double> ExpandState();

    void add(RigidBody& rb);
    void add(Beam& bm);
    void add(Spring& spr);
};

#endif