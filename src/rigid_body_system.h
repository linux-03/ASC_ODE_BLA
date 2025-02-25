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

    /**
    * @brief Computes the Gradient of the potential Energy for one body
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param body_index The index of the body to get the gradient fo
    * @return 12 Dimensional gradient
    */
    template<typename T>
    Vector<T> PotentialGradient(const VectorView<T> x, size_t body_index)
    {
      Vector<T> pot_grad(12);
      pot_grad.setConstant(0);

      pot_grad.segment(0, 3) = gravity_;

      return pot_grad;
    }

    /**
    * @brief Computes the Jacobien of the function of all constraint derived after a specific body
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param body_index The index of the body to derive after
    * @return Jacobina Matrix with heigth according to the number of beams and width 12
    */
    template<typename T>
    Matrix<T> JacobianConstraint(const VectorView<T> x, size_t body_index)
    {
      Matrix<T> G(12, this->Bodies(body_index).Beams().size() * dim_per_beam()/2);
      G.setConstant(0);
      
      for (size_t bm_index: this->Bodies(body_index).Beams())
      {
        Beam bm = Beams(bm_index);
        size_t index_a = bm.BodyIndexA();
        size_t index_b = bm.BodyIndexB();
        G.Block(0, bm_index*dim_per_beam()/2, 12, 3) += G_single_body_beam_n( bm, body_index, x.segment(index_a*dim_per_body(), 12), x.segment(index_b*dim_per_body(), 12));
        //std::cout << G << std::endl;
      }

      return G;
    }

    /**
    * @brief Computes the length constraint of a beam
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param beam_index The index of the beam to have the constraint calculated
    * @return Value T of the length difference
    */
    template<typename T>
    Vector<T> Constraint(const VectorView<T> x, size_t beam_index)
    {
      Beam& bm = this->beams_[beam_index];
      size_t bd_index_a = bm.BodyIndexA();
      size_t bd_index_b = bm.BodyIndexB();
      Vector<T> pos1 =  bm.PositionB( x.segment(bd_index_b*dim_per_body(), 12) );
      Vector<T> pos1_rel = bm.RelPosBToA();
      Vector<T> pos2 = x.segment(bd_index_a*dim_per_body(), 3) + ToMatrix(x.segment(bd_index_a*dim_per_body() + 3, 9))*pos1_rel;
      Vector<T> q = pos1 - pos2;

      
      return q;
    }

    /**
    * @brief Computes the derivative of a beam length constrain with respect to a specific body.
    *
    * @param bm The Beam of which you want to have the derivative of the constraint
    * @param bd_index The index of the body to derive after
    * @param q_a The 12-dimesnional vector representing the body position of body A of the beam
    * @param q_b The 12-dimesnional vector representing the body position of body B of
    * @return 12-dimesional gradient
    */
    template<typename T>
    Matrix<T> G_single_body_beam_n(Beam& bm, size_t bd_index, VectorView<T> q_a, VectorView<T> q_b) {
    
      Vector<T> res(dim_per_motion());
      Vector<AutoDiff<12, T>> q_a_diff = q_a;
      Vector<AutoDiff<12, T>> q_b_diff = q_b;

      if ((bd_index == bm.BodyIndexA()) && !(bm.ConnectorA().Fix())) {
        for(size_t i =0 ; i < 12; i++){
          q_a_diff(i).DValue(i) = 1;
        }
      } else {
        for(size_t i=0; i < 12; i++){
          q_b_diff(i).DValue(i) = 1;
        }
      }

      
      size_t bd_index_a = bm.BodyIndexA();
      size_t bd_index_b = bm.BodyIndexB();
      Vector<AutoDiff<12, T>> pos1 =  bm.PositionB( q_b_diff );
      Vector<AutoDiff<12, T>> pos1_rel = bm.RelPosBToA();
      Vector<AutoDiff<12, T>> pos2 = q_a_diff.segment(0, 3) + ToMatrix(q_a_diff.segment(3, 9))*pos1_rel;
      Vector<AutoDiff<12, T>> q = pos1 - pos2;

      //AutoDiff<12, T> g = q.squaredNorm();
      //Vector<AutoDiff<12, T>> q = q_a_diff.segment(0, 3) + ToMatrix(q_a_diff.segment(3, 9))*pos1 - q_b_diff.segment(0, 3) - ToMatrix(q_b_diff.segment(3, 9))*pos2;
      //Vector<AutoDiff<12, T>> q = bm.PositionA(q_a_diff) - bm.PositionB(q_b_diff);
      
      //AutoDiff<12, T> g = q*q - bm.Length()*bm.Length(); // Length of beam constraint
      Matrix<T> g_vec(12, 3);
      for (size_t j = 0; j < 3; j++) {
        for(size_t i = 0 ; i < 12; i++){
          g_vec(i, j) = q(j).DValue(i);
        }
      }
      
      return g_vec;
    }


    /**
    * @brief Computes the VelocityConstraint for a specific beam
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param beam_index The index of the beam to calculate for
    * @return T value of the constraint
    */
    template<typename T>
    Vector<T> Velocity_Constraint(const VectorView<T> x, size_t beam_index)
    {
      Matrix<T> G_a(3, 12); // Jacobian for one constraint and three translational coordinates
      Matrix<T> G_b(3, 12);

      G_a.setConstant(0);
      G_b.setConstant(0);

      Beam& bm = this->beams_[beam_index];

      size_t bd_index_a = bm.BodyIndexA();
      size_t bd_index_b = bm.BodyIndexB();

      G_a += G_single_body_beam_n(bm, bd_index_a, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12)).transpose();
      G_b += G_single_body_beam_n(bm, bd_index_b, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12)).transpose();
      
      Matrix<T> R_a = ToMatrix(x.segment(bd_index_a * dim_per_body() + 3, 9)) * vectorToSkewSymmetric( x.segment(bd_index_a * dim_per_body() + 18 + 3, 3) ) ;
      Matrix<T> R_b = ToMatrix(x.segment(bd_index_b * dim_per_body() + 3, 9)) * vectorToSkewSymmetric( x.segment(bd_index_b * dim_per_body() + 18 + 3, 3) ) ;
      
      
      Vector<T> c_a = G_a.transpose().Block(0, 0, 3, 3)*x.segment(bd_index_a * dim_per_body() + 18, 3) + 
                    G_a.Block(0, 3, 3, 9) * ToVector(R_a);

      Vector<T> c_b = G_b.transpose().Block(0, 0, 3, 3)*x.segment(bd_index_b * dim_per_body() + 18, 3) + 
                    G_b.Block(0, 3, 3, 9) * ToVector(R_b);

      return c_a + c_b;
    }

    size_t& NumBeams();
    size_t& NumSprings();
    size_t NumBodies();
    
    RigidBody& Bodies(size_t i);
    Beam& Beams(size_t i);
    std::vector<Beam> Beams();
    std::vector<RigidBody> Bodies();
    Spring& Springs(size_t i);

    Vector<double> connectorPosition(Connector c);

    void SaveState(const VectorView<double> x);
    void ManageConstraints(const VectorView<double> x);
    Vector<double> ExpandState();

    void add(RigidBody& rb);
    void add(Beam& bm);
    void add(Spring& spr);
};

#endif