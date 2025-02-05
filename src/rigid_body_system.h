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
      Matrix<T> G(12, this->Bodies(body_index).Beams().size());
      G.setConstant(0);
      
      for (size_t bm_index: this->Bodies(body_index).Beams())
      {
        Beam bm = Beams(bm_index);
        size_t index_a = bm.BodyIndexA();
        size_t index_b = bm.BodyIndexB();
        G.Col(bm_index).segment(0, 12) += G_single_body_beam( bm, body_index, x.segment(index_a*dim_per_body(), 12), x.segment(index_b*dim_per_body(), 12));
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
    Vector<T> G_single_body_beam(Beam& bm, size_t bd_index, VectorView<T> q_a, VectorView<T> q_b) {
    
      Vector<T> res(dim_per_motion());

      Vec<3, T> pos1 = bm.ConnectorA().RefPosition();
      Vec<3, T> pos2 = bm.ConnectorB().RefPosition();

      // if a is a fix
      if (bm.ConnectorA().Fix())  {
        for (size_t i = 0; i < 3; i++)  {
          T row = 2*(pos1(i) - q_b.Range(3 + i*3, 3 + i*3 + 3)*pos2 - q_b(i));

          res(i) = (-1)*row;

          for (size_t j = 0; j < 3; j ++)
          {
            res(3 + i*3 + j) =  (-1)*pos2(j)*row;
          }
        }
      }
      //if b is a fix
      else if (bm.ConnectorB().Fix())  {
        for (size_t i = 0; i < 3; i++)  {
          T row = 2*(q_a.Range(3 + i*3, 3 + i*3 + 3)*pos1 + q_a(i) - pos2(i));

          res(i) = row;
          res(dim_per_motion() + i) = 0;

          for (size_t j = 0; j < 3; j ++) {
            res(3 + i*3 + j) = pos1(j)*row;
          }
        }
      }
      else {
        // check for index to take derivative after - case that it is a
        if (bd_index == bm.ConnectorA().BodyIndex())
        {
          for (size_t i = 0; i < 3; i++)  {
            T row = 2*(q_a.Range(3 + i*3, 3 + i*3 + 3)*pos1 + q_a(i) - q_b.Range(3 + i*3, 3 + i*3 + 3)*pos2 - q_b(i));

            res(i) = row;

            for (size_t j = 0; j < 3; j ++) {
              res(3 + i*3 + j) = pos1(j)*row;
            }
          }
        }
        // case that it is b
        else
        {
          for (size_t i = 0; i < 3; i++)  {
            T row = 2*(q_a.Range(3 + i*3, 3 + i*3 + 3)*pos1 + q_a(i) - q_b.Range(3 + i*3, 3 + i*3 + 3)*pos2 - q_b(i));

            res(i) = (-1)*row;

            for (size_t j = 0; j < 3; j ++) {
              res(3 + i*3 + j) = (-1)*pos2(j)*row;
            }
          }
        }
      }
      return res;
    }


    /**
    * @brief Computes the VelocityConstraint for a specific beam
    *
    * @param x The full position vector containing each position and momentum of each body
    * @param beam_index The index of the beam to calculate for
    * @return T value of the constraint
    */
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

    G_a += G_single_body_beam(bm, bd_index_a, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12));
    G_b += G_single_body_beam(bm, bd_index_b, x.segment(bd_index_a*dim_per_body(), 12), x.segment(bd_index_b*dim_per_body(), 12));
    

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