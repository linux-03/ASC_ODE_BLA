#ifndef EQUATION_RB_CC
#define EQUATION_RB_CC

#include "equation_rb.h"
#include "Newton.h"
#include "autodiff.h"

namespace ASC_ode {

  RigidBodyEquation::RigidBodyEquation(RigidBody& rb, RigidBodySystem& rbs, double step_size): rb_(rb), rbs_(rbs), h_(step_size) {};

  size_t RigidBodyEquation::DimX() const {
    return this->dimX_;
  }

  size_t RigidBodyEquation::DimF() const {
    return this->dimF_;
  }

  template<typename T>
  Vector<T> RigidBodyEquation::FuncConstraint(const VectorView<T> x, bool old) const
  {
    Vector<T> res(12);
    res.setConstant(0);

    if (old)
    {
      for (size_t i = 0; i < rb_.Beams().size(); i++)
      {
        size_t bm_index = rb_.Beams()[i];
        //std::cout << rb_.Constraints().col(i);
        res +=  x(dim_per_body()*rbs_.NumBodies() + dim_per_beam()*bm_index) * rb_.Constraints().Col(i);
      }
    } else  {
      Matrix<T> constr = rbs_.JacobianConstraint(x, rb_.Index());
      for (size_t i = 0; i < rb_.Beams().size(); i++)
      {
        size_t bm_index = rb_.Beams()[i];
        res += x(dim_per_body()*rbs_.NumBodies() + dim_per_beam()*bm_index + 1) * constr.Col(i);
      }
    }

    return res;
  }

  void RigidBodyEquation::SetRotGradien(VectorView<double> x, VectorView<AutoDiff<24, double>> x_diff) const {

    if (x.norm() == 0) {
      //Gradinet for var a
      x_diff(5).DValue(3) = -1;
      x_diff(7).DValue(3) = 1;
      // Gradient for var b
      x_diff(2).DValue(4) = 1;
      x_diff(6).DValue(4) = -1;
      // Gradient for var c
      x_diff(1).DValue(5) = -1;
      x_diff(3).DValue(5) = 1;
    } else {
      Vector<AutoDiff<24, double>> axis_diff = x;

      AutoDiff<24, double> angle = x.norm();

      Matrix<AutoDiff<24, double>> H = vectorToSkewSymmetric(axis_diff);
  
      Matrix<AutoDiff<24, double>> res = Diagonal(3, 1) +  (sin(angle)/angle) * H + (1 - cos(angle))/(angle*angle) * H*H;

      x_diff = ToVector(res);
    }
  }


  void RigidBodyEquation::Evaluate(VectorView<double> x, VectorView<double> f) const
  {
    f.setConstant(0);

    Matrix<double> Rmean = 0.5*(ToMatrix(x.segment(3, 9)) + rb_.q());
    Matrix<double> Rnew = ToMatrix(x.segment(3,9));

    Vector<double> force_old = rb_.Force();
    Vector<double> vel_con_old = FuncConstraint(x, true);

    Vector<double> force_new = rbs_.PotentialGradient(x, rb_.Index());
    Vector<double> vel_con_new = FuncConstraint(x, false);

    // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1
    //std::cout << x << std::endl << std::endl;
    f.setConstant(0);
    
    f.segment(6, 3) = x.segment(0, 3) - rb_.q_trans() - h_*x.segment(24, 3); // q_trans_n - q_trans - h_*p
    f.segment(9, 3) = skewSymmetricToVector( ToMatrix(x.segment(3, 9)) - rb_.q() - h_*vectorToSkewSymmetric(x.segment(27, 3)));


    f.segment(18, 6) = x.segment(12, 6) - x.segment(24, 6);


    f.segment(0, 3) = x.segment(12, 3) - rb_.p_trans() + (h_/2)*(force_old.segment(0, 3) - vel_con_old.segment(0, 3));
    f.segment(3, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(x.segment(15, 3)) - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) + (h_/2)*(ToMatrix(force_old.segment(3, 9)) - ToMatrix(vel_con_old.segment(3, 9)))));


    f.segment(12, 3) = x.segment(18, 3) - x.segment(12, 3) + (h_/2)*(force_new.segment(0, 3) - vel_con_new.segment(0, 3)); // v - p + h/2 * f - f_v' 
    f.segment(15, 3) = skewSymmetricToVector(Rnew.transpose()*(Rnew*vectorToSkewSymmetric(x.segment(21, 3)) - Rmean*vectorToSkewSymmetric(x.segment(15, 3)) + (h_/2)*(ToMatrix(force_new.segment(3, 9)) - ToMatrix(vel_con_new.segment(3, 9)))));
  }

  void RigidBodyEquation::EvaluateDeriv(VectorView<double> x, MatrixView<double> df)  const {
    df.setConstant(0);

    Vector<AutoDiff<24, double>> x_diff = x;
    size_t prev_index = 0;
    Vector<AutoDiff<24, double>> f_diff(24);

    for (size_t bd_index = 0; bd_index < rbs_.NumBodies(); bd_index++) {

      x_diff(dim_per_body()*prev_index).DValue(0) = 0;
    x_diff(dim_per_body()*bd_index).DValue(0) = 1;

    x_diff(dim_per_body()*prev_index + 1).DValue(1) = 0;
    x_diff(dim_per_body()*bd_index + 1).DValue(1) = 1;

    x_diff(dim_per_body()*prev_index + 2).DValue(2) = 0;
    x_diff(dim_per_body()*bd_index + 2).DValue(2) = 1;
    

    for (size_t j = 0; j < 18; j++) {
      x_diff(dim_per_body()*prev_index + 12 + j).DValue(6 + j) = 0;
      x_diff(dim_per_body()*bd_index + 12 + j).DValue(6 + j) = 1;
    }
      prev_index = bd_index;

      SetRotGradien(x.segment(dim_per_body() * bd_index + 30, 3), x_diff.segment(dim_per_body() * bd_index + 3, 9));
      /*
      std::cout << x_diff(0) << std::endl;
      std::cout << x_diff(1) << std::endl;
      std::cout << x_diff(2) << std::endl;
      std::cout << x_diff(3) << std::endl;
      std::cout << x_diff(4) << std::endl;
      std::cout << x_diff(5) << std::endl;
      std::cout << x_diff(6) << std::endl;
      std::cout << x_diff(7) << std::endl;
      std::cout << x_diff(8) << std::endl;
      std::cout << x_diff(9) << std::endl;
      std::cout << x_diff(10) << std::endl;
      std::cout << x_diff(11) << std::endl;
      */
      Matrix<AutoDiff<24, double>> Rmean = 0.5*(ToMatrix(x_diff.segment(3, 9)) + rb_.q());
      Matrix<AutoDiff<24, double>> Rnew = ToMatrix(x_diff.segment(3,9));

      Vector<AutoDiff<24, double>> force_old = rb_.Force();
      Vector<AutoDiff<24, double>> vel_con_old = FuncConstraint(x_diff, true);

      Vector<AutoDiff<24, double>> force_new = rbs_.PotentialGradient(x_diff, rb_.Index());
      Vector<AutoDiff<24, double>> vel_con_new = FuncConstraint(x_diff, false);

      // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1
      //std::cout << x << std::endl << std::endl;
      
      
      f_diff.segment(0, 3) = x_diff.segment(12, 3) - rb_.p_trans() + (h_/2)*(force_old.segment(0, 3) - vel_con_old.segment(0, 3));

      f_diff.segment(3, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(x_diff.segment(15, 3)) - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) + (h_/2)*(ToMatrix(force_old.segment(3, 9)) - ToMatrix(vel_con_old.segment(3, 9)))));

      f_diff.segment(6, 3) = x_diff.segment(0, 3) - rb_.q_trans() - h_*x_diff.segment(24, 3);

      f_diff.segment(9, 3) = skewSymmetricToVector( ToMatrix(x_diff.segment(3, 9)) - rb_.q() - h_*vectorToSkewSymmetric(x_diff.segment(27, 3)));

      f_diff.segment(12, 3) = x_diff.segment(18, 3) - x_diff.segment(12, 3) + (h_/2)*(force_new.segment(0, 3) - vel_con_new.segment(0, 3));

      f_diff.segment(15, 3) = skewSymmetricToVector(Rnew.transpose()*(Rnew*vectorToSkewSymmetric(x_diff.segment(21, 3)) - Rmean*vectorToSkewSymmetric(x_diff.segment(15, 3)) + (h_/2)*(ToMatrix(force_new.segment(3, 9)) - ToMatrix(vel_con_new.segment(3, 9)))));

      f_diff.segment(18, 6) = x_diff.segment(12, 6) - x_diff.segment(24, 6);

      for (size_t j = 0; j < 24; j++) {
        df.Row(rb_.Index()*dim_per_body_rbs() + j).segment(dim_per_body_rbs()*bd_index, dim_per_body_rbs()) = f_diff(j); 
      }
    }

    prev_index = 0;

    Vector<AutoDiff<2, double>> x_bm_diff = x;
    Vector<AutoDiff<2, double>> f_bm_diff(24);
    
    for (size_t bm_index = 0; bm_index < rbs_.NumBeams(); bm_index++) {
      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*prev_index).DValue(0) = 0;
      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*prev_index + 1).DValue(1) = 0;

      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*bm_index).DValue(0) = 1;
      x_bm_diff(dim_per_body()*rbs_.NumBodies() + 2*bm_index + 1).DValue(1) = 1;
      prev_index = bm_index;


      Matrix<AutoDiff<2, double>> Rmean = 0.5*(ToMatrix(x_bm_diff.segment(3, 9)) + rb_.q());
      Matrix<AutoDiff<2, double>> Rnew = ToMatrix(x_bm_diff.segment(3,9));

      Vector<AutoDiff<2, double>> force_old = rb_.Force();
      Vector<AutoDiff<2, double>> vel_con_old = FuncConstraint(x_bm_diff, true);

      Vector<AutoDiff<2, double>> force_new = rbs_.PotentialGradient(x_bm_diff, rb_.Index());
      Vector<AutoDiff<2, double>> vel_con_new = FuncConstraint(x_bm_diff, false);

      // q_trans 3 q_rot 9 p_half 6 p 6 v 6 lambda 1 mu 1
      //std::cout << x << std::endl << std::endl;
      
      
      f_bm_diff.segment(0, 3) = x_bm_diff.segment(12, 3) - rb_.p_trans() + (h_/2)*(force_old.segment(0, 3) - vel_con_old.segment(0, 3));

      f_bm_diff.segment(3, 3) = skewSymmetricToVector( rb_.q().transpose() * ( Rmean * vectorToSkewSymmetric(x_bm_diff.segment(15, 3)) - rb_.q()* vectorToSkewSymmetric(rb_.p_skew()) + (h_/2)*(ToMatrix(force_old.segment(3, 9)) - ToMatrix(vel_con_old.segment(3, 9)))));
      
      f_bm_diff.segment(6, 3) = x_bm_diff.segment(0, 3) - rb_.q_trans() - h_*x_bm_diff.segment(24, 3);

      f_bm_diff.segment(9, 3) = skewSymmetricToVector( ToMatrix(x_bm_diff.segment(3, 9)) - rb_.q() - h_*vectorToSkewSymmetric(x_bm_diff.segment(27, 3)));

      f_bm_diff.segment(12, 3) = x_bm_diff.segment(18, 3) - x_bm_diff.segment(12, 3) + (h_/2)*(force_new.segment(0, 3) - vel_con_new.segment(0, 3));

      f_bm_diff.segment(15, 3) = skewSymmetricToVector(Rnew.transpose()*(Rnew*vectorToSkewSymmetric(x_bm_diff.segment(21, 3)) - Rmean*vectorToSkewSymmetric(x_bm_diff.segment(15, 3)) + (h_/2)*(ToMatrix(force_new.segment(3, 9)) - ToMatrix(vel_con_new.segment(3, 9)))));

      f_bm_diff.segment(18, 6) = x_bm_diff.segment(12, 6) - x_bm_diff.segment(24, 6);
    }

    for (size_t j = 0; j < 24; j++) {
      df.Row(rb_.Index()*dim_per_body_rbs() + j).segment(dim_per_body_rbs()*rbs_.NumBodies(), 2) = f_bm_diff(j); 
    }
  }
}

#endif