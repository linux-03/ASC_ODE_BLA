#ifndef EQUATION_RBS_CC
#define EQUATION_RBS_CC

#include "equation_rbs.h"
#include <chrono>

RigidBodySystemEquation::RigidBodySystemEquation(RigidBodySystem& rbs, double step_size): rbs_(std::ref(rbs)), h_(step_size)
{
  this->dimF_ = rbs.NumBodies()*dim_per_body_rbs() + rbs.NumBeams()*dim_per_beam();
  this->dimX_ = rbs.NumBodies()*dim_per_body_rbs() + rbs.NumBeams()*dim_per_beam();

  Vector<double> x = rbs.ExpandState();
  rbs.ManageConstraints(x);

  for (size_t i = 0; i < rbs.NumBodies(); i++)
  {
    std::shared_ptr<RigidBodyEquation> eq = std::make_shared<RigidBodyEquation>(rbs.Bodies(i), rbs, step_size);

    this->func_.addFunction(eq);
  }
}

size_t RigidBodySystemEquation::DimX() const {
  return this->dimX_;
}
size_t RigidBodySystemEquation::DimF() const
{
  return this->dimF_;
}
size_t RigidBodySystemEquation::BodyDimensions() const
{
  return this->rbs_.get().NumBodies()*dim_per_body_rbs();
}

void RigidBodySystemEquation::Evaluate (VectorView<double> x, VectorView<double> f) const
{ 
  Vector<double> x_large(rbs_.get().NumBodies()*dim_per_body() + rbs_.get().NumBeams()*dim_per_beam());

  ExpandToFullRotation(x, x_large, rbs_.get().NumBodies());

  func_.Evaluate(x_large, f.segment(0, this->BodyDimensions()));

  for (size_t i = 0; i < rbs_.get().NumBeams(); i++ )
  {
    f(this->BodyDimensions() + i*2) = this->rbs_.get().Constraint(x_large, i);
    f(this->BodyDimensions() + i*2 + 1) = this->rbs_.get().Velocity_Constraint(x_large, i);
  }
}

void RigidBodySystemEquation::EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
{

  Vector<double> x_large(rbs_.get().NumBodies()*dim_per_body() + rbs_.get().NumBeams()*dim_per_beam());

  ExpandToFullRotation(x, x_large, rbs_.get().NumBodies());

  df.setConstant(0);
  Matrix<double> df_block = df.Block(0, 0,  rbs_.get().NumBodies()*dim_per_body_rbs(), dimX_);
  func_.EvaluateDeriv(x_large, df_block);
  df.Block(0, 0,  rbs_.get().NumBodies()*dim_per_body_rbs(), dimX_) = df_block;

  Vector<AutoDiff<24, double>> x_diff = x_large;

  size_t prev_index = 0;

  Vector<AutoDiff<24, double>> f_diff(2*rbs_.get().NumBeams());

  for (size_t bd_index = 0; bd_index < rbs_.get().NumBodies(); bd_index++) {
    
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

    SetRotGradien(x.segment(dim_per_body_rbs()*bd_index + 3, 3), x_diff.segment(dim_per_body()*bd_index + 3, 9));

    for (size_t i = 0; i < rbs_.get().NumBeams(); i++ )
    {
      f_diff(i*2) = this->rbs_.get().Constraint(x_diff, i);
      f_diff(i*2 + 1) = this->rbs_.get().Velocity_Constraint(x_diff, i); 
    }

    for (size_t j = 0; j < rbs_.get().NumBeams(); j++) {
      df.Row(this->BodyDimensions() + 2*j).segment(dim_per_body_rbs()*bd_index, dim_per_body_rbs()) = f_diff(2*j); 
      df.Row(this->BodyDimensions() + 2*j + 1).segment(dim_per_body_rbs()*bd_index, dim_per_body_rbs()) = f_diff(2*j + 1); 
    }
  }
  //dNumeric(*(this), x, df);
}

void simulate_rbs(RigidBodySystem& rbs, double step_size, size_t steps_)
{
  RigidBodySystemEquation eqrb(rbs, step_size);

  Vector<double> x = rbs.ExpandState();

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < steps_; i++) {
    NewtonSolver(eqrb, x);

    rbs.SaveState(x);

    std::cout << x << std::endl;
    //std::cout << i << ": " << x(1)*9.81 + x.segment(18, 3).squaredNorm() << std::endl << std::endl; 
    //Energy : std::cout << x(1)*9.81 + 0.5*x.segment(6, 3).squaredNorm() << std::endl;
    // Pot Energy : std::cout << x(1)*9.81 << std::endl;
    // std::cout << 0.5*x.segment(18, 3).squaredNorm() << std::endl;
    //std::cout << x(2) << std::endl;
  }

  // Record the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
}

#endif