#ifndef RIGID_BODY_CC
#define RIGID_BODY_CC

#include "rigid_body.h"

using namespace ASC_bla;

RigidBody::RigidBody(): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3)
{
  q_trans_.setConstant(0);
  q_.setConstant(0);
  p_trans_.setConstant(0);
  p_skew_.setConstant(0);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);

  q_.diagonal(1);
}

RigidBody::RigidBody(VectorView<double> q, VectorView<double> p): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3)
{
  q_trans_ = q.segment(0, 3);
  q_ = ToMatrix(q.segment(3, 9));
  p_trans_ = p.segment(0, 3);
  p_skew_ = p.segment(3, 3);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);
}

VectorView<double> RigidBody::q_trans()
{
  return q_trans_;
}

MatrixView<double> RigidBody::q()
{
  return q_;
}

VectorView<double> RigidBody::p_half_trans()
{
  return p_half_trans_;
}

VectorView<double> RigidBody::p_half_skew()
{
  return p_half_skew_;
}

VectorView<double> RigidBody::v_trans()
{
  return v_trans_;
}

VectorView<double> RigidBody::v_skew()
{
  return v_skew_;
}

VectorView<double> RigidBody::p_trans()
{
  return p_trans_;
}

VectorView<double> RigidBody::p_skew()
{
  return p_skew_;
}

double& RigidBody::lambda()
{
  return lambda_;
}

double& RigidBody::mu()
{
  return mu_;
}
Vector<double> RigidBody::Vector_q()
{
  Vector q_vec(12);
  q_vec.segment(0, 3) = q_trans_;
  q_vec.segment(3, 9) = ToVector(q_);
  return q_vec;
}

size_t& RigidBody::Index()
{
  return this->body_index_;
}

std::vector<size_t>& RigidBody::Beams()
{
  return this->beams_;
}

void RigidBody::addBeam(size_t i)
{
  if (std::find(this->beams_.begin(), this->beams_.end(), i) == this->beams_.end())
  {
    this->beams_.push_back(i);
    this->constraints_.reset(new Matrix(12, beams_.size()));
    this->constraints_->setConstant(0);
  }
}

MatrixView<double> RigidBody::Constraints()
{
  return *(this->constraints_);
}

VectorView<double> RigidBody::Force()
{
  return this->force_;
}

#endif