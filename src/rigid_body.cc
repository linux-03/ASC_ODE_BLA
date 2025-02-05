#ifndef RIGID_BODY_CC
#define RIGID_BODY_CC

#include "rigid_body.h"

using namespace ASC_bla;

RigidBody::RigidBody(): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), axis_(3)
{
  q_trans_.setConstant(0);
  q_.setConstant(0);
  axis_.setConstant(0);
  p_trans_.setConstant(0);
  p_skew_.setConstant(0);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);

  q_.diagonal(1);
}

RigidBody::RigidBody(VectorView<double> q, VectorView<double> p): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), axis_(3)
{
  q_trans_ = q.segment(0, 3);
  q_ = ToMatrix(q.segment(3, 9));
  axis_ = MatrixToAxis(q_);
  p_trans_ = p.segment(0, 3);
  p_skew_ = p.segment(3, 3);
  p_half_trans_.setConstant(0);
  p_half_skew_.setConstant(0);
  v_trans_.setConstant(0);
  v_skew_.setConstant(0);
}

RigidBody::RigidBody(const RigidBody& other): force_(12), q_trans_(3), q_(3, 3), p_trans_(3),
p_skew_(3), p_half_trans_(3), p_half_skew_(3), v_trans_(3), v_skew_(3), axis_(3) {
  q_trans_ = other.q_trans_;
  q_ = other.q_;
  axis_ = other.axis_;
  p_trans_ = other.p_trans_;
  p_skew_ = other.p_skew_;
  p_half_trans_ = other.p_half_trans_;
  p_half_skew_ = other.p_half_skew_;
  v_trans_ = other.v_trans_;
  v_skew_ = other.v_skew_;
  lambda_ = other.lambda_;
  mu_ = other.mu_;
  beams_ = other.beams_;
  constraints_.reset(new Matrix(12, other.beams_.size()));
  *(constraints_) = *(other.constraints_);
  force_ = other.force_;
  vertices_ = other.vertices_;
  normals_ = other.normals_;
}

VectorView<double> RigidBody::q_trans()
{
  return q_trans_;
}

MatrixView<double> RigidBody::q()
{
  return q_;
}

VectorView<double> RigidBody::axis()
{
  return axis_;
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

// Add a vertex
    void RigidBody::add_vertex(const std::array<double, 3>& v) {
        vertices_.push_back(v);
    }

    // Add a normal
    void RigidBody::add_normal(const std::array<double, 3>& n) {
        normals_.push_back(n);
    }

    // Get vertices
    std::vector<std::array<double, 3>> RigidBody::vertices() {
        return vertices_;
    }

    // Get normals
    std::vector<std::array<double, 3>> RigidBody::normals() {
        return normals_;
    }

#endif