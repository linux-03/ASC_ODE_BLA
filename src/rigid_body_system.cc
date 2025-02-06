#ifndef RIGID_BODY_SYSTEM_CC
#define RIGID_BODY_SYSTEM_CC

#include "rigid_body_system.h"
#include "Beam.h"
#include "Connector.h"
#include <assert.h>
#include <memory.h>

using namespace ASC_bla;

RigidBodySystem::RigidBodySystem(): num_beams_(0), num_springs_(0) {}

template<typename T>
T RigidBodySystem::Potential(const VectorView<T> x)
{
  T pot = 0;
  for (size_t i = 0; i < NumBodies(); i++)
  {
    pot += x.segment(i*dim_per_body(), 3).dot(gravity_);
  }

  return pot;
}

size_t& RigidBodySystem::NumBeams() 
{
  return this->num_beams_;
}

size_t& RigidBodySystem::NumSprings()
{
  return this->num_springs_;
}

size_t RigidBodySystem::NumBodies()
{
  return this->bodies_.size();
}

RigidBody& RigidBodySystem::Bodies(size_t i)
{
  assert( i < this->bodies_.size());
  return this->bodies_[i].get();
}

std::vector<RigidBody> RigidBodySystem::Bodies() {
  std::vector<RigidBody> vec;
  for (std::reference_wrapper<RigidBody> rb: bodies_)
  {
    vec.push_back(rb.get());
  }
  return vec;
}

Vector<double> RigidBodySystem::connectorPosition(Connector c)
{
  if (c.Fix())
  {
    return c.RefPosition();
  }

  Vector<double> pos(3);
  pos.setConstant(0);

  pos += Bodies(c.BodyIndex()).q_trans() + Bodies(c.BodyIndex()).q()*c.RefPosition();

  return pos;
}

Beam& RigidBodySystem::Beams(size_t i)
{
  assert( i < this->beams_.size());
  return this->beams_[i].get();
}

std::vector<Beam> RigidBodySystem::Beams() {
  std::vector<Beam> vec;
  for (std::reference_wrapper<Beam> bm: beams_)
  {
    vec.push_back(bm.get());
  }
  return vec;
}

Spring& RigidBodySystem::Springs(size_t i)
{
  assert( i < this->springs_.size());
  return this->springs_[i].get();
}

void RigidBodySystem::SaveState(const VectorView<double> x)
{
  for (size_t i = 0; i < this->NumBodies(); i++)
  {
    RigidBody& rb = this->Bodies(i);
    size_t base = i*dim_per_body();
    rb.q_trans() = x.segment(base, 3);
    //std::cout << ToMatrix(x.segment(3, 9));
    rb.q() = ToMatrix(x.segment(base + 3, 9));
    rb.p_trans() = x.segment(base + 18, 3);
    rb.p_skew() = x.segment(base + 21, 3);
  }
  for (size_t i = 0; i < this->NumBeams(); i++)
  {
    Beam& bm = this->Beams(i);
    size_t base = this->NumBodies()*dim_per_body();

    bm.Lambda() = x(base + 2*i);
    bm.Mu() = x(base + 2*i + 1);
  }
  this->ManageConstraints(x);
}

void RigidBodySystem::ManageConstraints(const VectorView<double> x)
{
  for (size_t i = 0; i < this->NumBodies(); i++)
  {
    this->Bodies(i).Force() = this->PotentialGradient<double>(x, i);
    this->Bodies(i).Constraints() = this->JacobianConstraint<double>(x, i);    
  }
}

Vector<double> RigidBodySystem::ExpandState()
{
  Vector<double> x((this->NumBodies())*dim_per_body() + 2*(this->NumBeams()));
  x.setConstant(0);
  
  for (size_t i = 0; i < this->NumBodies(); i++)
  {
    RigidBody& rb = this->Bodies(i);
    size_t base = i*dim_per_body();

    x.segment(base, 3) = rb.q_trans();
    //std::cout << ToMatrix(x.segment(3, 9));
    x.segment(base + 3, 9) = ToVector(rb.q());
    x.segment(base + 18, 3) = rb.p_trans();
    x.segment(base + 21, 3) = rb.p_skew();
  }

  return x;
}

void RigidBodySystem::add(RigidBody& rb)
{
  rb.Index() = this->NumBodies();
  this->bodies_.push_back(std::ref(rb));
}

void RigidBodySystem::add(Beam& bm)
{
  this->beams_.push_back(std::ref(bm));
  bm.Index() = this->num_beams_;
  this->num_beams_++;
  this->Bodies(bm.BodyIndexA()).addBeam(bm.Index());
  this->Bodies(bm.BodyIndexB()).addBeam(bm.Index());
  double res = Norm(bm.PositionA(this->Bodies(bm.BodyIndexA()).Vector_q()) - bm.PositionB(this->Bodies(bm.BodyIndexB()).Vector_q()));
  bm.Length() = res;
}

void RigidBodySystem::add(Spring& spr)
{
  this->springs_.push_back(std::ref(spr));
  spr.Index() = this->num_springs_;
  this->num_springs_++;
}
#endif