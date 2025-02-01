#ifndef CONNECTOR_CC
#define CONNECTOR_CC

#include "Beam.h"

Beam::Beam(Connector connector_a, Connector connector_b)
{
  this->connector_a_ = connector_a;
  this->connector_b_ = connector_b;
  this->length_ = 0;
  this->index_ = -1;
}

size_t& Beam::Index()
{
  return this->index_;
}
double& Beam::Length()
{
  return this->length_;
}
size_t Beam::BodyIndexA()
{
  return this->connector_a_.BodyIndex();
}
size_t Beam::BodyIndexB()
{
  return this->connector_b_.BodyIndex();
}
Connector& Beam::ConnectorA()
{
  return this->connector_a_;
}
Connector& Beam::ConnectorB()
{
  return connector_b_;
}
double& Beam::Lambda()
{
  return this->lambda_;
}
double& Beam::Mu()
{
  return this->mu_;
}
#endif