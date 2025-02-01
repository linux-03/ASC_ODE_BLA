#ifndef SPRING_CC
#define SPRING_CC

#include "Spring.h"

Spring::Spring(Connector connector_a, Connector connector_b)
{
  this->connector_a_ = connector_a;
  this->connector_b_ = connector_b;
  this->stiffness_ = 1;
  this->index_ = -1;
}

size_t& Spring::Index()
{
  return this->index_;
}
double& Spring::Stiffness()
{
  return this->stiffness_;
}
Connector& Spring::ConnectorA()
{
  return this->connector_a_;
}
Connector& Spring::ConnectorB()
{
  return connector_b_;
}

#endif