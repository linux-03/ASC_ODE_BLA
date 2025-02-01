#ifndef CONNECTOR_CC
#define CONNECTOR_CC

#include "Connector.h"

Connector::Connector(Vector<double> pos, size_t body_index_a)
{
  this->fix_ = false;
  this->body_index_ = body_index_a;
  this->pos_ = pos;
}
Connector::Connector(Vector<double> pos)
{
  this->fix_ = true;
  this->body_index_ = 0;
  this->pos_ = pos;
}
size_t& Connector::BodyIndex()
{
  return this->body_index_;
}
bool& Connector::Fix()
{
  return this->fix_;
}
VectorView<double> Connector::RefPosition()
{
  return pos_;
}
double Connector::RefPosition(size_t i)
{
  return pos_(i);
}

#endif