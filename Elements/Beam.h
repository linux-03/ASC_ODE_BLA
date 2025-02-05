#ifndef BEAM_H
#define BEAM_H

#include <vector.h>
#include "Connector.h"

using namespace ASC_bla;

class Beam {
    size_t index_;
    double length_;
    double lambda_;
    double mu_;
    Connector connector_a_;
    Connector connector_b_;

  public:
    Beam (Connector connector_a, Connector connector_b);
    Beam (const Beam& other);

    size_t& Index();
    double& Length();
    size_t BodyIndexA();
    size_t BodyIndexB();
    double& Lambda();
    double& Mu();
    Connector& ConnectorA();
    Connector& ConnectorB();
    template<typename T>
    Vector<T> PositionA(const VectorView<T> q) {
      return this->connector_a_.Position(q);
    }
    template<typename T>
    Vector<T> PositionB(const VectorView<T> q) {
      return this->connector_b_.Position(q);
    }
};

#endif