#ifndef SPRING_H
#define SPRING_H

#include "Connector.h"
#include <vector.h>

class Spring{
    size_t index_;
    double stiffness_;
    Connector connector_a_;
    Connector connector_b_;

  public:
    Spring (Connector connector_a, Connector connector_b);

    size_t& Index();
    double& Stiffness();
    Connector& ConnectorA();
    Connector& ConnectorB();
};

#endif