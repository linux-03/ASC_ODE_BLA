#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "../helper/helper.h"
#include "vector.h"
#include "rigid_body.h"

using namespace ASC_bla;

enum struct ConnectorType {
    FIX,
    FREE,
    SPHERICAL
};

class Connector {
    ConnectorType type_;
    size_t body_index_;
    Vector<double> pos_ = {0,0,0};
    Vector<double> initial_pos_ = {0,0,0};
    Vector<double> initial_body_trans_ = {0,0,0};
    Matrix<double> initial_body_rot_ = Matrix<double>(3,3);


   
  public:
    Connector() = default;
    Connector(Vector<double> pos);
    Connector(Vector<double> pos, RigidBody& rb, ConnectorType type = ConnectorType::FREE);

    Connector(const Connector& other);

    size_t& BodyIndex();
    ConnectorType& Type();
    
    template<typename T>
    Vector<T> BodyFramePosition(const VectorView<T> q) {
      return q.segment(0, 3) + ToMatrix(q.segment(3, 9))*(this->pos_);
    }

    template<typename T>
    Vector<T> AbsPosition(const VectorView<T> q) {
      if (this->type_ == ConnectorType::FIX) {
        return this->pos_;
      }
      return q.segment(0, 3) + ToMatrix(q.segment(3, 9))*(this->pos_);
    }

    Vector<double> RefPosition();

    double RefPosition(size_t i);
    Vector<double> InitialPosition();
    Vector<double> InitialBodyTranslation();
    Matrix<double> InitialBodyRotation();
};

#endif