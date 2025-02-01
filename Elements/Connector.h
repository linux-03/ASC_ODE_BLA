#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "../helper/helper.h"
#include "vector.h"

using namespace ASC_bla;

class Connector {
    bool fix_;
    size_t body_index_;
    Vector<double> pos_ = {0,0,0};


   
  public:
    Connector() = default;
    Connector(Vector<double> pos);
    Connector(Vector<double> pos, size_t body_index_);

    size_t& BodyIndex();
    bool& Fix();
    
    template<typename T>
    Vector<T> Position(const VectorView<T> q) {
      if (this->fix_)
      {
        return this->pos_;
      }
      return q.segment(0, 3) + ToMatrix(q.segment(3, 9))*(this->pos_);
    }

    VectorView<double> RefPosition();

    double RefPosition(size_t i);
};

#endif