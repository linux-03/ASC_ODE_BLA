#ifndef EQUATION_RBS_H
#define EQUATION_RBS_H

#include <memory.h>
#include "rigid_body_system.h"
#include "rigid_body.h"
#include "equation_rb.h"
#include "nonlinfunc.h"
#include "helper.h"
#include "Newton.h"
#include "autodiff.h"

using namespace ASC_ode;

class RigidBodySystemEquation: public NonlinearFunction
{
    size_t dimF_;
    size_t dimX_;
    std::reference_wrapper<RigidBodySystem> rbs_;
    StackedFunction_large_input func_;
    Vector<double> x_;
    double h_;
    

  public:
    RigidBodySystemEquation(RigidBodySystem& rbs, double step_size);

    size_t DimX() const;
    size_t DimF() const;
    void Evaluate (VectorView<double> x, VectorView<double> f) const;
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const;

    void step();
    Vector<double>& x();

    size_t BodyDimensions() const;
};

RigidBodySystemEquation assemble(RigidBodySystem& rbs, double step_size);
void simulate_rbs(RigidBodySystem& rbs, double step_size, size_t steps_, std::function<void(int,double,VectorView<double>)> callback = nullptr);

#endif