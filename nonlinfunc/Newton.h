#ifndef Newton_h
#define Newton_h

#include <iomanip>
#include <matrix.h>
#include "nonlinfunc.h"
#include <vector.h>
#include <functional>

namespace ASC_ode {

  void NewtonSolver (NonlinearFunction& func, VectorView<double> x,
                    double tol = 1e-10, int maxsteps = 10,
                    std::function<void(int,double,VectorView<double>)> callback = nullptr);
                    
}

#endif
