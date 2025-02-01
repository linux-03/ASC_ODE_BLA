#ifndef Newton_cc
#define Newton_cc

#include "Newton.h"
#include "iostream"

namespace ASC_ode {

  void NewtonSolver (NonlinearFunction& func, VectorView<double> x, double tol, int maxsteps, std::function<void(int,double,VectorView<double>)> callback)
  {
    Vector<double> res(func.DimF());
    res.setConstant(0);
    Matrix<double> fprime(func.DimF(), func.DimX());
    fprime.setConstant(0);

    size_t count = 0;
    //std::cout << "x = " << x << std::endl;
    for (int i = 0; i < maxsteps; i++)
      {
        
        func.Evaluate(x, res);
        //std::cout << "f: " << res << std::endl << std::endl;
        //std::cout<< "eval" << std::endl;
        
        func.EvaluateDeriv(x, fprime);

        //std::cout << x << std::endl << std::endl;
        //std::cout << fprime << std::endl;
        Matrix<double> fprime_inv = inverse(fprime);

        
        
        x -= fprime_inv*res;

        //std::cout << fprime_inv << std::endl << std::endl << std::endl;
        //std::cout << x << std::endl;
        //throw std::domain_error("Newton did not converge");
        //std::cout << x << std::endl << std::endl;
        
        double err = res.norm();
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }
}

#endif
