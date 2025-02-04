#ifndef HELPER_cc
#define HELPER_cc

#include "helper.h"

size_t dim_per_body_rbs()
{
  return 24;
}

size_t dim_per_body()
{
  return 33;
}

void ExpandToFullRotation(VectorView<double> x, VectorView<double> x_exp, size_t num_bodies)
{
  for (size_t i = 0; i < num_bodies; i++)
  {
    x_exp.segment(i*dim_per_body(), 3) = x.segment(i*dim_per_body_rbs(), 3);
    x_exp.segment(i*dim_per_body() + 3, 9) = ToVector(AxisToMatrix(x.segment(i*dim_per_body_rbs() + 3, 3)));
    x_exp.segment(i*dim_per_body() + 12, 18) = x.segment(i*dim_per_body_rbs() + 6, 18);
    x_exp.segment(i*dim_per_body() + 30, 3) = x.segment(i*dim_per_body_rbs() + 3, 3);
  }

  size_t end = (x.Size() - num_bodies*dim_per_body_rbs())/2;

  if (end != 0)
  {
    x_exp.segment(num_bodies*dim_per_body(), end*2) = x.segment(num_bodies*dim_per_body_rbs(), end*2);
  }
}

size_t dim_per_motion()
{
  return 12;
}

size_t dim_per_beam()
{
  return 2;
}

Matrix<double> AxisToMatrix(Vector<double> ax) 
{
#ifdef DEBUG_MODE
  assert(ax.Size() == 3);
#endif

  double angle = ax.norm();

  Matrix<double> res = Diagonal(3, 1);

  if (angle == 0) {
    return res;
  }

  Matrix<double> H = vectorToSkewSymmetric(ax);
  
  res = res + (sin(angle)/angle) * H + (1 - cos(angle))/(angle*angle) * H*H;

  return res;
}

Vector<double> MatrixToAxis(const MatrixView<double> R) {
    Vector<double> rot(3);
    rot.setConstant(0);
    double trace = R(0,0) + R(1,1) + R(2,2);
    double angle = std::acos((trace - 1.0) / 2.0);

    if(angle == 0) return rot;
    
    double s = 2 * std::sin(angle);
    rot(0) = (R(2, 1)- R(1, 2)) / s;
    rot(1) = (R(0, 2) - R(2, 0)) / s;
    rot(2) = (R(1, 0) - R(0, 1)) / s;
    
    return angle*rot;
}


void SetRotGradien(VectorView<double> x, VectorView<AutoDiff<24, double>> x_diff)
{
  double norm = x.norm();
  if (norm == 0) {
    //Gradinet for var a
    x_diff(5).DValue(3) = -1;
    x_diff(7).DValue(3) = 1;
    // Gradient for var b
    x_diff(2).DValue(4) = 1;
    x_diff(6).DValue(4) = -1;
    // Gradient for var c
    x_diff(1).DValue(5) = -1;
    x_diff(3).DValue(5) = 1;
  } else {
    Vector<AutoDiff<24, double>> axis_diff = x;

    AutoDiff<24, double> angle = x.norm();

    Matrix<AutoDiff<24, double>> H = vectorToSkewSymmetric(axis_diff);

    Matrix<AutoDiff<24, double>> res = Diagonal(3, 1) +  (sin(angle)/angle) * H + (1 - cos(angle))/(angle*angle) * H*H;

    x_diff = ToVector(res);
  }
}

#endif