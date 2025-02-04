#ifndef HELPER_H
#define HELPER_H


#include <iostream>
#include <matrix.h>
#include <vector.h>
#include <assert.h>

using namespace ASC_bla;

template<typename T>
Matrix<T> ToMatrix(VectorView<T> vec) {
  //std::cout << vec << std::endl;
  Matrix<T> res(3, 3);
  res.Row(0) = vec.segment(0, 3);  // First row
  res.Row(1) = vec.segment(3, 3);  // Second row
  res.Row(2) = vec.segment(6, 3);  // Third row

  return res;
}

template<typename T>
Vector<T> ToVector(MatrixView<T> m) {
  Vector<T> res(9);
  res.setConstant(0);

  res.segment(0, 3) = m.Row(0);
  res.segment(3, 3) = m.Row(1);
  res.segment(6, 3) = m.Row(2);

  return res;
}

// Function to convert a 3D vector to a skew-symmetric matrix (hat operator)
template<typename T>
Matrix<T> vectorToSkewSymmetric(const VectorView<T> v) {
    Matrix<T> skewSymmetric(3, 3);
    skewSymmetric(0, 0) = 0;
    skewSymmetric(0, 1) = -v(2);
    skewSymmetric(0, 2) = v(1);
    skewSymmetric(1, 0) = v(2);
    skewSymmetric(1, 1) = 0;
    skewSymmetric(1, 2) = -v(0);
    skewSymmetric(2, 0) = -v(1);
    skewSymmetric(2, 1) = v(0);
    skewSymmetric(2, 2) = 0;
    
    return skewSymmetric;
}

// Function to convert a skew-symmetric matrix to a 3D vector (inverse hat, vee operator)
template<typename T>
auto skewSymmetricToVector(const MatrixExpr<T>& skewMatrix) -> Vector<decltype(skewMatrix(0, 0)*1.0)> {
    // Check if the matrix is skew-symmetric

    Vector<decltype(1.0*skewMatrix(0, 0))> res = {skewMatrix(2, 1), skewMatrix(0, 2), skewMatrix(1, 0)};
    return res;
}

Matrix<double> AxisToMatrix(Vector<double> ax);

void ExpandToFullRotation(VectorView<double> x, VectorView<double> x_exp, size_t num_bodies);

void SetRotGradien(VectorView<double> x, VectorView<AutoDiff<24, double>> x_diff);

Vector<double> MatrixToAxis(const MatrixView<double> R);

size_t dim_per_body();

size_t dim_per_body_rbs();

size_t dim_per_motion();

size_t dim_per_beam();

/*
size_t start_q_trans();
size_t size_q_trans();
size_t start_q_rot();
size_t size_q_rot();
size_t start_v_trans();
size_t size_v_trans();
size_t start_v_skew();
size_t size_v_skew();
size_t start_phalf_trans();
size_t size_phalf_trans();
size_t start_phalf_skew();
size_t size_phalf_skew();
*/


#endif