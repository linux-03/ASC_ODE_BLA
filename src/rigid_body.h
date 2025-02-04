#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include <iostream>
#include <vector.h>
#include <matrix.h>
#include <helper.h>
#include <cmath>
#include <algorithm>

using namespace std;


class RigidBody { 
    size_t body_index_ = -1;
    Vector<double> q_trans_;
    Matrix<double> q_;
    Vector<double> axis_;
    Vector<double> p_half_trans_; 
    Vector<double> p_half_skew_;   // Rotation matrix
    Vector<double> v_trans_;
    Vector<double> v_skew_;
    Vector<double> p_trans_;
    Vector<double> p_skew_;    // Momentum vector (6D: 3D translation + 3D rotation)
    double lambda_;
    double mu_;
    std::vector<size_t> beams_;
    std::unique_ptr<Matrix<double>> constraints_ = nullptr;
    Vector<double> force_;

  public:
    RigidBody();
    RigidBody(VectorView<double> q, VectorView<double> p);
    

    VectorView<double> q_trans();
    MatrixView<double> q();
    VectorView<double> axis();
    VectorView<double> p_half_trans();
    VectorView<double> p_half_skew();
    VectorView<double> v_trans();
    VectorView<double> v_skew();
    VectorView<double> p_trans();
    VectorView<double> p_skew();
    double& lambda();
    double& mu();
    Vector<double> Vector_q();

    size_t& Index();
    std::vector<size_t>& Beams();
    void addBeam(size_t);
    MatrixView<double> Constraints();
    VectorView<double> Force();
};

#endif