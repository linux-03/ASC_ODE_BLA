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
    Vector<double> p_half_trans_; 
    Vector<double> p_half_skew_;   // Rotation matrix
    Vector<double> v_trans_;
    Vector<double> v_skew_;
    Vector<double> p_trans_;
    Vector<double> p_skew_;    // Momentum vector (6D: 3D translation + 3D rotation)
    double lambda_;
    double mu_;
    size_t constraint_number_half_ = 0;
    std::vector<size_t> beams_;
    std::unique_ptr<Matrix<double>> constraints_ = nullptr;
    Vector<double> force_;
    Matrix<double> mass_matrix_;
    Matrix<double> mass_matrix_inv_;
    Matrix<double> inertia_;
    double mass_;
    

  public:
    RigidBody();
    RigidBody(VectorView<double> q, VectorView<double> p);
    RigidBody(double mass, MatrixView<double> inertia);    

    RigidBody(const RigidBody& other);

    VectorView<double> q_trans();
    MatrixView<double> q();
    VectorView<double> p_half_trans();
    VectorView<double> p_half_skew();
    VectorView<double> v_trans();
    VectorView<double> v_skew();
    VectorView<double> p_trans();
    VectorView<double> p_skew();
    MatrixView<double> MassMatrix();
    MatrixView<double> MassMatrixInv();
    MatrixView<double> Inertia();
    double Mass();
    double& lambda();
    double& mu();
    Vector<double> Vector_q();
    void recalcMassMatrix();

    size_t& Index();
    std::vector<size_t>& Beams();
    void addBeam(size_t i, size_t constraint_number);
    Matrix<double>& Constraints();
    Vector<double>& Force();
    size_t ConstraintNumberHalf() const;

    void add_vertex(const std::array<double, 3>& v);
    void add_normal(const std::array<double, 3>& n);
    std::vector<std::array<double, 3>> vertices();
    std::vector<std::array<double, 3>> normals();

    std::vector<std::array<double, 3>> vertices_;
    std::vector<std::array<double, 3>> normals_;
};

#endif