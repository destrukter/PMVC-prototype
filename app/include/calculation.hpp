#ifndef CALCULATION_HPP
#define CALCULATION_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/per_face_normals.h>
#include <fstream>
#include <string>

/*
/// Compute Point Mean Value Coordinates (PMVC) for a single query point.
Eigen::VectorXd computePMVC(
    const Eigen::MatrixXd &V,   // cage vertices (Nx3)
    const Eigen::MatrixXi &F,   // cage faces (Mx3)
    const Eigen::RowVector3d &p // query point (1x3)
);
*/

/*
/// Compute PMVC for an entire mesh of query points.
Eigen::MatrixXd computePMVCForMesh(
    const Eigen::MatrixXd &Vmesh, // query mesh vertices (Px3)
    const Eigen::MatrixXd &Vcage, // cage vertices (Nx3)
    const Eigen::MatrixXi &Fcage  // cage faces (Mx3)
);
*/

/// Apply deformation using PMVC weights and deformed cage.
Eigen::MatrixXd applyDeformation(
    const Eigen::MatrixXd &weights,      // PMVC weights (PxN)
    const Eigen::MatrixXd &VdeformedCage // deformed cage vertices (Nx3)
);

/// Utility to check if a file exists.
bool fileExists(const std::string &path);

/// Compute Mean Value Coordinates (MVC) for a single vertex in a robust way.
void computeMVCForOneVertexSimple(
    const Eigen::MatrixXd &C,  // cage vertices (Nx3)
    const Eigen::MatrixXi &CF, // cage faces (Mx3)
    const Eigen::Vector3d eta, // query point (3,)
    Eigen::VectorXd &weights,  // output weights (Nx1)
    Eigen::VectorXd &w_weights // intermediate weights (Nx1)
);

/// Compute MVC for multiple query points.
void computeMVC(
    const Eigen::MatrixXd &C,     // cage vertices (Nx3)
    const Eigen::MatrixXi &CF,    // cage faces (Mx3)
    const Eigen::MatrixXd &eta_m, // query points (Px3)
    Eigen::MatrixXd &phi          // output weights (NxP)
);

#endif // CALCULATION_HPP
