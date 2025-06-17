#ifndef CALCULATION_HPP
#define CALCULATION_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/per_face_normals.h>
#include <fstream>
#include <string>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <vector>
#include <array>

struct Face
{
    std::array<int, 3> vi;
};

struct Mesh
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
};

Eigen::MatrixXd applyDeformation(
    const Eigen::MatrixXd &weights,
    const Eigen::MatrixXd &VdeformedCage);

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

bool is_face_visible(const Eigen::Vector3d &v, const Face &f, const Mesh &cage);

std::vector<Face> clip_face_along_visibility(const Eigen::Vector3d &v, const Face &f, const Mesh &cage);

Mesh construct_temporary_mesh(const std::vector<Face> &faces);

Mesh build_visibility_frustum(const Eigen::Vector3d &v, const Eigen::Vector3d &face_center);

Mesh clip_face_along_visibility(
    const Eigen::Vector3d &v_pos,
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    int face_idx);

void compute_pmvc(
    const Eigen::MatrixXd &obj_V,
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    Eigen::MatrixXd &mvc_coords);

#endif // CALCULATION_HPP
