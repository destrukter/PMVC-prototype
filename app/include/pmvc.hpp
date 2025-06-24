#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/read_triangle_mesh.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/Hit.h>

void computePMVC_GPU_lipman(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::MatrixXd &obj_V,
    Eigen::MatrixXd &pmvc_coords);

void computePMVC_GPU_ours(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::MatrixXd &obj_V,
    Eigen::MatrixXd &pmvc_coords);

void checkVisibility(
    const Eigen::MatrixXd &cage_V,
    const Eigen::VectorXd &mesh_Vertex,
    const Eigen::MatrixXi &cage_F,
    Eigen::MatrixXd &visible_V);