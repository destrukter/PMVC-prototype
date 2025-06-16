#include "calculation.hpp"

#include <igl/read_triangle_mesh.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/edges.h>
#include <igl/writeOBJ.h>

int main(int argc, char *argv[])
{
    Eigen::MatrixXd Vmesh;
    Eigen::MatrixXi Fmesh;

    Eigen::MatrixXd Vcage;
    Eigen::MatrixXi Fcage;
    Eigen::MatrixXi Ecage;

    Eigen::MatrixXd VcageDeformed;
    Eigen::MatrixXi FcageDeformed;
    Eigen::MatrixXi EcageDeformed;

    std::string meshPath = "";
    std::string cagePath = "";
    std::string deformedCagePath = "";

    if (argc == 4)
    {
        meshPath = argv[1];
        cagePath = argv[2];
        deformedCagePath = argv[3];
    }

    const std::string defaultMesh = "/home/destukter/Thesis/PMVC-test/meshes/armadilloman.obj";
    const std::string defaultCage = "/home/destukter/Thesis/PMVC-test/meshes/armadilloman_cages_triangulated.obj";
    const std::string defaultDeformedCage = "/home/destukter/Thesis/PMVC-test/meshes/armadilloman_cages_triangulated_deformed_1.obj";

    std::string finalMeshPath = fileExists(meshPath) ? meshPath : defaultMesh;
    std::string finalCagePath = fileExists(cagePath) ? cagePath : defaultCage;
    std::string finalDeformedCagePath = fileExists(deformedCagePath) ? deformedCagePath : defaultDeformedCage;

    igl::read_triangle_mesh(finalMeshPath, Vmesh, Fmesh);
    igl::read_triangle_mesh(finalCagePath, Vcage, Fcage);
    igl::read_triangle_mesh(finalDeformedCagePath, VcageDeformed, FcageDeformed);

    igl::edges(Fcage, Ecage);
    igl::edges(FcageDeformed, EcageDeformed);

    // TODO add calculation for PMVC
    Eigen::MatrixXd weights;
    computeMVC(Vcage, Fcage, Vmesh, weights);
    // Eigen::MatrixXd weights = computePMVCForMesh(Vmesh, Vcage, Fcage);

    Eigen::MatrixXd VmeshDeformed = applyDeformation(weights, VcageDeformed);

    double offset = (Vmesh.col(0).maxCoeff() - Vmesh.col(0).minCoeff()) * 1.5;
    VmeshDeformed.col(0).array() += offset;
    VcageDeformed.col(0).array() += offset;

    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(Vmesh, Fmesh);
    viewer.data().set_face_based(false);

    viewer.append_mesh();
    // viewer.data(1).set_points(Vcage, Eigen::RowVector3d(1, 0, 0));
    viewer.data(1).set_edges(Vcage, Ecage, Eigen::RowVector3d(0, 0, 1));

    viewer.append_mesh();
    // viewer.data(2).set_points(VcageDeformed, Eigen::RowVector3d(1, 0, 0));
    viewer.data(2).set_edges(VcageDeformed, EcageDeformed, Eigen::RowVector3d(0, 0, 1));

    viewer.append_mesh();
    viewer.data(3).set_mesh(VmeshDeformed, Fmesh);
    viewer.data(3).set_face_based(false);
    viewer.data(3).set_colors(Eigen::RowVector3d(1, 0, 0));

    std::cout << "Starting viewer..." << std::endl;
    viewer.launch();

    return 0;
}
