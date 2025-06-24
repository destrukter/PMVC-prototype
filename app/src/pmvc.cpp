#include "pmvc.hpp"

void computePMVC_GPU_lipman(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::MatrixXd &obj_V,
    Eigen::MatrixXd &pmvc_coords)
{
}

void checkVisibility(
    const Eigen::MatrixXd &cage_V,
    const Eigen::VectorXd &mesh_Vertex,
    const Eigen::MatrixXi &cage_F,
    Eigen::MatrixXd &visible_V)
{
    Eigen::Vector3d origin(mesh_Vertex[0], mesh_Vertex[1], mesh_Vertex[2]);

    std::vector<Eigen::Vector3d> visible_points;

    for (int i = 0; i < cage_V.rows(); i++)
    {
        Eigen::Vector3d dir = cage_V.row(i) - origin;
        dir.normalize();

        std::vector<igl::Hit<double>> hits;

        bool intersected = igl::ray_mesh_intersect(origin, dir, cage_V, cage_F, hits);

        if (!intersected)
        {
            visible_points.push_back(cage_V.row(i));
        }
        else
        {
            double dist_to_point = (cage_V.row(i) - origin).norm();
            double dist_to_hit = hits[0].t;
            if (dist_to_hit > dist_to_point + 1e-6)
            {
                visible_points.push_back(cage_V.row(i));
            }
        }
    }

    visible_V.resize(visible_points.size(), 3);
    for (int i = 0; i < (int)visible_points.size(); i++)
    {
        visible_V.row(i) = visible_points[i];
    }
}