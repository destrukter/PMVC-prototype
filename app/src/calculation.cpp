#include "calculation.hpp"

Eigen::MatrixXd
applyDeformation(
    const Eigen::MatrixXd &weights,
    const Eigen::MatrixXd &VdeformedCage)
{
    return weights.transpose() * VdeformedCage;
}

bool fileExists(const std::string &path)
{
    std::ifstream file(path);
    return file.good();
}

/// A simple but more robust scheme to compute MVC supposed in https://github.com/superboubek/QMVC
void computeMVCForOneVertexSimple(Eigen::MatrixXd const &C, Eigen::MatrixXi const &CF,
                                  Eigen::Vector3d eta, Eigen::VectorXd &weights, Eigen::VectorXd &w_weights)
{
    double epsilon = 0.000000001;

    auto const num_vertices_cage = C.rows();
    auto const num_faces_cage = CF.rows();

    w_weights.setZero();
    weights.setZero();
    double sumWeights = 0.0;

    Eigen::VectorXd d(num_vertices_cage);
    d.setZero();
    Eigen::MatrixXd u(num_vertices_cage, 3);

    for (unsigned int v = 0; v < num_vertices_cage; ++v)
    {
        const Eigen::Vector3d cage_vertex = C.row(v);
        d(v) = (eta - cage_vertex).norm();
        if (d(v) < epsilon)
        {
            weights(v) = 1.0;
            return;
        }
        u.row(v) = (cage_vertex - eta) / d(v);
    }

    unsigned int vid[3];
    double l[3], theta[3], w[3];

    for (unsigned int t = 0; t < num_faces_cage; ++t)
    {
        // the Norm is CCW :
        for (unsigned int i = 0; i <= 2; ++i)
        {
            vid[i] = CF(t, i);
        }

        for (unsigned int i = 0; i <= 2; ++i)
        {
            const Eigen::Vector3d v_0 = u.row(vid[(i + 1) % 3]);
            const Eigen::Vector3d v_1 = u.row(vid[(i + 2) % 3]);
            l[i] = (v_0 - v_1).norm();
        }

        for (unsigned int i = 0; i <= 2; ++i)
        {
            const Eigen::Vector3d v_0 = u.row(vid[(i + 1) % 3]);
            const Eigen::Vector3d v_1 = u.row(vid[(i + 2) % 3]);
            theta[i] = 2. * asin((v_0 - v_1).norm() * .5);
        }

        // test in original MVC paper: (they test if one angle psi is close to 0: it is "distance sensitive" in the sense that it does not
        // relate directly to the distance to the support plane of the triangle, and the more far away you go from the triangle, the worse it is)
        // In our experiments, it is actually not the good way to do it, as it increases significantly the errors we get in the computation of weights and derivatives,
        // especially when evaluating Hfx, Hfy, Hfz which can be of norm of the order of 10^3 instead of 0 (when specifying identity on the cage, see paper)

        // simple test we suggest:
        // the determinant of the basis is 2*area(T)*d( eta , support(T) ), we can directly test for the distance to support plane of the triangle to be minimum

        const Eigen::Vector3d c_0 = C.row(vid[0]);
        const Eigen::Vector3d c_1 = C.row(vid[1]);
        const Eigen::Vector3d c_2 = C.row(vid[2]);
        double determinant = (c_0 - eta).dot((c_1 - c_0).cross(c_2 - c_0));
        double sqrdist = determinant * determinant / (4 * ((c_1 - c_0).cross(c_2 - c_0)).squaredNorm());
        double dist = std::sqrt(sqrdist);

        if (dist < epsilon)
        {
            // then the point eta lies on the support plane of the triangle
            double h = (theta[0] + theta[1] + theta[2]) * .5;
            if (M_PI - h < epsilon)
            {
                // eta lies inside the triangle t , use 2d barycentric coordinates :
                for (unsigned int i = 0; i <= 2; ++i)
                {
                    w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
                }
                sumWeights = w[0] + w[1] + w[2];

                w_weights.setZero();
                weights(vid[0]) = w[0] / sumWeights;
                weights(vid[1]) = w[1] / sumWeights;
                weights(vid[2]) = w[2] / sumWeights;
                return;
            }
        }

        Eigen::Vector3d pt[3], N[3];
        for (unsigned int i = 0; i < 3; ++i)
        {
            pt[i] = C.row(CF(t, i));
        }
        for (unsigned int i = 0; i < 3; ++i)
        {
            N[i] = (pt[(i + 1) % 3] - eta).cross(pt[(i + 2) % 3] - eta);
        }

        for (unsigned int i = 0; i <= 2; ++i)
        {
            w[i] = 0.0;
            for (unsigned int j = 0; j <= 2; ++j)
            {
                w[i] += theta[j] * N[i].dot(N[j]) / (2.0 * N[j].norm());
            }
            w[i] /= determinant;
        }

        sumWeights += (w[0] + w[1] + w[2]);
        w_weights(vid[0]) += w[0];
        w_weights(vid[1]) += w[1];
        w_weights(vid[2]) += w[2];
    }

    for (unsigned int v = 0; v < num_vertices_cage; ++v)
    {
        weights(v) = w_weights(v) / sumWeights;
    }
}

void computeMVC(const Eigen::MatrixXd &C, const Eigen::MatrixXi &CF, Eigen::MatrixXd const &eta_m,
                Eigen::MatrixXd &phi)
{
    phi.resize(C.rows(), eta_m.rows());
    Eigen::VectorXd w_weights(C.rows());
    Eigen::VectorXd weights(C.rows());

    for (int eta_idx = 0; eta_idx < eta_m.rows(); ++eta_idx)
    {
        const Eigen::Vector3d eta = eta_m.row(eta_idx);
        computeMVCForOneVertexSimple(C, CF, eta, weights, w_weights);
        phi.col(eta_idx) = weights;
    }
}

Mesh build_visibility_frustum(const Eigen::Vector3d &v, const Eigen::Vector3d &face_center)
{
    Mesh frustum;

    Eigen::Vector3d dir = face_center - v;
    dir.normalize();

    Eigen::Vector3d up = dir.unitOrthogonal();
    Eigen::Vector3d right = dir.cross(up);

    double distance = 1e5;
    double size = 1e5;
    Eigen::Vector3d far_center = v + dir * distance;

    frustum.V.resize(5, 3);
    frustum.V.row(0) = v;
    frustum.V.row(1) = far_center + size * (up + right);
    frustum.V.row(2) = far_center + size * (-up + right);
    frustum.V.row(3) = far_center + size * (-up - right);
    frustum.V.row(4) = far_center + size * (up - right);

    frustum.F.resize(6, 3);
    frustum.F << 0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 1,
        1, 2, 3,
        1, 3, 4;

    return frustum;
}

Mesh build_closed_prism_around_face(
    const Eigen::MatrixXd &FV,
    double d)
{
    Mesh prism;

    // Compute face normal (convert rows to cols)
    Eigen::Vector3d v0 = FV.row(0).transpose();
    Eigen::Vector3d v1 = FV.row(1).transpose();
    Eigen::Vector3d v2 = FV.row(2).transpose();

    Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0).normalized();

    // Extrude vertices along normal by distance d
    Eigen::MatrixXd top(3, 3);
    for (int i = 0; i < 3; ++i)
        top.row(i) = FV.row(i) + d * normal.transpose();

    // Build vertices: base + top
    prism.V.resize(6, 3);
    prism.V << FV, top;

    prism.F.resize(8, 3);
    int f = 0;

    // Base face
    prism.F.row(f++) << 0, 1, 2;

    // Top face (reverse orientation)
    prism.F.row(f++) << 5, 4, 3;

    for (int i = 0; i < 3; ++i)
    {
        int i_next = (i + 1) % 3;

        prism.F.row(f++) << i, i_next, i_next + 3;
        prism.F.row(f++) << i, i_next + 3, i + 3;
    }

    return prism;
}

Mesh clip_face_along_visibility(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::Vector3d &v_pos,
    int face_idx)
{
    using namespace igl::copyleft::cgal;

    Eigen::MatrixXd FV(3, 3);
    for (int i = 0; i < 3; ++i)
        FV.row(i) = cage_V.row(cage_F(face_idx, i));

    Eigen::MatrixXi FF(1, 3);
    FF << 0, 1, 2;

    Eigen::Vector3d face_center = FV.colwise().mean();
    Mesh frustum = build_visibility_frustum(v_pos, face_center);

    // extrude by a small amount, say 0.01
    Mesh prism = build_closed_prism_around_face(FV, 0.1);

    // Then do boolean with prism instead of FV, FF
    Eigen::MatrixXd CV;
    Eigen::MatrixXi CF;
    igl::copyleft::cgal::mesh_boolean(
        prism.V, prism.F,
        frustum.V, frustum.F,
        igl::MESH_BOOLEAN_TYPE_INTERSECT,
        CV, CF);

    Mesh result;
    result.V = CV;
    result.F = CF;
    return result;
}

void computePMVC_CPU(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::MatrixXd &obj_V,
    Eigen::MatrixXd &pmvc_coords)
{
    pmvc_coords.resize(obj_V.rows(), cage_V.rows());

    for (int vi = 0; vi < obj_V.rows(); ++vi)
    {
        const Eigen::Vector3d v_pos = obj_V.row(vi);
        std::vector<Eigen::MatrixXd> clip_Vs;
        std::vector<Eigen::MatrixXi> clip_Fs;
        int offset = 0;

        for (int fi = 0; fi < cage_F.rows(); ++fi)
        {
            Mesh clipped = clip_face_along_visibility(cage_V, cage_F, v_pos, fi);

            if (clipped.F.rows() > 0)
            {
                clip_Vs.push_back(clipped.V);
                clip_Fs.push_back(clipped.F.array() + offset);
                offset += clipped.V.rows();
            }
        }

        if (clip_Vs.empty())
        {
            pmvc_coords.row(vi).setZero();
            continue;
        }

        Eigen::MatrixXd temp_V(offset, 3);
        Eigen::MatrixXi temp_F;
        {
            int vpos = 0, fpos = 0;
            std::vector<Eigen::MatrixXi> F_all;

            for (size_t i = 0; i < clip_Vs.size(); ++i)
            {
                temp_V.block(vpos, 0, clip_Vs[i].rows(), 3) = clip_Vs[i];
                vpos += clip_Vs[i].rows();
                F_all.push_back(clip_Fs[i]);
                fpos += clip_Fs[i].rows();
            }
            F_all.push_back(temp_F);
        }

        Eigen::VectorXd weights;
        Eigen::VectorXd w_weights;
        computeMVCForOneVertexSimple(temp_V, temp_F, v_pos, weights, w_weights);
        pmvc_coords.row(vi) = weights.transpose();
    }
}

void computePMVC_GPU_lipman(
    const Eigen::MatrixXd &cage_V,
    const Eigen::MatrixXi &cage_F,
    const Eigen::MatrixXd &obj_V,
    Eigen::MatrixXd &pmvc_coords)
{
}