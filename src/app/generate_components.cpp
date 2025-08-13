#include "holonomy/core/common.h"
#include "holonomy/core/viewer.h"
#include "feature/feature/error.h"
#include "feature/core/io.h"
#include "feature/surgery/cut_mesh_layout.h"
#include "feature/interface.h"

#include "util/vf_mesh.h"

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include <CLI/CLI.hpp>
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

using namespace Penner;
using namespace Penner::Holonomy;
using namespace Penner::Feature;



int main(int argc, char* argv[])
{
    // Get command line arguments
    CLI::App app{"View a quad mesh"};
    std::string mesh_filename = "";

    // IO Parameters
    app.add_option("--mesh", mesh_filename, "Mesh filepath")->check(CLI::ExistingFile)->required();
    CLI11_PARSE(app, argc, argv);

    spdlog::set_level(spdlog::level::debug);

    // Get input mesh
    Eigen::MatrixXd V, uv, N;
    Eigen::MatrixXi F, FT, FN;
    spdlog::info("Using mesh at {}", mesh_filename);
    igl::readOBJ(mesh_filename, V, uv, N, F, FT, FN);

    std::vector<VertexEdge> E = load_mesh_edges(mesh_filename);

    Eigen::MatrixXi F_is_seam = find_seams(F, FT);
    auto [V_seams, E_seams] = generate_edges(V, F, F_is_seam);

    // get components of layout
    Eigen::VectorXi components;
    igl::facet_components(FT, components);
    int num_components = components.maxCoeff() + 1;
    std::vector<std::vector<int>> component_faces(num_components, std::vector<int>({}));
    int num_faces = F.rows();
    for (int fi = 0; fi < num_faces; ++fi)
    {
        int ci = components[fi];
        component_faces[ci].push_back(fi);
    }

    for (int ci = 0; ci < num_components; ++ci)
    {
        int num_component_faces = component_faces[ci].size();
        Eigen::MatrixXi Fc(num_component_faces, 3);
        Eigen::MatrixXi FTc(num_component_faces, 3);
        for (int i = 0; i < num_component_faces; ++i)
        {
            int fi = component_faces[ci][i];
            Fc.row(i) = F.row(fi);
            FTc.row(i) = FT.row(fi);
        }
        Eigen::MatrixXi NF, NFT, NFN, TF;
        Eigen::MatrixXd NV, Nuv, NN, TV;
        Eigen::VectorXi I, J;
        igl::remove_unreferenced(V, Fc, NV, NF, I, J);
        igl::remove_unreferenced(uv, FTc, Nuv, NFT, I, J);

        int num_component_vertices = Nuv.rows();
        Eigen::MatrixXd uv_embed = Eigen::MatrixXd::Zero(num_component_vertices, 3);
        for (int vi = 0; vi < num_component_vertices; ++vi)
        {
            uv_embed(vi, 0) = Nuv(vi, 0);
            uv_embed(vi, 1) = Nuv(vi, 1);
        }

        spdlog::info("{}", Nuv);
        std::string output_filename = "component_" + std::to_string(ci) + ".obj";
        igl::writeOBJ(output_filename, NV, NF, NN, NFN, Nuv, NFT);
        output_filename = "layout_" + std::to_string(ci) + ".obj";
        igl::writeOBJ(output_filename, uv_embed, NFT, NN, NFN, Nuv, NFT);
    }
}
