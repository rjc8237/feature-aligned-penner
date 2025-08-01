#include "holonomy/core/common.h"
#include "feature/core/viewer.h"
#include "feature/core/quads.h"
#include "util/io.h"

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

using namespace Penner;
using namespace Penner::Feature;

int main(int argc, char* argv[])
{
    // Get command line arguments
    CLI::App app{"View a quad mesh"};
    std::string mesh_filename = "";
    std::string output_filename = "./quad_with_uv.obj";

    // IO Parameters
    app.add_option("--mesh", mesh_filename, "Mesh filepath")->check(CLI::ExistingFile)->required();
    app.add_option("--output", output_filename, "Output filepath");
    CLI11_PARSE(app, argc, argv);

    spdlog::set_level(spdlog::level::debug);

    // Get input mesh
    Eigen::MatrixXd V, _uv, N;
    Eigen::MatrixXi F, _FT, FN;
    spdlog::info("Using mesh at {}", mesh_filename);
    igl::readOBJ(mesh_filename, V, _uv, N, F, _FT, FN);

    int num_faces = F.rows();
    Eigen::MatrixXd uv(4, 2);
    Eigen::MatrixXi FT(num_faces, 4);
    uv.row(0) << 0., 0.;
    uv.row(1) << 1., 0.;
    uv.row(2) << 1., 1.;
    uv.row(3) << 0., 1.;
    for (int f = 0; f < num_faces; ++f)
    {
        FT.row(f) << 0, 1, 2, 3;
    }
    igl::writeOBJ(output_filename, V, F, N, FN, uv, FT);

    //view_quad_mesh(V, F);
}
