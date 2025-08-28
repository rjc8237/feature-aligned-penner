#include "holonomy/core/common.h"
#include "feature/core/viewer.h"
#include "feature/interface.h"
#include "feature/feature/features.h"
#include "feature/surgery/cut_metric_generator.h"
#include "feature/core/io.h"
#include "holonomy/field/frame_field.h"
#include "holonomy/core/viewer.h"

#include <igl/readOBJ.h>
#include <CLI/CLI.hpp>

using namespace Penner;
using namespace Penner::Feature;
using namespace Penner::Holonomy;

int main(int argc, char* argv[])
{
    // Get command line arguments
    CLI::App app{"Modify frame field in interactive viewer"};
    std::string mesh = "";
    std::string input_dir = "./";

    // IO Parameters
    app.add_option("--name", mesh, "Mesh name (without obj suffix, e.g., fandisk)")->required();
    app.add_option("-i,--input", input_dir, "Input directory")->check(CLI::ExistingDirectory)->required();
    CLI11_PARSE(app, argc, argv);

    spdlog::set_level(spdlog::level::debug);

    std::string mesh_filename = join_path(input_dir, mesh + ".obj");
    std::string feature_filename = join_path(input_dir, mesh + "_features");
    std::string hard_feature_filename = join_path(input_dir, mesh + "_hard_features");
    std::string field_filename = join_path(input_dir, mesh + ".ffield");

    // Get input mesh
    Eigen::MatrixXd V, uv, N;
    Eigen::MatrixXi F, FT, FN;
    spdlog::info("optimizing mesh at {}", mesh_filename);
    igl::readOBJ(mesh_filename, V, uv, N, F, FT, FN);

    // Get features and field
    std::vector<VertexEdge> feature_edges, hard_feature_edges;
    feature_edges = load_feature_edges(feature_filename);
    hard_feature_edges = load_feature_edges(hard_feature_filename);
    auto [reference_field, theta, kappa, period_jump] = load_frame_field(field_filename);

    FeatureFinder feature_finder(V, F);
    feature_finder.mark_features(feature_edges);
    auto [V_cut, F_cut, V_map, F_is_feature] = feature_finder.generate_feature_cut_mesh();

    // generate cut metric with fields
    MarkedMetricParameters marked_metric_params;
    marked_metric_params.use_log_length = true;
    CutMetricGenerator cut_metric_generator(V_cut, F_cut, marked_metric_params, {});
    cut_metric_generator.set_fields(F_cut, reference_field, theta, kappa, period_jump);
    auto [marked_metric, vtx_reindex, face_reindex, rotation_form, Th_hat] = cut_metric_generator.get_union_metric( marked_metric_params);

    // view field
    IntrinsicNRosyField field_generator;
    field_generator.initialize(marked_metric);
    field_generator.set_field(marked_metric, vtx_reindex, F_cut, face_reindex, theta, kappa, period_jump);
    Eigen::MatrixXd V_disp = displace_cut_faces(V_cut, F_cut);
    field_generator.collapse_adjacent_cones(marked_metric);
    field_generator.view(marked_metric, vtx_reindex, V_disp);

    view_feature_cross_field(V_cut, F_cut, reference_field, theta, kappa, period_jump);
}
