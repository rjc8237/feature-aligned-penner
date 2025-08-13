#include <igl/readOBJ.h>
#include "feature/interface.h"
#include "feature/core/io.h"
#include "holonomy/field/frame_field.h"
#include "util/io.h"
#include "util.h"
#include "holonomy/core/viewer.h"
#include "feature/surgery/cut_metric_generator.h"

#include <CLI/CLI.hpp>
#include <igl/bounding_box_diagonal.h>


using namespace Penner;
using namespace Penner::Optimization;
using namespace Penner::Holonomy;
using namespace Penner::Feature;

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);

    // Get command line arguments
    CLI::App app{"Generate a feature aligned parametrization."};
    std::string mesh = "";
    std::string input_dir = "./";
    std::string output_dir = "./";

    // IO Parameters
    bool use_existing_field = false;
    bool show_parameterization = false;
    app.add_option("--name", mesh, "Mesh name (without obj suffix, e.g., fandisk)")->required();
    app.add_option("-i,--input", input_dir, "Input directory")->check(CLI::ExistingDirectory)->required();
    app.add_option("-o,--output", output_dir, "Output directory");
    app.add_flag("--use_existing_field", use_existing_field, "Use precomputed field at the input directory");
    app.add_flag("--show_parameterization", show_parameterization, "Show aligned parameterization");

    // Marked Metric Parameters
    NewtonParameters alg_params;
    add_newton_parameters(app, alg_params);
    CLI11_PARSE(app, argc, argv);

    std::filesystem::create_directory(output_dir);

    // create filepaths for input data
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
    Eigen::MatrixXd reference_field;
    Eigen::VectorXd theta;
    Eigen::MatrixXd kappa;
    Eigen::MatrixXi period_jump;
    if (use_existing_field)
    {
        spdlog::info("loading feature edges");
        feature_edges = load_feature_edges(feature_filename);
        hard_feature_edges = load_feature_edges(hard_feature_filename);

        spdlog::info("loading constraints");
        std::tie(reference_field, theta, kappa, period_jump) = load_frame_field(field_filename);
    }
    else {
        // refine input mesh
        std::tie(V, F, feature_edges, hard_feature_edges) = generate_refined_feature_mesh(V, F, false);
        FeatureFinder feature_finder(V, F);
        feature_finder.mark_features(feature_edges);
        auto[V_cut, F_cut, V_map, F_is_feature] = feature_finder.generate_feature_cut_mesh();

        int radius = 5;
        Scalar rel_anisotropy=0.9;
        Scalar abs_anisotropy=0.2;
        Scalar bb_diag = igl::bounding_box_diagonal(V);
        auto [direction, is_fixed_direction] = compute_field_direction(
            V_cut,
            F_cut,
            radius,
            abs_anisotropy / bb_diag,
            rel_anisotropy);
        MarkedMetricParameters marked_metric_params;
        marked_metric_params.remove_trivial_torus = false; // FIXME
        marked_metric_params.use_log_length = true;
        marked_metric_params.use_initial_zero = false;
        CutMetricGenerator cut_metric_generator(V_cut, F_cut, marked_metric_params, {});
        cut_metric_generator.generate_fields(V_cut, F_cut, V_map, direction, is_fixed_direction);
        std::tie(reference_field, theta, kappa, period_jump) = cut_metric_generator.get_field();
    }

    // get optimized metric
    spdlog::info("projecting to feature constraints");
    alg_params.output_dir = output_dir;
    alg_params.error_eps = 1e-10;
    alg_params.solver = "ldlt";
    MarkedMetricParameters marked_metric_params;
    AlignedMetricGenerator aligned_metric_generator(
        V,
        F,
        feature_edges,
        hard_feature_edges,
        reference_field,
        theta,
        kappa,
        period_jump,
        marked_metric_params);
    aligned_metric_generator.optimize_relaxed(alg_params);
    aligned_metric_generator.parameterize(false);
    auto [V_r, F_r, uv_r, FT_r, fn_to_f_r, endpoints_r] = aligned_metric_generator.get_parameterization();

    if (show_parameterization) view_seamless_parameterization(V_r, F_r, uv_r, FT_r, "refined mesh", true);

    std::string output_filename = join_path(output_dir, mesh+"_opt.obj");
    write_obj_with_uv(output_filename, V_r, F_r, uv_r, FT_r);

    //std::string output_filename = join_path(output_dir, "optimized_corner_coords");
    //write_matrix(opt_corner_coords, output_filename, " ");

}
