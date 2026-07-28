####################################################################
# run_sams_mse.jl
#
# Reproducible repeated-subsampling analysis for the SAMS study.
# Produces the manuscript table for the interaction between SAMS
# status and fish-oil supplementation.
#
# Usage:
#   julia --project=. run_sams_mse.jl       # R = 100
#   julia --project=. run_sams_mse.jl 1     # quick test
####################################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Printf

include(joinpath(@__DIR__, "gibbs_src", "construct_sams_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_all_covariates.jl"))
include(joinpath(@__DIR__, "gibbs_src", "subsampling_analysis.jl"))

const TRAIN_SIZES = [20, 30, 40, 50, 60]
const DEFAULT_REPETITIONS = 100
const SPLIT_SEED = 2026
const GIBBS_SEED = 1000
const FOCAL_COVARIATE = "Interaction SAMS-Fish Oil"

function parse_repetitions(args)
    R = isempty(args) ? DEFAULT_REPETITIONS : parse(Int, first(args))
    @assert R > 0 "Number of repetitions must be positive"
    return R
end

function latex_booktabs_table(df::DataFrame)
    io = IOBuffer()

    println(io, raw"\begin{tabular}{cccccc}")
    println(io, raw"\toprule")
    println(
        io,
        "Training Size & Reference Size & MatrixLM MSE & " *
        "Bayesian MSE & Reference MSE Ratio & " *
        "Error Reduction (\\%) \\\\"
    )
    println(io, raw"\midrule")

    for row in eachrow(df)
        @printf(
            io,
            "%d & %d & %.4f & %.4f & %.4f & %.1f \\\\\n",
            row.n_train,
            row.n_test,
            row.mean_mse_mlm,
            row.mean_mse_bayes,
            row.ratio_of_mean_mses,
            row.percent_reduction_from_mean_mses
        )
    end

    println(io, raw"\bottomrule")
    println(io, raw"\end{tabular}")

    return String(take!(io))
end

R = parse_repetitions(ARGS)

table_dir = joinpath(@__DIR__, "results", "tables")
mkpath(table_dir)

println("Constructing SAMS matrices...")
obj = construct_sams_matrices()

println(
    "Preparing repeated subsamples for $(length(TRAIN_SIZES)) " *
    "training sizes with R=$R..."
)

splits = make_subsample_splits(
    obj.Y,
    obj.X;
    train_sizes = TRAIN_SIZES,
    R = R,
    seed = SPLIT_SEED
)

fit_results = fit_repeated_subsamples(
    splits,
    obj.Y,
    obj.X,
    obj.Z,
    obj.subclass_of_met,
    obj.H;
    mu0 = 0.0,
    s0 = 1.0,
    halfcauchy_scale = 1.0,
    n_iter = 5000,
    burnin = 1000,
    thin = 1,
    seed0 = GIBBS_SEED,
    keep_results = false
)

repetition_metrics, summary = summarize_reference_agreement(
    fit_results;
    covariate_names = String.(obj.coef_names)
)

@assert FOCAL_COVARIATE in summary.covariate "Missing focal covariate: $FOCAL_COVARIATE"

interaction_summary = subset(
    summary,
    :covariate => ByRow(==(FOCAL_COVARIATE))
)

interaction_table = select(
    interaction_summary,
    :n_train,
    :n_test,
    :mean_mse_mlm,
    :mean_mse_bayes,
    :ratio_of_mean_mses,
    :percent_reduction_from_mean_mses
)

sort!(interaction_table, :n_train)

all_reps_path = joinpath(
    table_dir,
    "sams_rmr_all_reps.csv"
)

summary_path = joinpath(
    table_dir,
    "sams_rmr_summary.csv"
)

interaction_path = joinpath(
    table_dir,
    "sams_rmr_interaction_fish_oil.csv"
)

CSV.write(all_reps_path, repetition_metrics)
CSV.write(summary_path, summary)
CSV.write(interaction_path, interaction_table)

println()
println("SAMS results: $FOCAL_COVARIATE")
println(latex_booktabs_table(interaction_table))

println("Saved:")
println(all_reps_path)
println(summary_path)
println(interaction_path)
println("Done.")
