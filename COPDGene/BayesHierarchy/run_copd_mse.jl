############################################################################
# run_copd_mse.jl
#
# Reproducible repeated-subsampling analysis for the COPDGene
# application. Produces covariate-specific reference MSE metrics.
#
# Usage:
#   julia --project=. run_copd_mse.jl       # R = 100
#   julia --project=. run_copd_mse.jl 1     # quick test
############################################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames

include(joinpath(@__DIR__, "gibbs_src", "construct_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_all_covariates.jl"))
include(joinpath(@__DIR__, "gibbs_src", "subsampling_analysis.jl"))

const TRAIN_SIZES = [
    40, 60, 80, 100, 125, 150,
    200, 250, 300, 400, 500, 600
]

const DEFAULT_REPETITIONS = 100
const SPLIT_SEED = 2026
const GIBBS_SEED = 1000

function parse_repetitions(args)
    R = isempty(args) ? DEFAULT_REPETITIONS : parse(Int, first(args))
    @assert R > 0 "Number of repetitions must be positive"
    return R
end

R = parse_repetitions(ARGS)

table_dir = joinpath(@__DIR__, "results", "tables")
mkpath(table_dir)

println("Constructing COPDGene matrices...")
obj = construct_copdgene_matrices()

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

covariates_of_interest = ["Age", "BMI", "COPD: 1"]

missing_covariates = setdiff(
    covariates_of_interest,
    unique(summary.covariate)
)

@assert isempty(missing_covariates) "Missing covariates: $missing_covariates"

selected_covariates = subset(
    summary,
    :covariate => ByRow(x -> x in covariates_of_interest)
)

selected_covariates = select(
    selected_covariates,
    :covariate,
    :n_train,
    :n_test,
    :mean_mse_mlm,
    :mean_mse_bayes,
    :ratio_of_mean_mses,
    :percent_reduction_from_mean_mses,
    :R
)

sort!(selected_covariates, [:covariate, :n_train])

all_reps_path = joinpath(
    table_dir,
    "copdgene_rmr_all_reps.csv"
)

summary_path = joinpath(
    table_dir,
    "copdgene_rmr_summary.csv"
)

selected_path = joinpath(
    table_dir,
    "copdgene_rmr_selected_covariates.csv"
)

CSV.write(all_reps_path, repetition_metrics)
CSV.write(summary_path, summary)
CSV.write(selected_path, selected_covariates)

println()
println("Saved:")
println(all_reps_path)
println(summary_path)
println(selected_path)
println("Done.")
