################################################################################
# run_sams_plot.jl
#
# Purpose:
#   Run the Bayesian hierarchical model on the SAMS data and 
#   create plot for comparing Bayesian and MatrixLM estimates for one covariate.
################################################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "gibbs_src", "construct_sams_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "posterior_summaries.jl"))
include(joinpath(@__DIR__, "gibbs_src", "make_plots.jl"))

mkpath(joinpath(@__DIR__, "results", "figures"))

println("Constructing SAMS matrices...")
obj = construct_sams_matrices()

println("Fitting MatrixLM...")
mlm_res = fit_copdgene_matrixlm(obj.X, obj.Y, obj.Z)

cov_choose = "Fish Oil"

println("Selecting target covariate: ", cov_choose)
k_cov = obj.coef_lookup[cov_choose]

b_vec = vec(mlm_res.coef[k_cov, :])
se_vec = vec(mlm_res.se[k_cov, :])

println("Running Bayesian hierarchical model...")

subclass_of_met = obj.subclass_of_met
H = obj.H

seeds = [101, 102, 103, 104]

res_chains = [
    gibbs_meta_hier_traces_one(
        b_vec,
        se_vec,
        subclass_of_met,
        H;
        mu0 = 0.0,
        s0 = 1.0,
        halfcauchy_scale = 1.0,
        n_iter = 5000,
        burnin = 1000,
        thin = 1,
        seed = s
    )
    for s in seeds
]

println("Creating posterior summaries...")

cmp_df = make_comparison_dataframe(
    res_chains,
    b_vec,
    se_vec
)

println("Creating plots...")

plot_name = "sams_" * replace(cov_choose, " " => "_") * "_comparison.png"

cmp_plot_path = joinpath(
    @__DIR__,
    "results",
    "figures",
    plot_name
)

make_comparison_plot(cmp_df, cmp_plot_path)

println("Saved comparison plot to: ", cmp_plot_path)

println("Done.")