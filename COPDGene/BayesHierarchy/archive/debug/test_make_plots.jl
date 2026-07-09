using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "construct_matrices.jl"))
include(joinpath(@__DIR__, "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "posterior_summaries.jl"))
include(joinpath(@__DIR__, "make_plots.jl"))

obj = construct_copdgene_matrices()
mlm_res = fit_copdgene_matrixlm(obj.X, obj.Y, obj.Z)

k_age = findfirst(==("Age"), obj.coef_names)

b_vec = vec(mlm_res.coef[k_age, :])
se_vec = vec(mlm_res.se[k_age, :])

m = length(b_vec)
subclass_of_met = ones(Int, m)
H = 1

seeds = [101, 102]

res_chains = [
    gibbs_meta_hier_traces_one(
        b_vec,
        se_vec,
        subclass_of_met,
        H;
        n_iter = 1000,
        burnin = 200,
        thin = 1,
        seed = s
    )
    for s in seeds
]

cmp_df = make_comparison_dataframe(
    res_chains,
    b_vec,
    se_vec
)

mkpath(joinpath(@__DIR__, "results", "figures"))

make_se_comparison_plot(
    cmp_df,
    joinpath(@__DIR__, "results", "figures", "test_se_comparison.png")
)

make_estimate_comparison_plot(
    cmp_df,
    joinpath(@__DIR__, "results", "figures", "test_estimate_comparison.png")
)

println("Plots saved in results/figures/")