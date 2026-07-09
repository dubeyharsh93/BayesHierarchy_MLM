using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "construct_matrices.jl"))
include(joinpath(@__DIR__, "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_sampler.jl"))

obj = construct_copdgene_matrices()
mlm_res = fit_copdgene_matrixlm(obj.X, obj.Y, obj.Z)

# Choose Age coefficient
k_age = findfirst(==("Age"), obj.coef_names)

b_vec = vec(mlm_res.coef[k_age, :])
se_vec = vec(mlm_res.se[k_age, :])

# Temporary grouping for testing:
# one group for all metabolites.
# Later we will replace this with DB subclass grouping.
m = length(b_vec)
subclass_of_met = ones(Int, m)
H = 1

res = gibbs_meta_hier_traces_one(
    b_vec,
    se_vec,
    subclass_of_met,
    H;
    mu0 = 0.0,
    s0 = 1.0,
    halfcauchy_scale = 1.0,
    n_iter = 1000,
    burnin = 200,
    thin = 1,
    seed = 42
)

println("theta_draws size: ", size(res.theta_draws))
println("beta_draws size: ", size(res.beta_draws))
println("draws_scalar size: ", size(res.draws_scalar))
println("Posterior mean theta0: ", mean(res.draws_scalar.theta0))
println("Posterior mean tau_w2: ", mean(res.draws_scalar.tau_w2))
println("Posterior mean tau_v2: ", mean(res.draws_scalar.tau_v2))