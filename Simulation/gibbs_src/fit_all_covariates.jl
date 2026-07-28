############################################################
# fit_all_covariates.jl
#
# Purpose:
#   Fit the one-level Bayesian hierarchical model for all 
#   covariates using MatrixLM estimates and standard errors.  
############################################################

using Statistics

function summarize_theta_draws(theta_draws::AbstractMatrix{<:Real})
    m, n_keep = size(theta_draws)
    mean_theta = vec(mean(theta_draws; dims=2))
    sd_theta   = vec(std(theta_draws; dims=2, corrected=true))
    return mean_theta, sd_theta
end

function fit_bayes_all_covariates_one(
    B_obs::AbstractMatrix{<:Real},
    SE_obs::AbstractMatrix{<:Real},
    subclass_of_met::Vector{Int},
    H::Int;
    mu0::Float64 = 0.0,
    s0::Float64 = 1.0,
    halfcauchy_scale::Float64 = 1.0,
    n_iter::Int = 5000,
    burnin::Int = 1000,
    thin::Int = 1,
    seed0::Int = 42,
    keep_results::Bool = false
)
    p, m = size(B_obs)
    @assert size(SE_obs) == (p, m) "SE_obs must have same shape as B_obs"

    B_bayes  = Array{Float64}(undef, p, m)
    SE_bayes = Array{Float64}(undef, p, m)

    res_list = keep_results ? Vector{Any}(undef, p) : Any[]

    for k in 1:p
        b_vec  = vec(Float64.(B_obs[k, :]))
        se_vec = vec(Float64.(SE_obs[k, :]))

        res_k = gibbs_meta_hier_traces_one(
            b_vec,
            se_vec,
            subclass_of_met,
            H;
            mu0 = mu0,
            s0 = s0,
            halfcauchy_scale = halfcauchy_scale,
            n_iter = n_iter,
            burnin = burnin,
            thin = thin,
            seed = seed0 + k
        )

        @assert res_k.theta_draws !== nothing "theta_draws missing"

        mean_theta, sd_theta = summarize_theta_draws(res_k.theta_draws)

        @inbounds begin
            B_bayes[k, :]  .= mean_theta
            SE_bayes[k, :] .= sd_theta
        end

        if keep_results
            res_list[k] = res_k
        end

    end

    return B_bayes, SE_bayes, res_list
end
