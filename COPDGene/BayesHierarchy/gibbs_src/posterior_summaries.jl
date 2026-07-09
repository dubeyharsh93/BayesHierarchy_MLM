############################################################
# posterior_summaries.jl
#
# Purpose:
#   Summarize Bayesian posterior draws and compare them
#   with MatrixLM estimates and standard errors.
############################################################

using DataFrames
using Statistics
using StatsBase

function combine_theta_draws(res_chains)
    return hcat([r.theta_draws for r in res_chains]...)
end

function combine_beta_draws(res_chains)
    return hcat([r.beta_draws for r in res_chains]...)
end

function summarize_draw_matrix(draws::AbstractMatrix)
    n_params = size(draws, 1)

    posterior_mean = Vector{Float64}(undef, n_params)
    posterior_sd = Vector{Float64}(undef, n_params)
    lower_95 = Vector{Float64}(undef, n_params)
    upper_95 = Vector{Float64}(undef, n_params)

    for j in 1:n_params
        x = vec(draws[j, :])
        posterior_mean[j] = mean(x)
        posterior_sd[j] = std(x)
        lower_95[j] = quantile(x, 0.025)
        upper_95[j] = quantile(x, 0.975)
    end

    return DataFrame(
        index = 1:n_params,
        posterior_mean = posterior_mean,
        posterior_sd = posterior_sd,
        lower_95 = lower_95,
        upper_95 = upper_95
    )
end

function make_theta_summary(res_chains)
    theta_all = combine_theta_draws(res_chains)
    return summarize_draw_matrix(theta_all)
end

function make_beta_summary(res_chains)
    beta_all = combine_beta_draws(res_chains)
    return summarize_draw_matrix(beta_all)
end

function make_comparison_dataframe(
    res_chains,
    b_obs::Vector{Float64},
    se_obs::Vector{Float64};
    metabolite_names = nothing
)

    theta_summary = make_theta_summary(res_chains)

    m = length(b_obs)

    if metabolite_names === nothing
        metabolite_names = string.("met_", 1:m)
    end

    df = DataFrame(
        metabolite_index = 1:m,
        metabolite = metabolite_names,
        matrixlm_estimate = b_obs,
        matrixlm_se = se_obs,
        bayes_estimate = theta_summary.posterior_mean,
        bayes_se = theta_summary.posterior_sd,
        bayes_lower_95 = theta_summary.lower_95,
        bayes_upper_95 = theta_summary.upper_95
    )

    df.abs_matrixlm_estimate = abs.(df.matrixlm_estimate)
    df.abs_bayes_estimate = abs.(df.bayes_estimate)
    df.se_ratio_matrixlm_to_bayes = df.matrixlm_se ./ df.bayes_se

    return df
end