############################################################
# gibbs_sampler.jl
#
# Purpose:
#   Gibbs sampler for the one-level Bayesian hierarchical
#   model using MatrixLM estimates and standard errors.
############################################################

using Random
using Distributions
using DataFrames
using Statistics


# Half-Cauchy distribution truncated to positive values
half_cauchy(s) = truncated(Cauchy(0.0, s), 0.0, Inf)

# Function for normal update from precision (for conjugate updates)
@inline function normal_update_from_prec(lik_prec, lik_sum, prior_prec, prior_mean)
    post_prec = lik_prec + prior_prec
    post_mean = (lik_sum + prior_prec * prior_mean) / post_prec
    post_var  = 1.0 / post_prec
    return post_mean, post_var
end

# Function to sample from Inverse-Gamma distribution
@inline ig_sample(shape, scale) = rand(InverseGamma(shape, scale))

# Gibbs Sampler for one level hierarchy 
function gibbs_meta_hier_traces_one(
    b_obs::Vector{Float64},           # length m
    se_obs::Vector{Float64},          # length m
    subclass_of_met::Vector{Int},     # length m, values in 1..H (Total DB class)
    H::Int;                           # number of DB categories
    mu0::Float64 = 0.0,
    s0::Float64 = 1.0,
    halfcauchy_scale::Float64 = 1.0,
    n_iter::Int = 2000,
    burnin::Int = 500,
    thin::Int = 1,
    seed::Int = 1234,
)
    @assert length(b_obs) == length(se_obs) "b_obs and se_obs must have same length"
    @assert length(subclass_of_met) == length(b_obs) "subclass_of_met must have length m"
    @assert maximum(subclass_of_met) == H "subclass_of_met must be in 1..H"
    @assert burnin < n_iter "burnin must be < n_iter"
    @assert thin ≥ 1

    Random.seed!(seed)
    m = length(b_obs)
    lam_j = 1.0 ./ (se_obs .^ 2)  # observation precisions

    # --- initialize state ---
    theta = copy(b_obs)                                     # length m (metabolite-level effects)
    beta  = [mean(theta[subclass_of_met .== h]) for h in 1:H]  # class means
    theta0 = mean(beta)                                     # global mean over classes

    tau_w2 = 1.0   # within-class variance for theta_j | beta_h
    tau_v2 = 1.0   # between-class variance for beta_h | theta0

    lambda_w = 1.0    # IG mixture auxiliaries for half-Cauchy
    lambda_v = 1.0

    # --- precompute group indices ---
    idx_by_sub = [findall(==(h), subclass_of_met) for h in 1:H]

    # --- storage sizes ---
    n_keep = floor(Int, (n_iter - burnin) ÷ thin)

    # scalars per draw (for quick monitoring)
    draws_scalar = DataFrame(
        theta0 = Vector{Float64}(undef, n_keep),
        tau_w2 = Vector{Float64}(undef, n_keep),
        tau_v2 = Vector{Float64}(undef, n_keep),
    )

    # vectors per draw
    beta_draws  = Array{Float64}(undef, H, n_keep)   # columns = kept iters
    theta_draws = Array{Float64}(undef, m, n_keep)

    keep_idx = 0

    for it in 1:n_iter
        #######################
        # 1) update theta_j   #
        #######################
        for j in 1:m
            h = subclass_of_met[j]
            lik_prec   = lam_j[j]
            lik_sum    = lam_j[j] * b_obs[j]
            prior_prec = 1.0 / tau_w2
            prior_mean = beta[h]
            mean_th = (lik_sum + prior_prec * prior_mean) / (lik_prec + prior_prec)
            var_th  = 1.0 / (lik_prec + prior_prec)
            theta[j] = rand(Normal(mean_th, sqrt(var_th)))
        end

        #######################
        # 2) update beta_h    #
        #######################
        for h in 1:H
            J = idx_by_sub[h]
            mh = length(J)

            # likelihood: theta_j | beta_h ~ N(beta_h, tau_w2)
            lik_prec   = mh / tau_w2
            lik_sum    = (1.0 / tau_w2) * sum(theta[J])

            # prior: beta_h | theta0 ~ N(theta0, tau_v2)
            prior_prec = 1.0 / tau_v2
            prior_mean = theta0

            mean_b = (lik_sum + prior_prec * prior_mean) / (lik_prec + prior_prec)
            var_b  = 1.0 / (lik_prec + prior_prec)
            beta[h] = rand(Normal(mean_b, sqrt(var_b)))
        end

        #######################
        # 3) update theta0    #
        #######################
        # beta_h | theta0 ~ N(theta0, tau_v2)
        lik_prec   = H / tau_v2
        lik_sum    = (1.0 / tau_v2) * sum(beta)

        # prior: theta0 ~ N(mu0, s0^2)
        prior_prec = 1.0 / (s0^2)
        prior_mean = mu0

        mean_t0 = (lik_sum + prior_prec * prior_mean) / (lik_prec + prior_prec)
        var_t0  = 1.0 / (lik_prec + prior_prec)
        theta0  = rand(Normal(mean_t0, sqrt(var_t0)))

        ############################################
        # 4) update tau_w2 and tau_v2 via IG mix   #
        ############################################

        # tau_w2: within-class spread of theta_j around beta_{class(j)}
        ssw = sum((theta .- beta[subclass_of_met]).^2)
        tau_w2   = rand(InverseGamma((m + 1.0) / 2.0, 0.5 * ssw + 1.0 / lambda_w))
        lambda_w = rand(InverseGamma(1.0, 1.0 / (halfcauchy_scale^2) + 1.0 / tau_w2))

        # tau_v2: between-class spread of beta_h around theta0
        ssv = sum((beta .- theta0).^2)
        tau_v2   = rand(InverseGamma((H + 1.0) / 2.0, 0.5 * ssv + 1.0 / lambda_v))
        lambda_v = rand(InverseGamma(1.0, 1.0 / (halfcauchy_scale^2) + 1.0 / tau_v2))

        ###################################
        # 5) store post-burn (with thin) #
        ###################################
        if it > burnin && ((it - burnin) % thin == 0)
            keep_idx += 1
            draws_scalar.theta0[keep_idx] = theta0
            draws_scalar.tau_w2[keep_idx] = tau_w2
            draws_scalar.tau_v2[keep_idx] = tau_v2

            beta_draws[:, keep_idx]  .= beta
            theta_draws[:, keep_idx] .= theta
        end
    end

    return (
        draws_scalar = draws_scalar,
        beta_draws = beta_draws,
        theta_draws = theta_draws,
        last_state = (
            theta0 = theta0,
            tau_w2 = tau_w2,
            tau_v2 = tau_v2,
            beta = copy(beta),
            theta = copy(theta),
        )
    )

end