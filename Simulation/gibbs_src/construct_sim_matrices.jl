
############################################################
# construct_sim_matrices.jl
#
# Purpose:
#   Construct simulation matrices for Bayesian hierarchical 
#   model.
############################################################

using Random
using Distributions
using LinearAlgebra
using Statistics
using DataFrames

function make_sim_coef_lookup(coef_names::Vector{String})
    return Dict(nm => i for (i, nm) in enumerate(coef_names))
end

"""
simulate_bilinear_identity_1level

Simulates X, B_true, and Y under

    Y = X * B_true + E

with a one-level subclass hierarchy across metabolites.
"""

function simulate_bilinear_identity_1level(;
    n::Int = 98,
    m::Int = 770,
    p::Int = 10,
    H::Int = 4,
    prop_sub = fill(1/H, H),
    theta0_scale::Real = 0.2f0,
    tau_v::Real = 0.15f0,
    tau_w::Real = 0.10f0,
    sigma_y::Real = 1.0f0,
    seed::Int = 123,
    T = Float32
)
    @assert length(prop_sub) == H
    @assert abs(sum(prop_sub) - 1) < 1e-6

    Random.seed!(seed)

    theta0_scaleT = T(theta0_scale)
    tau_vT = T(tau_v)
    tau_wT = T(tau_w)
    sigma_yT = T(sigma_y)

    X = randn(T, n, p)

    counts = round.(Int, m .* collect(prop_sub))
    counts[end] += m - sum(counts)

    subclass_of_met = Vector{Int}(undef, m)
    idx = 1
    for h in 1:H
        for _ in 1:counts[h]
            subclass_of_met[idx] = h
            idx += 1
        end
    end

    perm = randperm(m)
    subclass_of_met = subclass_of_met[perm]

    theta0 = rand.(Normal(T(0), theta0_scaleT), p)

    beta = Matrix{T}(undef, H, p)
    for h in 1:H, k in 1:p
        beta[h, k] = rand(Normal(theta0[k], tau_vT))
    end

    theta = Matrix{T}(undef, m, p)
    for j in 1:m
        h = subclass_of_met[j]
        for k in 1:p
            theta[j, k] = rand(Normal(beta[h, k], tau_wT))
        end
    end

    B_true = permutedims(theta)

    Y_mean = X * B_true
    Y = Y_mean .+ randn(T, n, m) .* sigma_yT

    Z = Matrix{Float64}(I, m, m)
    coef_names = ["x$(k)" for k in 1:p]
    coef_lookup = Dict(nm => i for (i, nm) in enumerate(coef_names))

    return (
        X = Matrix{Float64}(X),
        Y = Matrix{Float64}(Y),
        Z = Z,
        B_true = Matrix{Float64}(B_true),
        theta0 = theta0,
        beta = beta,
        subclass_of_met = Int.(subclass_of_met),
        perm = perm,
        H = H,
        coef_names = coef_names,
        pseudo_coef_names = coef_names,
        coef_lookup = coef_lookup
    )
end

function simulate_bilinear_identity_1level2(;
    n::Int,
    m::Int,
    p::Int,
    H::Int,
    prop_sub::Vector,
    theta0_scale = 0.2f0,
    tau_v = 0.12f0,
    tau_w = 0.08f0,
    sigma_y = 1.0f0,
    T = Float32,
    seed::Int = 1
)

    @assert length(prop_sub) == H
    @assert isapprox(sum(prop_sub), 1.0; atol = 1e-6)

    rng = MersenneTwister(seed)

    X = randn(rng, T, n, p)

    # Assign metabolites to subclasses
    counts = floor.(Int, m .* prop_sub)
    counts[end] += m - sum(counts)

    subclass_of_met = Int[]
    for h in 1:H
        append!(subclass_of_met, fill(h, counts[h]))
    end

    subclass_of_met = subclass_of_met[1:m]

    # Global effect per covariate
    theta0 = randn(rng, T, p) .* T(theta0_scale)

    # Subclass effects beta[k,h]
    beta = Matrix{T}(undef, p, H)
    for k in 1:p
        for h in 1:H
            beta[k, h] = theta0[k] + T(tau_v) * randn(rng, T)
        end
    end

    # Metabolite-level true coefficients B[k,j]
    B_true = Matrix{T}(undef, p, m)
    for k in 1:p
        for j in 1:m
            h = subclass_of_met[j]
            B_true[k, j] = beta[k, h] + T(tau_w) * randn(rng, T)
        end
    end

    E = T(sigma_y) .* randn(rng, T, n, m)

    Y = X * B_true + E

    Z = Matrix{T}(I, m, m)

    coef_names = ["x$(k)" for k in 1:p]
    coef_lookup = make_sim_coef_lookup(coef_names)

    return (
        X = Matrix{Float64}(X),
        Y = Matrix{Float64}(Y),
        Z = Matrix{Float64}(Z),
        B_true = Matrix{Float64}(B_true),
        theta0 = Vector{Float64}(theta0),
        beta = Matrix{Float64}(beta),
        subclass_of_met = Int.(subclass_of_met),
        H = H,
        coef_names = coef_names,
        pseudo_coef_names = coef_names,
        coef_lookup = coef_lookup,
        n = n,
        m = m,
        p = p,
        sigma_y = Float64(sigma_y),
        tau_v = Float64(tau_v),
        tau_w = Float64(tau_w)
    )
end