############################################################
# run_sim_mse.jl
#
# Purpose:
#   Run simulation study and produce results to compare MSE 
# of Bayesian hierarchical model MatrixLM.
############################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Statistics
using PrettyTables

include(joinpath(@__DIR__, "gibbs_src", "construct_sim_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "train_test_split.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_all_covariates.jl"))

mkpath(joinpath(@__DIR__, "results", "tables"))

function test_mse(Y_te::AbstractMatrix{<:Real},
                  X_te::AbstractMatrix{<:Real},
                  B::AbstractMatrix{<:Real})

    Yhat = X_te * B
    R = Y_te .- Yhat

    return mean(abs2, R)
end

function run_one_simulation_rep_1level(;
    n::Int,
    m::Int,
    p::Int,
    H::Int,
    prop_sub,
    theta0_scale,
    tau_v,
    tau_w,
    sigma_y,
    T,
    split_frac::Float64,
    sim_seed::Int,
    split_seed::Int,
    bayes_seed0::Int,
    bayes_mu0::Float64,
    bayes_s0::Float64,
    bayes_halfcauchy_scale::Float64,
    bayes_n_iter::Int,
    bayes_burnin::Int,
    bayes_thin::Int
)

    sim = simulate_bilinear_identity_1level(
        n = n,
        m = m,
        p = p,
        H = H,
        prop_sub = prop_sub,
        theta0_scale = theta0_scale,
        tau_v = tau_v,
        tau_w = tau_w,
        sigma_y = sigma_y,
        T = T,
        seed = sim_seed
    )

    split = train_test_split_rows(
        sim.Y,
        sim.X;
        train_frac = split_frac,
        seed = split_seed
    )

    mlm_train = fit_copdgene_matrixlm(
        split.X_train,
        split.Y_train,
        sim.Z
    )

    B_bayes_train, SE_bayes_train, res_train =
        fit_bayes_all_covariates_one(
            mlm_train.coef,
            mlm_train.se,
            sim.subclass_of_met,
            sim.H;
            mu0 = bayes_mu0,
            s0 = bayes_s0,
            halfcauchy_scale = bayes_halfcauchy_scale,
            n_iter = bayes_n_iter,
            burnin = bayes_burnin,
            thin = bayes_thin,
            seed0 = bayes_seed0,
            keep_results = false
        )

    B_mlm = Float64.(mlm_train.coef)
    B_bys = Float64.(B_bayes_train)
    B_true = Float64.(sim.B_true)

    coef_mse_mlm = mean(abs2, B_mlm .- B_true)
    coef_mse_bys = mean(abs2, B_bys .- B_true)

    test_mse_mlm = test_mse(split.Y_test, split.X_test, B_mlm)
    test_mse_bys = test_mse(split.Y_test, split.X_test, B_bys)

    return (
        coef_mse_mlm = coef_mse_mlm,
        coef_mse_bayes = coef_mse_bys,
        ratio_coef_mse = coef_mse_mlm / coef_mse_bys,
        test_mse_mlm = test_mse_mlm,
        test_mse_bayes = test_mse_bys,
        ratio_test_mse = test_mse_mlm / test_mse_bys
    )
end

function run_simstudy_grid_1level_2hetero(;
    nrep::Int = 10,

    nm_grid = [
        (n = 60,  m = 300),
        (n = 60,  m = 770),
        (n = 200, m = 300),
        (n = 200, m = 770),
        (n = 400, m = 770),
    ],

    sim_base_kwargs = (
        p = 10,
        H = 4,
        prop_sub = [0.55, 0.25, 0.15, 0.05],
        theta0_scale = 0.2f0,
        sigma_y = 1.0f0,
        T = Float32
    ),

    hetero_regimes = [
        (name = "moderate", tau_v = 0.12f0, tau_w = 0.08f0),
        (name = "large",    tau_v = 0.30f0, tau_w = 0.12f0),
    ],

    split_frac::Float64 = 0.70,
    seed_base::Int = 1,

    bayes_s0::Float64 = 1.0,
    bayes_halfcauchy_scale::Float64 = 1.0,
    bayes_n_iter::Int = 2000,
    bayes_burnin::Int = 100,
    bayes_thin::Int = 1,
    bayes_mu0::Float64 = 0.0
)

    rows = NamedTuple[]

    for (rid,regime) in enumerate(hetero_regimes)
        for (sid, nm) in enumerate(nm_grid)
            for rep in 1:nrep

                println("Running regime=$(regime.name), n=$(nm.n), m=$(nm.m), rep=$rep / $nrep")

                sim_seed = seed_base + 10_000_000 * rid + 1_000_000 * sid + 10_000 * rep
                split_seed = sim_seed + 777
                bayes_seed0 = sim_seed + 2000

                #sim_seed = seed_base + 10_000 * rep + 100 * nm.n + nm.m
                #split_seed = seed_base + rep
                #bayes_seed0 = 1000 + rep

                res = run_one_simulation_rep_1level(
                    n = nm.n,
                    m = nm.m,
                    p = sim_base_kwargs.p,
                    H = sim_base_kwargs.H,
                    prop_sub = sim_base_kwargs.prop_sub,
                    theta0_scale = sim_base_kwargs.theta0_scale,
                    tau_v = regime.tau_v,
                    tau_w = regime.tau_w,
                    sigma_y = sim_base_kwargs.sigma_y,
                    T = sim_base_kwargs.T,
                    split_frac = split_frac,
                    sim_seed = sim_seed,
                    split_seed = split_seed,
                    bayes_seed0 = bayes_seed0,
                    bayes_mu0 = bayes_mu0,
                    bayes_s0 = bayes_s0,
                    bayes_halfcauchy_scale = bayes_halfcauchy_scale,
                    bayes_n_iter = bayes_n_iter,
                    bayes_burnin = bayes_burnin,
                    bayes_thin = bayes_thin
                )

                push!(
                    rows,
                    (
                        heterogeneity = regime.name,
                        n = Float64(nm.n),
                        m = Float64(nm.m),
                        rep = rep,
                        coef_mse_mlm = res.coef_mse_mlm,
                        coef_mse_bayes = res.coef_mse_bayes,
                        ratio_coef_mse = res.ratio_coef_mse,
                        test_mse_mlm = res.test_mse_mlm,
                        test_mse_bayes = res.test_mse_bayes,
                        ratio_test_mse = res.ratio_test_mse
                    )
                )
            end
        end
    end

    df = DataFrame(rows)

    df_summary = combine(
        groupby(df, [:heterogeneity, :n, :m]),
        :ratio_coef_mse => mean => :ratio_coef_mse_mean,
        :ratio_test_mse => mean => :ratio_test_mse_mean
    )

    hetero_order = Dict("moderate" => 1, "large" => 2)
    sort!(
        df_summary,
        [:heterogeneity, :n, :m],
        by = x -> x isa String ? hetero_order[x] : x
    )

    return df, df_summary
end

println("Running simulation study...")

df_all, df_summary = run_simstudy_grid_1level_2hetero()

CSV.write(
    joinpath(@__DIR__, "results", "tables", "simulation_1level_2hetero_all_reps.csv"),
    df_all
)

CSV.write(
    joinpath(@__DIR__, "results", "tables", "simulation_1level_2hetero_summary.csv"),
    df_summary
)

df_rounded_sel = select(
    df_summary,
    :heterogeneity => :Heterogeneity,
    :n,
    :m,
    :ratio_coef_mse_mean,
    :ratio_test_mse_mean
)

df_rounded_sel.ratio_coef_mse_mean =
    round.(df_rounded_sel.ratio_coef_mse_mean; digits = 3)

df_rounded_sel.ratio_test_mse_mean =
    round.(df_rounded_sel.ratio_test_mse_mean; digits = 3)

println()
println("Simulation summary:")
println(df_rounded_sel)

println()
println("LaTeX table:")
pretty_table(
    df_rounded_sel;
    backend = :latex,
    table_format = LatexTableFormat(
        borders = LatexTableBorders(
            top_line = "\\toprule",
            header_line = "\\midrule",
            merged_header_cell_line = "",
            middle_line = "",
            bottom_line = "\\bottomrule"
        )
    ),
    alignment = :c
)

println()
println("Saved:")
println(joinpath(@__DIR__, "results", "tables", "simulation_1level_2hetero_all_reps.csv"))
println(joinpath(@__DIR__, "results", "tables", "simulation_1level_2hetero_summary.csv"))
println("Done.")