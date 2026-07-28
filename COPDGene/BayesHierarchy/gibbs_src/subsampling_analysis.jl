############################################################
# subsampling_analysis.jl
#
# Purpose:
#   Repeated COPDGene subsampling analysis and calculation
#   of reference MSE metrics for all covariates.
############################################################

using DataFrames
using Random
using Statistics

"""
Create repeated complementary training/reference splits for
each requested training size.

Only row indices are stored to avoid retaining many copies of
the full X and Y matrices in memory.
"""
function make_subsample_splits(
    Y::AbstractMatrix,
    X::AbstractMatrix;
    train_sizes = [
        40, 60, 80, 100, 125, 150,
        200, 250, 300, 400, 500, 600
    ],
    R::Int = 100,
    seed::Int = 2026
)
    n = size(Y, 1)

    @assert size(X, 1) == n "X and Y must have the same number of subjects"
    @assert R > 0 "R must be positive"
    @assert all(
        (0 .< train_sizes) .& (train_sizes .< n)
    ) "Training sizes must lie between 1 and n - 1"

    rng = MersenneTwister(seed)
    splits = NamedTuple[]

    for n_train in train_sizes
        for repetition in 1:R
            permutation = randperm(rng, n)

            push!(
                splits,
                (
                    n_train = n_train,
                    n_test = n - n_train,
                    repetition = repetition,
                    train_idx = permutation[1:n_train],
                    test_idx = permutation[(n_train + 1):end]
                )
            )
        end
    end

    return splits
end

"""
Fit MatrixLM on each training and reference subset and fit the
Bayesian hierarchy using the training estimates.
"""
function fit_repeated_subsamples(
    splits,
    Y::AbstractMatrix,
    X::AbstractMatrix,
    Z::AbstractMatrix,
    subclass_of_met::Vector{Int},
    H::Int;
    mu0::Float64 = 0.0,
    s0::Float64 = 1.0,
    halfcauchy_scale::Float64 = 1.0,
    n_iter::Int = 5000,
    burnin::Int = 1000,
    thin::Int = 1,
    seed0::Int = 1000,
    keep_results::Bool = false
)
    fit_results = NamedTuple[]
    training_sizes = unique(
        [split.n_train for split in splits]
    )

    repetitions_by_size = Dict(
        n_train => maximum(
            split.repetition for split in splits
            if split.n_train == n_train
        )
        for n_train in training_sizes
    )

    for split in splits
        if split.repetition == 1
            println(
                "Starting n_train=$(split.n_train), " *
                "n_reference=$(split.n_test), " *
                "repetitions=$(repetitions_by_size[split.n_train])"
            )
        end

        X_train = X[split.train_idx, :]
        Y_train = Y[split.train_idx, :]
        X_test = X[split.test_idx, :]
        Y_test = Y[split.test_idx, :]

        mlm_train = fit_copdgene_matrixlm(X_train, Y_train, Z)
        mlm_test = fit_copdgene_matrixlm(X_test, Y_test, Z)

        gibbs_seed =
            seed0 +
            10_000 * split.n_train +
            split.repetition

        B_bayes_train, SE_bayes_train, res_by_cov =
            fit_bayes_all_covariates_one(
                mlm_train.coef,
                mlm_train.se,
                subclass_of_met,
                H;
                mu0 = mu0,
                s0 = s0,
                halfcauchy_scale = halfcauchy_scale,
                n_iter = n_iter,
                burnin = burnin,
                thin = thin,
                seed0 = gibbs_seed,
                keep_results = keep_results
            )

        push!(
            fit_results,
            (
                n_train = split.n_train,
                n_test = split.n_test,
                repetition = split.repetition,
                B_mlm_train = Matrix(mlm_train.coef),
                B_mlm_test = Matrix(mlm_test.coef),
                B_bayes_train = B_bayes_train
            )
        )

        n_repetitions = repetitions_by_size[split.n_train]
        if split.repetition == 1 ||
           split.repetition % 10 == 0 ||
           split.repetition == n_repetitions
            println(
                "  Completed repetition $(split.repetition) / $n_repetitions"
            )
        end
    end

    return fit_results
end

"""
Calculate repetition-level reference MSEs and summarize them
by training size and covariate.
"""
function summarize_reference_agreement(
    fit_results;
    covariate_names = nothing
)
    @assert !isempty(fit_results) "fit_results cannot be empty"

    first_B = fit_results[1].B_mlm_train
    n_covariates = size(first_B, 1)

    if covariate_names === nothing
        covariate_names =
            ["Covariate_$k" for k in 1:n_covariates]
    end

    @assert length(covariate_names) == n_covariates

    repetition_table = DataFrame(
        n_train = Int[],
        n_test = Int[],
        repetition = Int[],
        covariate_index = Int[],
        covariate = String[],
        mse_mlm = Float64[],
        mse_bayes = Float64[],
        ratio_mlm_bayes = Float64[],
        percent_reduction = Float64[]
    )

    for result in fit_results
        for k in 1:n_covariates
            mlm_train = vec(result.B_mlm_train[k, :])
            mlm_test = vec(result.B_mlm_test[k, :])
            bayes_train = vec(result.B_bayes_train[k, :])

            keep =
                isfinite.(mlm_train) .&
                isfinite.(mlm_test) .&
                isfinite.(bayes_train)

            @assert any(keep) "No finite paired estimates remain for covariate $k"

            mse_mlm = mean(
                (mlm_test[keep] .- mlm_train[keep]).^2
            )

            mse_bayes = mean(
                (mlm_test[keep] .- bayes_train[keep]).^2
            )

            push!(
                repetition_table,
                (
                    result.n_train,
                    result.n_test,
                    result.repetition,
                    k,
                    string(covariate_names[k]),
                    mse_mlm,
                    mse_bayes,
                    mse_mlm / mse_bayes,
                    100 * (1 - mse_bayes / mse_mlm)
                )
            )
        end
    end

    grouped = groupby(
        repetition_table,
        [:n_train, :n_test, :covariate_index, :covariate]
    )

    summary_table = combine(
        grouped,
        :mse_mlm => mean => :mean_mse_mlm,
        :mse_bayes => mean => :mean_mse_bayes,
        :mse_mlm => std => :sd_mse_mlm,
        :mse_bayes => std => :sd_mse_bayes,
        :ratio_mlm_bayes => mean => :mean_repetition_ratio,
        :ratio_mlm_bayes => median => :median_repetition_ratio,
        :ratio_mlm_bayes => std => :sd_repetition_ratio,
        :ratio_mlm_bayes =>
            (x -> std(x) / sqrt(length(x))) =>
            :se_repetition_ratio,
        [:mse_mlm, :mse_bayes] =>
            ((mse_mlm, mse_bayes) ->
                mean(mse_mlm) / mean(mse_bayes)
            ) =>
            :ratio_of_mean_mses,
        [:mse_mlm, :mse_bayes] =>
            ((mse_mlm, mse_bayes) ->
                100 * (
                    1 -
                    mean(mse_bayes) / mean(mse_mlm)
                )
            ) =>
            :percent_reduction_from_mean_mses,
        nrow => :R
    )

    sort!(
        summary_table,
        [:covariate_index, :n_train]
    )

    return repetition_table, summary_table
end
