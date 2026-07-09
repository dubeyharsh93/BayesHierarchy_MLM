####################################################################
# run_pans_mse.jl
#
# Purpose:
#   Run the Bayesian hierarchical model on the PANSTEATITIS data and 
#   compute MSE ratios for predictions and estimate reproducibility.
####################################################################

using Pkg
Pkg.activate(@__DIR__)

using Statistics

include(joinpath(@__DIR__, "gibbs_src", "construct_pans_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "train_test_split.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_all_covariates.jl"))

function test_mse(Y_te::AbstractMatrix{<:Real},
                  X_te::AbstractMatrix{<:Real},
                  B::AbstractMatrix{<:Real})

    n_te, m = size(Y_te)
    @assert size(X_te, 1) == n_te
    p = size(X_te, 2)
    @assert size(B) == (p, m)

    Yhat = X_te * B
    R = Y_te .- Yhat

    mse_all = mean(abs2, R)
    mse_met = vec(mean(abs2, R; dims = 1))
    mse_ind = vec(mean(abs2, R; dims = 2))

    return mse_all, mse_met, mse_ind
end

println("Constructing PANSTEATITIS matrices...")
obj = construct_pans_matrices()

println("Creating train/test split...")
split = train_test_split_rows(
    obj.Y,
    obj.X;
    train_frac = 0.70,
    seed = 5
)

println("Train size: ", size(split.Y_train))
println("Test size:  ", size(split.Y_test))

println("Fitting MatrixLM on training data...")
mlm_train = fit_copdgene_matrixlm(
    split.X_train,
    split.Y_train,
    obj.Z
)

println("Running one-level Bayesian model on training estimates...")
B_bayes_train, SE_bayes_train, res_train =
    fit_bayes_all_covariates_one(
        mlm_train.coef,
        mlm_train.se,
        obj.subclass_of_met,
        obj.H;
        mu0 = 0.0,
        s0 = 1.0,
        n_iter = 5000,
        burnin = 1000,
        thin = 1,
        seed0 = 1000,
        keep_results = false
    )

Y_test_f = Float64.(split.Y_test)
X_test_f = Float64.(split.X_test)

B_mlm = Float64.(mlm_train.coef)
B_bys = Float64.(B_bayes_train)

mse_mlm, mse_mlm_met, mse_mlm_ind =
    test_mse(Y_test_f, X_test_f, B_mlm)

mse_bys, mse_bys_met, mse_bys_ind =
    test_mse(Y_test_f, X_test_f, B_bys)

println()
println("Prediction results")
println("------------------")
println("Test MSE (MatrixLM): ", mse_mlm)
println("Test MSE (Bayes):    ", mse_bys)
println("Ratio of Test MSE MLM and Bayes: ", mse_mlm / mse_bys)

println()
println("Fitting MatrixLM on test data...")
mlm_test = fit_copdgene_matrixlm(
    split.X_test,
    split.Y_test,
    obj.Z
)

println("Running one-level Bayesian model on test estimates...")
B_bayes_test, SE_bayes_test, res_test =
    fit_bayes_all_covariates_one(
        mlm_test.coef,
        mlm_test.se,
        obj.subclass_of_met,
        obj.H;
        mu0 = 0.0,
        s0 = 1.0,
        n_iter = 5000,
        burnin = 1000,
        thin = 1,
        seed0 = 1000,
        keep_results = false
    )

mse_bys_te = mean(abs2, B_bayes_train .- B_bayes_test)
mse_mlm_te = mean(abs2, mlm_train.coef .- mlm_test.coef)

println()
println("Estimate reproducibility results")
println("--------------------------------")
println("MSE between Bayes estimates (train vs test): ", mse_bys_te)
println("MSE between MatrixLM estimates (train vs test): ", mse_mlm_te)
println("Ratio of MSEs (MatrixLM / Bayes): ", mse_mlm_te / mse_bys_te)

println()
println("Done.")