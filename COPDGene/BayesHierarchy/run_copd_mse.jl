############################################################################
# run_copd_mse.jl
#
# Reproducible COPDGene Bayesian hierarchical analysis to
# produce the prediction and estimate reproducibility results 
# comparing MatrixLM and Bayesian hierarchical model.
############################################################################

using Pkg
Pkg.activate(@__DIR__)

using Statistics

include(joinpath(@__DIR__,"gibbs_src", "construct_matrices.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_matrixlm.jl"))
include(joinpath(@__DIR__, "gibbs_src", "gibbs_sampler.jl"))
include(joinpath(@__DIR__, "gibbs_src", "train_test_split.jl"))
include(joinpath(@__DIR__, "gibbs_src", "fit_all_covariates.jl"))

function test_mse(Y_test, X_test, B)
    Yhat = X_test * B
    R = Y_test .- Yhat
    return mean(abs2, R)
end

println("Constructing COPDGene matrices...")
obj = construct_copdgene_matrices()

println("Creating train/test split...")
split = train_test_split_rows(obj.Y, obj.X; train_frac = 0.70, seed = 2026)

println("Train size: ", size(split.Y_train))
println("Test size:  ", size(split.Y_test))

println("Fitting MatrixLM on training data...")
mlm_train = fit_copdgene_matrixlm(split.X_train, split.Y_train, obj.Z)

println("Running Bayesian model on training estimates...")
B_bayes_train, SE_bayes_train, res_by_cov = fit_bayes_all_covariates_one(
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

B_mlm = Float64.(mlm_train.coef)
B_bys = Float64.(B_bayes_train)

mse_mlm = test_mse(Float64.(split.Y_test), Float64.(split.X_test), B_mlm)
mse_bys = test_mse(Float64.(split.Y_test), Float64.(split.X_test), B_bys)

println()
println("Prediction results")
println("------------------")
println("MatrixLM test MSE: ", mse_mlm)
println("Bayesian test MSE: ", mse_bys)
println("Prediction ratio MatrixLM / Bayes: ", mse_mlm / mse_bys)

println()
println("Fitting MatrixLM on test data...")
mlm_test = fit_copdgene_matrixlm(split.X_test, split.Y_test, obj.Z)

println("Running Bayesian model on test estimates...")
B_bayes_test, SE_bayes_test, res_by_cov_test = fit_bayes_all_covariates_one(
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

mse_bys_te = mean(abs2, B_bys .- B_bayes_test)
mse_te_tr  = mean(abs2, mlm_test.coef .- mlm_train.coef)

println()
println("Estimate reproducibility results")
println("--------------------------------")
println("MatrixLM estimate MSE: ", mse_te_tr)
println("Bayesian estimate MSE: ", mse_bys_te)
println("Estimate reproducibility ratio MatrixLM / Bayes: ",
        mse_te_tr / mse_bys_te)

println()
println("Done.")