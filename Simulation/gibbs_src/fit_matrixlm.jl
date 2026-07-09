############################################################
# fit_matrixlm.jl
#
# Purpose:
#   Fit MatrixLM for COPDGene using X, Y, Z.
#
# Input:
#   X : individual-level predictor matrix
#   Y : metabolite response matrix
#   Z : metabolite feature matrix, identity for one-level model
#
# Output:
#   MatrixLM fit, coefficients, t-statistics, standard errors
############################################################

using MatrixLM
using LinearAlgebra
using Statistics

using MatrixLM

function fit_copdgene_matrixlm(X, Y, Z)
    fit = mlm(
        RawData(Response(Y), Predictors(X, Z)),
        addXIntercept = false,
        addZIntercept = false
    )

    coef = MatrixLM.coef(fit)
    tstat = MatrixLM.t_stat(fit)
    se = abs.(coef ./ tstat)

    return (
        fit = fit,
        coef = coef,
        tstat = tstat,
        se = se
    )
end