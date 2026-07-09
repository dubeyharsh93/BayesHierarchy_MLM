############################################################
# make_plots.jl
#
# Purpose:
#   Make SAMS Bayesian vs MatrixLM comparison plots.
############################################################

using Plots
using DataFrames
using Statistics

function make_comparison_plot(cmp_df::DataFrame, output_path::String)

    gr()  

    p1 = scatter(
        cmp_df.matrixlm_estimate, cmp_df.bayes_estimate,
        xlabel = "MatrixLM Estimate",
        ylabel = "Bayesian Estimate",
        title  = "",
        legend = false
    )
    plot!(
        p1,
        [minimum(cmp_df.matrixlm_estimate), maximum(cmp_df.matrixlm_estimate)],
        [minimum(cmp_df.matrixlm_estimate), maximum(cmp_df.matrixlm_estimate)],
        l = :dash, c = :black
    )

    p2 = scatter(
        log.(cmp_df.matrixlm_se), log.(cmp_df.bayes_se),
        xlabel = "MatrixLM std error (log)",
        ylabel = "Bayesian std error (log)",
        title  = "",
        legend = false
    )
    plot!(
        p2,
        [minimum(log.(cmp_df.matrixlm_se)), maximum(log.(cmp_df.matrixlm_se))],
        [minimum(log.(cmp_df.matrixlm_se)), maximum(log.(cmp_df.matrixlm_se))],
        l = :dash, c = :black
    )

    # Combine side by side
    p = plot(p1, p2, layout = (1, 2), size = (700, 400))

    savefig(p, output_path)

    return p
end