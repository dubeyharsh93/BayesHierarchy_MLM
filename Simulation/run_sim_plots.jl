############################################################
# run_sim_plots.jl
#
# Purpose:
#   Generate manuscript simulation figures from the saved
#   three-regime simulation summary.
############################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using LaTeXStrings
using Plots

gr()

table3 = CSV.read(
    joinpath(
        @__DIR__,
        "results",
        "tables",
        "simulation_1level_3hetero_summary.csv"
    ),
    DataFrame
)

mkpath(joinpath(@__DIR__, "results", "figures"))

heterogeneity_labels = Dict(
    "low" =>
        L"\mathrm{Low}\ (\tau_{\mathrm{total}} = 0.072)",
    "moderate" =>
        L"\mathrm{Moderate}\ (\tau_{\mathrm{total}} = 0.144)",
    "high" =>
        L"\mathrm{High}\ (\tau_{\mathrm{total}} = 0.288)"
)

default(
    titlefontsize = 12,
    guidefontsize = 11,
    tickfontsize = 9,
    legendfontsize = 10,
    linewidth = 2.5,
    markersize = 6,
    framestyle = :box,
    grid = false,
    dpi = 300
)

function make_ratio_plot(
    table::DataFrame,
    ratio_column::Symbol;
    ylabel::String,
    output_name::String
)
    p = plot(
        xlabel = "Sample size (n)",
        ylabel = ylabel,
        title = "",
        legend = :topright,
        size = (650, 500),
        margin = 5Plots.mm,
        guidefont = font(11, "Helvetica"),
        tickfont = font(9, "Helvetica")
    )

    for regime in ["low", "moderate", "high"]
        d = subset(
            table,
            :heterogeneity => ByRow(==(regime))
        )

        sort!(d, :n)

        plot!(
            p,
            d.n,
            d[!, ratio_column];
            marker = :circle,
            label = heterogeneity_labels[regime]
        )
    end

    hline!(
        p,
        [1.0];
        color = :black,
        linestyle = :dash,
        linewidth = 1.5,
        label = "Equal performance"
    )

    savefig(
        p,
        joinpath(@__DIR__, "results", "figures", output_name)
    )

    return p
end

p_coef = make_ratio_plot(
    table3,
    :ratio_coef_mse_mean;
    ylabel = "Estimation MSE Ratio (MatrixLM / Bayesian)",
    output_name = "sim_est_ratio.pdf"
)

p_pred = make_ratio_plot(
    table3,
    :ratio_test_mse_mean;
    ylabel = "Prediction MSE Ratio (MatrixLM / Bayesian)",
    output_name = "sim_pred_ratio.pdf"
)

display(p_coef)
display(p_pred)
