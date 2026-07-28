############################################################
# make_plots.jl
#
# Purpose:
#   Generate the COPDGene Reference MSE Ratio figure.
############################################################

using DataFrames
using Plots

function make_copdgene_rmr_plot(
    summary_table::DataFrame,
    output_path::String;
    myfont::String = "Helvetica",
    mytitlefontsize::Int = 12
)
    covariates_of_interest = ["Age", "BMI", "COPD: 1"]

    display_names = Dict(
        "Age" => "Age",
        "BMI" => "BMI",
        "COPD: 1" => "COPD Status"
    )

    default(
        titlefontsize = mytitlefontsize,
        guidefontsize = 11,
        tickfontsize = 9,
        legendfontsize = 10,
        linewidth = 2.5,
        markersize = 6,
        framestyle = :box,
        grid = false,
        dpi = 300
    )

    p = plot(
        xlabel = "Training sample size",
        ylabel = "Reference MSE Ratio\n(MatrixLM / Bayesian)",
        title = "",
        legend = :topright,
        size = (650, 500),
        margin = 5Plots.mm,
        guidefont = font(11, myfont),
        tickfont = font(9, myfont)
    )

    for covariate_name in covariates_of_interest
        d = subset(
            summary_table,
            :covariate => ByRow(==(covariate_name))
        )

        @assert !isempty(d) "No results found for $covariate_name"
        sort!(d, :n_train)

        plot!(
            p,
            d.n_train,
            d.ratio_of_mean_mses;
            marker = :circle,
            label = display_names[covariate_name]
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

    savefig(p, output_path)

    return p
end
