#################################################################################
# run_copd_plot.jl
#
# Generate the COPDGene Reference MSE Ratio manuscript figure
# from the saved subsampling summary.
#################################################################################

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames

include(joinpath(@__DIR__, "gibbs_src", "make_plots.jl"))

table_path = joinpath(
    @__DIR__,
    "results",
    "tables",
    "copdgene_rmr_selected_covariates.csv"
)

figure_dir = joinpath(@__DIR__, "results", "figures")
mkpath(figure_dir)

figure_path = joinpath(
    figure_dir,
    "copdgene_rmr.pdf"
)

println("Reading COPDGene Reference MSE Ratio summary...")
summary_table = CSV.read(table_path, DataFrame)

println("Generating COPDGene manuscript figure...")
make_copdgene_rmr_plot(
    summary_table,
    figure_path;
    myfont = "Helvetica",
    mytitlefontsize = 12
)

println("Saved: ", figure_path)
println("Done.")
