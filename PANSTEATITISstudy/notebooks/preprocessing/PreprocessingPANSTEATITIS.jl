# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# # Preprocessing Step
# ---

# This notebook carried out the preprocessing steps for the metabolomics data:    
# - Imputation
# - Normalization
# - Log2 Transformation

# ## Input

# ### Libraries

# To use RCall for the first time, one needs to 
# the location of the R home directory.
firstTimeRCall = false
if firstTimeRCall
    using Pkg
    io = IOBuffer()
    versioninfo(io)
    if occursin("Windows", String(take!(io)))
        ENV["R_HOME"] = "C:/PROGRA~1/R/R-43~1.1" # from R.home() in R
    else 
        ENV["R_HOME"] = "/usr/lib/R"

    end
    Pkg.build("RCall")
end      

using CSV, DataFrames, Missings#, CategoricalArrays
using StatsBase, Statistics#, MultivariateStats
using FreqTables#, Plots, StatsPlots
using RCall 

# ### Ext. Functions

include(joinpath(@__DIR__,"..","..","src","preprocessing.jl" ));
include(joinpath(@__DIR__,"..","..","src","wrangling_utils.jl" ));

# ### Load data

# #### Reference file

# Get reference metabolite file
fileRef = joinpath(@__DIR__,"..","..","data","processed","refMeta.csv");
dfRef = CSV.read(fileRef, DataFrame);
print_df_size(dfRef)

# #### Metabolite signatures

# Get negative metabolite file
fileMetabo = realpath(joinpath(@__DIR__,"..","..","data","processed","Metabo.csv"));
dfMetabo = CSV.read(fileMetabo, DataFrame);
println("The negative metabolite dataset contains $(size(dfMetabo, 1)) samples and $(size(dfMetabo, 2)-1) metabolites.")

# ## Imputation

# Check if imputation is needed:

summary_variables_missing(dfMetabo)

# ## Normalization
# ----

# ### Probabilistic Quotient Normalization
#
# > 1. Perform an integral normalization (typically a constant
# integral of 100 is used).
# > 2. Choose/calculate the reference spectrum (the best approach
# is the calculation of the median spectrum of control samples).
# > 3. Calculate the quotients of all variables of interest of the test
# spectrum with those of the reference spectrum.
# > 4. Calculate the median of these quotients.
# > 5. Divide all variables of the test spectrum by this median.
#

df = pqnorm(dfMetabo, startCol = 2);
first(df, 3)

# ## Transformation
# ---
#
# A simple and widely used transformation to make data more symmetric and homoscedastic is the log-transformation.

df = log2tx(df, startCol = 2);
first(df, 2)

# + [markdown] kernel="SoS"
# ## Save pretreatments

# + kernel="Julia 1.5.3"
fileMeta = joinpath(@__DIR__,"..","..","data","processed","nl2_Meta.csv");
df |> CSV.write(fileMeta)
# -

versioninfo()

R"""
sessionInfo()
"""


