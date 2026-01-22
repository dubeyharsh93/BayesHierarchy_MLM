# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# # Preprocessing Step
# ---

# This notebook carrieds out the preprocessing steps for the metabolomics data:    
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
    ENV["R_HOME"] = "C:/PROGRA~1/R/R-42~1.1" # from R.home() in R
    Pkg.build("RCall")
end     

using CSV, DataFrames
using RCall

using CSV, DataFrames, Missings #, CategoricalArrays
using StatsBase, Statistics #, MultivariateStats#, RCall
using FreqTables #, Plots, StatsPlots

# ### Ext. Functions

include(joinpath(@__DIR__,"..","..","src","preprocessing.jl" ));
include(joinpath(@__DIR__,"..","..","src","wrangle_utils.jl" ));

# ### Load data

# Get reference metabolite file
fileRef = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","refMeta.csv");
dfRef = CSV.read(fileRef, DataFrame);
println("The reference metabolite dataset contains $(size(dfRef, 1)) different metabolites.")

# Get negative metabolite file
dfNegMetabo = readCOPDdata(realpath(joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","negMeta.csv"))) 
println("The negative metabolite dataset contains $(size(dfNegMetabo, 1)) samples and $(size(dfNegMetabo, 2)-1) metabolites.")

# Get polar metabolite file
dfPolarMetabo = readCOPDdata(joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","polarMeta.csv"));
println("The polar metabolite dataset contains $(size(dfPolarMetabo, 1)) samples and $(size(dfPolarMetabo, 2)-1) metabolites.")

# Get positive early metabolite file
dfPosEarlyMetabo = readCOPDdata(joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","posEarlyMeta.csv"));
println("The positive early metabolite dataset contains $(size(dfPosEarlyMetabo, 1)) samples and $(size(dfPosEarlyMetabo, 2)-1) metabolites.")

# Get positive late metabolite file
dfPosLateMetabo = readCOPDdata(joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","posLateMeta.csv"));
println("The positive late metabolite dataset contains $(size(dfPosLateMetabo, 1)) samples and $(size(dfPosLateMetabo, 2)-1) metabolites.")

# ### Join dataframes

# The dataframe `dfPolarMetabo` contains Integer type instead of Float type, which produces an error during the imputation. We need to convert the values to Integer type.

# + tags=[]
typeof(dfPolarMetabo[:, 2])
# -

dfPolarMetabo[!,2:end] .= convert.(Union{Missing, Float64}, dfPolarMetabo[:, 2:end]);
typeof(dfPolarMetabo[:, 2])

df = leftjoin(dfNegMetabo, dfPolarMetabo, on = :SampleID)
leftjoin!(df, dfPosEarlyMetabo, on = :SampleID)
leftjoin!(df, dfPosLateMetabo, on = :SampleID)
size(df)

# ## Imputation

names(dfRef)

df = imputeSPIROMICS(df, dfRef);

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

df[!,2:end] .= convert.(Union{Missing, Float64}, df[:, 2:end]);
df = pqnorm(df, startCol = 2);

# ## Transformation
# ---
#
# A simple and widely used transformation to make data more symmetric and homoscedastic is the log-transformation.

df = log2tx(df, startCol = 2);

first(df)

# + [markdown] kernel="SoS"
# ## Save pretreatments

# + kernel="Julia 1.5.3"
fileMeta = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","inl2_Meta.csv");
df |> CSV.write(fileMeta)
# -

versioninfo()

R"""
sessionInfo()
"""


