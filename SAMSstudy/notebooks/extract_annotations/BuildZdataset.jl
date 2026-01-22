# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# # Buidling data set to generate Z matrix, based on Goslin
# ---

# This notebook generates data sets to build our Z matrices based on our model:       
#     
#         𝐘 = 𝐗 𝛃 𝐙' + 𝐄 , where    
#
# - 𝐘 is the outcome data matrix (trigliceride levels)
# - 𝐗 is a design matrix constructed from individual covariates (phenotype status, whether they were taking fish oil)
# - 𝛃 is the coefficient matrices to estimate
# - 𝐙 is a design matrix constructed from outcome covariates (e.g., information about the triglicerides)
# - 𝐄 is random error assumed to have mean zero, independent across individuals, but possibly correlated across outcomes
#
# Goslin is a Lipid Shorthand Nomenclature Grammars and uses ANTLRv4 compatible context-free EBNF grammars; we use the R package `rgoslin` to extract lipids information such as oxidation, number of double bonds... 

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
        ENV["R_HOME"] = "C:/PROGRA~1/R/R-44~1.0" # from R.home() in R
    else 
        ENV["R_HOME"] = "/usr/lib/R"

    end
    Pkg.build("RCall")
end         

using CSV, DataFrames, DataFramesMeta, Missings, CategoricalArrays
using StatsBase
# using StatsBase, Statistics, MultivariateStats, MatrixLM, Random, Distributions
using RCall
# using LinearAlgebra, RCall
# using FreqTables, Plots
# using StatsPlots

# ### Ext. Functions

include(joinpath(@__DIR__,"..", "..","src","utils.jl" ));

# ### Load data

# Load look up table for the positive lipids.

# Load look up table `dfXref` contains the negative and the positve lipid ID (e.g., posLip1762) 
# and the original names of the potential identified lipids molecule (e.g. 4_TG(22:4_22:4_22:4)+NH4)
lipidsXref = realpath((@__DIR__)*"/../../data/data_processed/inl2b_Lipids_Xref.csv")
dfXref = DataFrame(CSV.File(lipidsXref));

# ## Extract first candidate of each lipids ID

# Create a column containing the first cleaned candidate.

dfXref.Lipids = string.(map(x -> getLipName(x), dfXref.OriginalNames));

# Let create extra columns that include information about oxydation ("Ox") and plasmogen information ("plasmanyl" or "plasmenyl"):

# Create oxydation and plasmogen column
dfXref.Oxidation = repeat(["no"], size(dfXref)[1]);
dfXref.Plasmalogen = repeat(["no"], size(dfXref)[1]);

# Example:

# dfXref[sample(1:size(dfXref)[1], 3, replace= false),:] # 3 random lipids
dfXref[[495, 496, 696],:]

# We need to convert the nomenclature to be compatible with `GOSLIN` package nomenclature and fill in `Oxidation` and `Plasmalogen` columns:

# +
# Clean Names 
for i in 1:size(dfXref)[1]
    dfXref.Lipids[i] = replace(dfXref.Lipids[i], "-H"=>"")
    dfXref.Lipids[i] = replace(dfXref.Lipids[i], "SPLASH_"=>"")
    dfXref.Lipids[i] = replace(dfXref.Lipids[i], "SPLASH/"=>"")
    if occursin("AcCar", dfXref.Lipids[i])
        dfXref.Lipids[i] = replace(dfXref.Lipids[i], "AcCar"=>"CAR" )
    else
        if occursin(r"(?i)plasmenyl-", dfXref.Lipids[i])
            dfXref.Lipids[i] = replace(dfXref.Lipids[i], r"(?i)plasmenyl-"=>"" )
            dfXref.Plasmalogen[i] = "plasmenyl"
        else
            if occursin(r"(?i)plasmanyl-", dfXref.Lipids[i])
            dfXref.Lipids[i] = replace(dfXref.Lipids[i], r"(?i)plasmanyl-"=>"" )
            dfXref.Plasmalogen[i] = "plasmanyl"
            else
                if occursin("Cer-NS", dfXref.Lipids[i])
                    dfXref.Lipids[i] = replace(dfXref.Lipids[i], "Cer-NS"=>"Cer" )
                else
                    if occursin("Ox", dfXref.Lipids[i])
                        dfXref.Lipids[i] = replace(dfXref.Lipids[i], "Ox"=>"")
                        dfXref.Oxidation[i] = "yes"
                    else
                        if occursin("AcCa(", dfXref.Lipids[i])
                                dfXref.Lipids[i] = replace(dfXref.Lipids[i], "AcCa("=>"CAR(" )
                        end
                    end
                end
            end
        end
    end
    
    
end
# -

replace("TG(18:8_20:2_20:5)", "_" => "/")

# ### Save intermediate results

# +
# realpath((@__DIR__)*"/../../data/dataprocessed/")

dfXref |> CSV.write("../../data/data_processed/dfXrefTest.csv")
# -

# ## Build Z data sets

# To build our raw Z dataset:
# * we check the validity of the nomenclature grammar and filter invalid ones
# * we filter the class of lipids of interest, for example we can keep only Triglycerides by filtering lipids containing the string "TG". 

# ### Build Z data set containing only Triglycerides

R"""
suppressMessages(library(rgoslin))
suppressMessages(library(tidyverse));
"""

@rput dfXref;

# #### Lipids goup  preselection: triglycerides

R"""
# check validity
dfXref$Valid <- suppressWarnings(sapply(dfXref$Lipids, isValidLipidName))

# filter for valid and TG lipids
dfXref <- dfXref %>% filter(Valid) %>% filter(str_detect(Lipids, "TG"))

# get lipids information: total C and Double Bonds 
dfGoslin <- as_tibble(parseLipidNames(dfXref$Lipids))[, c("Original.Name", "Total.C", "Total.DB", "Lipid.Maps.Main.Class")];

# change columns name
names(dfGoslin) <- c("Lipids", "Total_C", "Total_DB", "Class")

# dfGoslin <- right_join(subset(dfXref, select = -c(Valid)), dfGoslin) # remove Valid col
dfGoslin <- cbind(subset(dfXref, select = -c(Valid)), subset(dfGoslin, select = c("Total_C", "Total_DB", "Class")))

"""
@rget dfGoslin;

# We added two additional information: total number of carbon and total number of double bound.

# +
# Parse Total_C and Total_DB
# dfGoslin.Total_C =parse.(Int64, dfGoslin.Total_C);
# dfGoslin.Total_DB =parse.(Int64, dfGoslin.Total_DB);
# -

# #### Save raw triglycerides Z matrix

dfGoslin |> CSV.write("../../data/data_processed/ZmatRawTG.csv")

# ### Build Z matrix-PE, PC, PA

R"""
suppressMessages(library(rgoslin))
suppressMessages(library(tidyverse));
"""

@rput dfXref;

# #### Lipids goup  preselection: phospholipids

R"""
# check validity
dfXref$Valid <- suppressWarnings(sapply(dfXref$Lipids, isValidLipidName))

# filter for valid and Phospholipids
dfXref <- dfXref %>% filter(Valid) %>% filter(str_detect(Lipids, "PE|PC|PA")) %>% 
            filter(!str_detect(Lipids, "TG")) %>% # remove Plasmalogen TG 
            filter(!str_detect(Lipids, "L")) # remove Lysophospholipids

# get lipids information: total C and Double Bonds 
dfGoslin <- as_tibble(parseLipidNames(dfXref$Lipids))[, c("Original Name", "Total C", "Total DB", "Lipid Maps Main Class")];

# change columns name
names(dfGoslin) <- c("Lipids", "Total_C", "Total_DB", "Class")

dfGoslin <- right_join(subset(dfXref, select = -c(Valid)), dfGoslin) # remove Valid col

"""
@rget dfGoslin;

# We added two additional information: total number of carbon and total number of double bound.

# Parse Total_C and Total_DB
dfGoslin.Total_C =parse.(Int64, dfGoslin.Total_C);
dfGoslin.Total_DB =parse.(Int64, dfGoslin.Total_DB);

# #### Save raw phospholipids Z matrix

dfGoslin |> CSV.write("../../data/data_processed/ZmatRawPhos.csv")

# ### Build Z matrix- all lipids

R"""
suppressMessages(library(rgoslin))
suppressMessages(library(tidyverse));
"""

@rput dfXref;

# #### Lipids goup  preselection: all lipids

# + jupyter={"outputs_hidden": true}
R"""
# check validity
dfXref$Valid <- sapply(dfXref$Lipids, isValidLipidName)

print(dim(dfXref))

# filter for valid Lipids nomenclature
dfXref <- dfXref %>% filter(Valid)

print(dim(dfXref))
# get lipids information: total C and Double Bonds 
parseLipidNames(dfXref$Lipids[c(1:1)])#[, c("Original Name", "Total C", "Total DB", "Lipid Maps Main Class")];

# change columns name
# names(dfGoslin) <- c("Lipids", "Total_C", "Total_DB", "Class")
#print(dim(dfXref))


"""
# @rget dfGoslin;
# -

names(dfXref)

@rget dfXref

sum(dfXref.Valid)

R"""
# check validity
dfXref$Valid <- suppressWarnings(sapply(dfXref$Lipids, isValidLipidName))
# filter for valid Lipids nomenclature
dfXref <- dfXref %>% filter(Valid)
# get lipids information: total; C and Double Bonds 
dfGoslin <- as_tibble(parseLipidNames(dfXref$Lipids))[, c("Original Name", "Total C", "Total DB", "Lipid Maps Main Class")];
# change columns name
names(dfGoslin) <- c("Lipids", "Total_C", "Total_DB", "Class")
dfGoslin <- right_join(subset(dfXref, select = -c(Valid)), dfGoslin) # remove Valid col

"""
@rget dfGoslin;

# We added two additional information: total number of carbon and total number of double bound.

# Parse Total_C and Total_DB
dfGoslin.Total_C =parse.(Int64, dfGoslin.Total_C);
dfGoslin.Total_DB =parse.(Int64, dfGoslin.Total_DB);

# #### Save raw Z matrix for all lipids

dfGoslin |> CSV.write("../../data/data_processed/ZmatRawAll.csv")

countmap(dfGoslin.Class)

versioninfo()

R"""
sessionInfo()
"""


