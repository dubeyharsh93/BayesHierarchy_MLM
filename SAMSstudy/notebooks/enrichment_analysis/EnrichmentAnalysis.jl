# -*- coding: utf-8 -*-
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

# # Enrichment Analysis using R
# ____
# ### Gregory Farage, Śaunak Sen *et al.*
#     gfarage@uthsc.edu / sen@uthsc.edu
#     Division of Biostatistics
#     Department of Preventive Medicine
#     University of Tennessee Health Science Center
#     Memphis, TN
#
# This notebook carries out encrichment analysis.

# ## Background
#
# For decades, statins have been effective and widely popular cholesterol-lowering agents with substantial benefits for preventing and treating cardiovascular disease.
# However, some patients can not tolerate statins. Statin intolerance is usually associated with muscle pain as side effects, also known as statin-associated muscle symptoms (SAMS). 
#
# **SAMS are particularly difficult to treat:**
# * no validated biomarkers or tests to confirm patient self-reports of SAMS
# * some self-reported patients have non-specific muscle pain not attributable to statin therapy
#
#
# ## Project Summary
#
# This project seeks to identify biomarkers and pathway's related to susceptibility to SAMS based on metabolomic and lipidomic studies.   
#
# **Comparisons of metabolomic/lipidomic profiling are made for:**
# * CS-baseline group vs CS group (paired samples)
# * CN group vs CS group (unpaired samples)
#
# **Group’s acronym definition:**
# * CS-baseline group: Patients with documented SAMS off statin    
#
# * CS group: Patients with documented SAMS following rechallenge with statin (Rechallenge) for up to 4 weeks or until symptomatic.     
#        
# * CN group: Controls patients who have no history of SAMS and who are currently treated with statin.

# ## Libraries

using Pkg 

Pkg.activate("../..")

using CSV, DataFrames, DataFramesMeta, Missings, CategoricalArrays
using StatsBase, Statistics, MatrixLM
using Random, Distributions, StatsModels#, MultivariateStats
using LinearAlgebra, Latexify
using FreqTables, Plots, StatsPlots, Images, FileIO
using Plots.PlotMeasures

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

# ## External functions 

include(joinpath(@__DIR__, "..","..","src","mLinearModel.jl" ));
include(joinpath(@__DIR__,"..","..", "src","myPlots.jl" ));

# ## Load data

# +
# Load look up table for lipids
LipidsXref = realpath((@__DIR__)*"/../../data/data_processed/inl2b_Lipids_Xref.csv")
dfLipidsXref = DataFrame(CSV.File(LipidsXref));

# Load data set
fileLipids = realpath((@__DIR__)*"/../../data/data_processed/inl2b_Lipids.csv")
dfLipids = DataFrame(CSV.File(fileLipids));
# -

# ### Filter cases

# Select only CN (control - no SAMS history) and CS (cases) patients on statin.

# +
# true => unpaired
# false => paired (rechallenge)
isunpairFlag = true;
dfInd, dfMeta = getCases_new(dfLipids, isunpaired = isunpairFlag);
slctFishOil = true

# Check to filter none fish oil users 
hasFishOil = true;
# -

# ### Standardize

# standardize lipids mass-to-charge ratio
funStandardize!(dfMeta, isunpaired = isunpairFlag);

# ## Individual characteristics
#
# Two-way contingency table considering `Group` and `FishOil` variables: 
#

fullTable = freqtable(dfInd[:, [:FishOil, :Group]], :FishOil, :Group)
matFullfreq = collect(fullTable[:,:]);
matFullfreq = vcat(matFullfreq, sum(matFullfreq, dims= 1))
matFullfreq = hcat(matFullfreq, sum(matFullfreq, dims= 2))
dfFreqFull = DataFrame(matFullfreq, :auto)
rename!(dfFreqFull, Symbol.(vcat(names(fullTable,2), "Total")))
insertcols!(dfFreqFull, 1, :_ => vcat("Fish Oil: ".*(names(fullTable,1)), Symbol("Total")))

# ## Lipid characteristics
#
# Each lipid ID corresponds to the nomenclature of a molecule that has been profiled in our lipidomic: 
#

dfLipidsXref[sample(1:size(dfLipidsXref)[1],3, replace=false, ordered=true), :]

# Frequency table of the present lipids main class according to the Grammar of Succinct Lipid Nomenclature (Goslin): 

# +
# Get Zmat 
dfZrawAll =  DataFrame(CSV.File("../../data/data_processed/ZmatRawAll.csv"))

classLipids = freqtable(dfZrawAll.Class)
DataFrame(Class = names(classLipids, 1), Count = collect(classLipids[:]))

# -

# Since triglycerides (TG) are one of the important constituents of the lipid fraction of the human body, and represent the main lipid component of our dietary fat and fat depots, we are interested to study variation in triglyceride composition.


# ## Using matrix linear models to study variation in triglyceride composition
#
# Matrix linear models provide a framework for studying associations in high-throughput data using bilinear models.
#
# $$Y= X B Z^\prime + E,$$
#
# where
# - ``Y`` is the outcome data matrix (trigliceride levels)
# - ``X`` is a design matrix constructed from individual covariates (phenotype status, whether they were taking fish oil)
# - ``Z`` is a design matrix constructed from outcome covariates (information about the triglicerides)
# - ``E`` is random error assumed to have mean zero, independent across individuals, but possibly correlated across outcomes
#
# $$V(vec(E)) = \Sigma \otimes I.$$
#

load("../../images/matrixlinearmodel.png")

# The model can also be viewed as a bilinear model of $y$ in terms of the $x$ and $z$ covariates:
#
# $$y_{ij} = \sum_{i{=}1}^n \sum_{j{=}1}^m b_{kl} x_{ik} z_{jl} + e_{ij}.$$
#
# The regression coefficients ($b$'s) can be interpreted as interactions between the corresponding $x$ and $z$ covariates.
#

# ## Data analysis

# ### Plotting Attributes

myfont = "Helvetica"
mytitlefontsize = 12; 

# ### X matrix

names(dfInd)

# +
contrasts_sams = Dict(
            :Group => EffectsCoding(base = sort(unique(dfInd.Group))[1]),
            :FishOil => EffectsCoding(base = sort(unique(dfInd.FishOil))[1]),
)

frml_sams = "1 +  FishOil +Group"
formulaX = eval(Meta.parse(string("@formula(0 ~ ", frml_sams, ").rhs")))


if occursin("+", frml_sams) || occursin("*", frml_sams)
    vCovarNames = collect(string.(formulaX))
else
    vCovarNames = [string(formulaX)]
end

if ["1"] ⊆ vCovarNames
    vCovarNames[findall(vCovarNames .== "1")] .= "Intercept"
end

idx2change = findall(occursin.("&",vCovarNames))
vCovarNames[idx2change]  .= replace.(vCovarNames[idx2change], " & " => "Ξ")

mX = modelmatrix(formulaX, dfInd, hints = contrasts_sams);
X = mX

vFrmlNames = Vector{String}()

if frml_sams == "1" #frml_c[] == "1"
    # vFrmlNames = ["(Intercept)"]
    vPseudoFrmlNames = ["(Intercept)"]
else
    sch_copd = schema(formulaX, dfInd, contrasts_sams)
    vFrmlNames = apply_schema(formulaX, sch_copd) |> coefnames
end

function fix_covar_name(s::String)
    s = replace(s, "("=>"", ")"=>"", ": "=>"_")
    s = replace(s, " & " => "Ξ")
end

vPseudoFrmlNames = fix_covar_name.(vFrmlNames)

show(vFrmlNames)


# -

vCovarNames

# ### Y matrix

# Subset triglycerides:

# +
# Get Zmat for triglycerides
rdZ = "Triglycerides"
fZrawTG = "../../data/data_processed/ZmatRawTG.csv"
dfZraw =  DataFrame(CSV.File(fZrawTG));
# need to filter lipids in dfZrawTG because it was created before imputation
lip_ID_TG = intersect(names(dfMeta), dfZraw.lipID);
filter!(:lipID => x -> x in lip_ID_TG, dfZraw)

m =  length(dfZraw.lipID);
# -

dfY = select(dfMeta, dfZraw.lipID);
mY = Matrix(dfY);
Y = mY;

# ### Matrix Linear Model with Z=I 

# +
# Select X design model
# if is4ways is true, covariates are
# CN - No Fish Oil, CN - Fish Oil, CS - No Fish Oil, CS - Fish Oil (CN: Control, CS: Statin Intolerant)
#if is4ways is false, covariates are
# Intercept, Group, Fish Oil, Interaction
is4ways = false;

# ZI is identity matrix 
ZI = convert(Array{Float64, }, collect(I(m)));
CoefZI, CIZI, TstatZI, varZI = getCoefs(mY, mX,ZI);

# CoefZI, CIZI, TstatZI, varZI = getCoefs(df, ZI, responseSelection = dfZraw.lipID,
#                                  isunpaired = isunpairFlag, 
#                                  hasFishOil = hasFishOil, 
#                                  is4ways = is4ways);

# Get X and Y 

# if hasFishOil
#     if is4ways
#         X, Y = getXY4ways(df, responseSelection = dfZraw.lipID,
#                  isunpaired = isunpairFlag);
#         covarSelection = ["1"=> "CN-No Fish Oil (intercept)", "2"=> "CN-Fish Oil",
#                       "3"=> "CS-No Fish Oil", "4"=> "CS-Fish Oil"]
#     else
#         X, Y = getXY(df, responseSelection = dfZraw.lipID,
#                  isunpaired = isunpairFlag);
        covarSelection = ["1"=> "Intercept", "2"=> "SAMS status",
                      "3"=> "Fish Oil", "4"=> "Interaction SAMS-Fish Oil"]
#     end
# else
#     X, Y = getXYnoFishOil(df, responseSelection = dfZraw.lipID,
#                  isunpaired = isunpairFlag);
#     covarSelection = ["1"=> "Intercept", "2"=> "SAMS status"]
# end

# # nameX contains labels for plotting
# if  isunpairFlag
#     if is4ways
#         namesX = ["CN-No Fish Oil" "CN-Fish Oil" "CS-No Fish Oil" "CS-Fish Oil"]
#     else
#         namesX = ["Intercept" "SAMS status" "Fish Oil" "Interaction SAMS-Fish Oil"]
#     end
# else
#     namesX = ["Intercept" "Interaction SAMS-Fish Oil"]
# end;

namesX = vCovarNames
# -

# ### X matrix (individual characteristics)

# We assigned a Group-value of 1 to patients in the CS group and -1 to those in the CN group. Similarly, we use a FishOil-value of 1 to indicate that a patient is taking a fish oil supplement, and -1 to indicate that they are not taking one.

DataFrame(X, vCovarNames) |> (x->first(x, 3))

# ### Y matrix (outcome matrix, lipid concentrations)

dfY = DataFrame(Y, dfZraw.lipID) 
first(dfY, 7)

# ### Z matrix (lipid characteristics)
#
# Triglycerides'usefull information to build our Z matrix:
#

dfZraw[1:3,Not(2)]

# ### Lipid characteristic distributions

histogram2d(dfZraw.Total_DB, dfZraw.Total_C,
            c= :matter, aspectratio = 0.35, xlim = [0,15],
            colorbar_title = "Number of Triglycerides",
            xlabel= "Total Double Bonds", ylabel =  "Total Carbons",
            title = string("Number of ", rdZ, " by the number of carbon \n atoms and the number of unsaturated bonds \n"),
            fontfamily = myfont,
            titlefontsize = mytitlefontsize,
)


# #### Scatterplots of lipid effects analyzed individually

# covarSelection
vCovarNames

# slct1 = covarSelection[3][1];
idxColXmat = 2 #parse(Int64,slct1);

a = 0.25
zlim2 = maximum(abs.(TstatZI[idxColXmat,:]))
p_scatter = scatter(
    dfZraw.Total_DB .+ rand(Uniform(-a,a),length(dfZraw.Total_DB)),
    dfZraw.Total_C,
    zcolor = TstatZI[idxColXmat,:], colorbar = true, m = (:bluesreds, 0.8), 
    colorbar_title = string("\nT-statistics of ",namesX[idxColXmat]),
    xticks=0:1:14, clims = (-zlim2, zlim2), 
    xlabel= "Total Double Bonds", ylabel =  "Total Carbons",
    legend = false, grid = false, tickfontsize=8, 
    fontfamily = myfont,
    right_margin = 10mm,  format = :svg,
    # title = string("T-statistics at the individual level \n for the ", namesX[idxColXmat]),
    titlefontsize = mytitlefontsize,
)

# ## Modeling decisions
#
# We modeled the MLM according to the length of Carbon chain and the degree of unstaturation of the Triglycerides.

# ### Total Carbon 

# Build the Z matrix based on carbon chain length:

# +
dfRefTG = copy(dfZraw)

vbrks_c = [0,40,45,50,55, 60, 65, 70]
lvls_c = vcat(
    string.(vbrks_c[1:end-1]).*" ≤ Total Carbon < ".*string.(vbrks_c[2:end]),
    "Total Carbon ≥ ".*string.(vbrks_c[end]));
vcatc = cut(
            dfRefTG.Total_C, vbrks_c; 
            labels = lvls_c, 
            extend = true
)

dfRefTG.Total_C_cat = vcatc
namesZtc =lvls_c;
# levelsC = lvls_c

# -

lvls_c

# +
# Create a DataFrame from tStats_diff and append column names from df_baseline
dfTstatsZI = DataFrame(
    hcat(permutedims(TstatZI), names(dfY)[1:end]),
    vcat(["Intercept", "FishOil", "Group"], ["lipID"])
);

# Join dfTstatsZI with dfRef to add SuperClassID and SubClassID using CHEM_ID1 as the key
dfTstatsZI = leftjoin(
    dfTstatsZI, 
    dfRefTG[:, [:lipID, :Total_C_cat]], on = :lipID
    # dfRefTG[:, [:lipID, :Total_C_cat , :Total_DB_cat]], on = :lipID
)

# Generate names for covariate figures based on indices
nameCovarFig = "FishOil"

# Group data by SuperClassID
gdf = groupby(dfTstatsZI, :Total_C_cat);

# Calculate the mean T-statistics for each super class and create a new DataFrame
dfMeanTst = DataFrames.combine(gdf, Symbol(nameCovarFig) => mean => Symbol(nameCovarFig)) 
# Sort the DataFrame by SuperClassID
sort!(dfMeanTst, :Total_C_cat);

# Create a dot plot of T-statistics by super class
p_dot = eval(Meta.parse("@df dfTstatsZI dotplot(string.(:Total_C_cat), :$(nameCovarFig), legend = false, markersize = 4)"))
# Overlay a scatter plot on the dot plot with mean values
eval(Meta.parse("@df dfMeanTst scatter!(string.(:Total_C_cat), :$(nameCovarFig), legend = false, color = :orange)"))
# Add a horizontal line at T=0 for reference
hline!([0], color= :red, 
    label = "",
    xlabel = "Total C category", xrotation = 45,
    ylabel = string("T-statistics ", "Treatment"),
    title = "T-statistics per class",
    titlefontsize = mytitlefontsize,
    fontfamily = myfont, grid = false,
)

# Display the plot
plot(p_dot)
# -

# Apply MLM processing:

# +
# Design matrix Z
mZtc_cat = MatrixLM.design_matrix(
	@mlmformula(1 + Total_C_cat),
	dfRefTG,
	 # Dict(:Total_C_cat => StatsModels.FullDummyCoding())
    Dict(:Total_C_cat => StatsModels.DummyCoding(;
            base = lvls_c[1],
            levels = lvls_c
        )
    )
)


# Generate Z matrix
# mZtc_cat = modelmatrix(@formula(y ~ 0 + Total_C_cat).rhs, 
#                     dfRefTG, 
#                     hints = Dict(
#                             :Total_C_cat => StatsModels.FullDummyCoding()));
# mZtc_cat |> x -> hcat(Matrix(dfRefTG[:, [:Total_C_cat]]), x) 

CoefZtc, CIZtc, TstatZtc, varZtc = getCoefs(
    mY, mX,mZtc_cat; 
    hasXIntercept = true, hasZIntercept = true
);
# -

# Plot results:

idx = 1
idxColXmatCI = idxColXmat;
p_tc_ci = confidenceplot(
        CoefZtc[idxColXmatCI, idx:end],
        namesZtc[idx:end],
        CIZtc[idxColXmatCI,idx:end],
        xlabel = namesX[idxColXmatCI]*" Effect Size", legend = false,
        fontfamily = myfont,
        title = string("Confidence interval for \n the ", "carbon count") , 
        titlefontsize = mytitlefontsize,
        size=(400,300),  
)

# ### Double Bonds 

dfZdbcat = DataFrame(intercept = ones(Float64, size(dfZraw)[1]),
        TotalDB_3_6 = ((dfZraw.Total_DB .>= 3) .& (dfZraw.Total_DB .< 6))*1,
        TotalDB_6_9 = ((dfZraw.Total_DB .>= 6) .& (dfZraw.Total_DB .< 9))*1,
        TotalDB_9 = (dfZraw.Total_DB .>= 9)*1	);

# +
vbrks_db = [0,3,6,9];
lvls_db = vcat(
    string.(vbrks_db[1:end-1]).*" ≤ Double Bonds < ".*string.(vbrks_db[2:end]),
    "Double Bonds  ≥ ".*string.(vbrks_db[end]));
    # string.(vbrks_db[end]).*" ≤ Double Bonds");
vcatdb = cut(
            dfRefTG.Total_DB, vbrks_db; 
            labels = lvls_db, 
            extend = true)

dfRefTG.Total_DB_cat = vcatdb

namesZdb = lvls_db;



# +
# Create a DataFrame from tStats_diff and append column names from df_baseline
dfTstatsZI = DataFrame(
    hcat(permutedims(TstatZI), names(dfY)[1:end]),
    vcat(["Intercept", "FishOil", "Group"], ["lipID"])
);

# Join dfTstatsZI with dfRef to add SuperClassID and SubClassID using CHEM_ID1 as the key
dfTstatsZI = leftjoin(
    dfTstatsZI, 
    dfRefTG[:, [:lipID, :Total_DB_cat]], on = :lipID
    # dfRefTG[:, [:lipID, :Total_C_cat , :Total_DB_cat]], on = :lipID
)

# Generate names for covariate figures based on indices
nameCovarFig = "FishOil"

# Group data by SuperClassID
gdf = groupby(dfTstatsZI, :Total_DB_cat);

# Calculate the mean T-statistics for each super class and create a new DataFrame
dfMeanTst = DataFrames.combine(gdf, Symbol(nameCovarFig) => mean => Symbol(nameCovarFig)) 
# Sort the DataFrame by SuperClassID
sort!(dfMeanTst, :Total_DB_cat);

# Create a dot plot of T-statistics by super class
p_dot = eval(Meta.parse("@df dfTstatsZI dotplot(string.(:Total_DB_cat), :$(nameCovarFig), legend = false, markersize = 4)"))
# Overlay a scatter plot on the dot plot with mean values
eval(Meta.parse("@df dfMeanTst scatter!(string.(:Total_DB_cat), :$(nameCovarFig), legend = false, color = :orange)"))
# Add a horizontal line at T=0 for reference
hline!([0], color= :red, 
    label = "",
    xlabel = "Total DB category", xrotation = 45,
    ylabel = string("T-statistics ", "Treatment"),
    title = "T-statistics per class",
    titlefontsize = mytitlefontsize,
    fontfamily = myfont, grid = false,
)

# Display the plot
plot(p_dot)

# +
# Design matrix Z
mZdb_cat = MatrixLM.design_matrix(
	@mlmformula(1 + Total_DB_cat),
	dfRefTG,
	 # Dict(:Total_C_cat => StatsModels.FullDummyCoding())
    Dict(:Total_DB_cat => StatsModels.DummyCoding(;
            base = lvls_db[1],
            levels = lvls_db
        )
    )
)

# mZdb_cat = modelmatrix(@formula(y ~ 0 + Total_DB_cat).rhs, 
#                     dfRefTG, 
#                     hints = Dict(:Total_DB_cat => StatsModels.FullDummyCoding()));

CoefZdb, CIZdb, TstatZdb, varZdb  = getCoefs(
    mY, mX,mZdb_cat; 
    hasXIntercept = true, hasZIntercept = true
);
# -

idx = 1
idxColXmatCI = idxColXmat;
p_db_ci = confidenceplot(
        CoefZdb[idxColXmatCI, idx:end],
        namesZdb[idx:end],
        CIZdb[idxColXmatCI,idx:end],
        xlabel = namesX[idxColXmatCI]*" Effect Size", legend = false,
        fontfamily = myfont,
        title = string("Confidence interval for \n the ", "double bonds") , 
        titlefontsize = mytitlefontsize,
        size=(400,300), 
)

# ### Adjusted 

# +
# Design matrix Z

levelsCdb_cat = vcat(
    ["Intercept"], 
    lvls_c[2:end],
    lvls_db[2:end]
)

mZcdb_cat = MatrixLM.design_matrix(
	@mlmformula(1 + Total_C_cat + Total_DB_cat),
	dfRefTG,
    Dict(
       :Total_DB_cat => StatsModels.DummyCoding(;
            base = lvls_db[1],
            levels = lvls_db
        ),
        :Total_C_cat => StatsModels.DummyCoding(;
            base = lvls_c[1],
            levels = lvls_c
        )
    )
) 

# mZcdb_cat = modelmatrix(
#     @formula(y ~ 1 + Total_C_cat + Total_DB_cat).rhs, 
#     dfRefTG, 
#     hints = Dict(
#         :Total_C_cat => StatsModels.DummyCoding(;
#             base = lvls_c[1],
#             levels = lvls_c
#             ),
#         :Total_DB_cat => StatsModels.DummyCoding(;
#             base = lvls_db[1],
#             levels = lvls_db
#             )
#         )
# );


CoefZcdbcat, CIZcdbcat, TstatZcdbcat, varZcdbcat = getCoefs(
    mY, mX, mZcdb_cat;
    hasXIntercept = true, hasZIntercept = true
);
# -

# Join unadjusted results
namesZtcdb = vcat(namesZtc[2:end], namesZdb[2:end])
CoefZtcdb  = hcat(CoefZtc[:, 2:end] , CoefZdb[:, 2:end])
CIZtcdb    = hcat(CIZtc[:, 2:end] , CIZdb[:, 2:end])
TstatZtcdb = hcat(TstatZtc[:, 2:end] , TstatZdb[:, 2:end])
varZtcdb   = hcat(varZtc[:, 2:end] , varZdb[:, 2:end]);

# +
idx = 1
idxColXmatCI = idxColXmat;
p_ci_unadjusted = confidenceplot(
                    CoefZtcdb[idxColXmatCI, :],	
                    namesZtcdb,
                    CIZtcdb[idxColXmatCI, :],
                    xlabel = namesX[idxColXmatCI]*" Effect Size", legend = false,
                    fontfamily = myfont,
                    title = string("Unadjusted") , 
                    titlefontsize = mytitlefontsize,
                    yaxis = true, 
                    right_margin = -10mm, 
                    size=(400,300),  
            )
ptg_cdbcat_adjusted = confidenceplot(
                        vec(permutedims(CoefZcdbcat[idxColXmatCI, 2:end])), 
                        levelsCdb_cat[2:end],
                        vec(permutedims(CIZcdbcat[idxColXmatCI, 2:end])),
                        xlabel = namesX[idxColXmatCI]*" Effect Size", 
                        # xlim = (-0.25,0.25),
                        title = string("Adjusted") , 
                        legend = false,
                        fontfamily = myfont,
                        titlefontsize = mytitlefontsize,
                        yaxis = false, 
                        left_margin = -20mm,
)

p_ua = plot(p_ci_unadjusted, ptg_cdbcat_adjusted, size = (750, 500)) 

savefig("../../images/statin_mlmCI_unadjusted_Fish_Oil_carbon_count_and_degree_of_unsaturation_fomodel.svg")

p_ua
# -

# ## Enrichment Analysis

# ### Over Representation Analysis

# #### Compute P-values for each triglycerides

TstatZI

tstat_fishoil = TstatZI[2,:];

sampleN = 98
pval_fishoil = ccdf.(TDist(sampleN-1), abs.(TstatZI[2,:])).*2;

histogram(pval_fishoil, bin= :rice, legend = false, grid = false)

# +
# Hypergeometric(s, f, n)
#Determine
# number of metbolites of interest in a set 
# total of metbolites in a set 
# total of metabolites not in a set
# total number of metabolites of interest
# -

function ora_results(predictor_pval, myZgroups, mysubsets, dfMetabolo, dfAnnotation; thresh = 0.05)

    # Initiate results dataframe
    dfORA = DataFrame(
        Subset = mysubsets, 
        Pval = Vector{Float64}(undef, length(mysubsets))
    );

    # Find index TG with pval less than threshold value , default is 0.05
    idx_sgnf_tg = findall(predictor_pval .< thresh); 
    # Get names of the TGs
    set_sgnf_tg = names(dfMetabolo)[idx_sgnf_tg]

    # Get metbolites population size
    size_metabolites = size(dfMetabolo, 2)
    # Number of success in population
    size_significant = idx_sgnf_tg |> length

    for i in 1:length(mysubsets)
        idx_subset = findall(myZgroups .== mysubsets[i]);
        subset_tg = dfAnnotation.lipID[idx_subset]

        # Sample size
        size_subset = idx_subset |> length
        # Number of successes in sample
        size_significant_subset = intersect(set_sgnf_tg, subset_tg) |> length

        # Hypergeometric distribution
        # mydist = Hypergeometric(
        #             size_subset, 
        #             size_metabolites-size_subset,
        #             size_significant
        #         )
    
        mydist = Hypergeometric(
            size_significant,        
            size_metabolites - size_significant,
            size_subset
            )
        # Calculate p-value
        dfORA.Pval[i] = ccdf(mydist, size_significant_subset)
    end

    return dfORA
end

# #### Carbon

dfORA_tc = ora_results(pval_fishoil, dfRefTG.Total_C_cat, lvls_c, dfY, dfRefTG)

# #### Double Bonds

dfORA_db = ora_results(pval_fishoil, dfRefTG.Total_DB_cat, lvls_db, dfY, dfRefTG)

pval_fishoil_cat_adjusted = ccdf.(TDist(sampleN-1), abs.(TstatZcdbcat[2,:])).*2

dfORA_compare = vcat(dfORA_tc[2:end, :], dfORA_db[2:end, :])
dfORA_compare.Pval_MLM_adjusted = pval_fishoil_cat_adjusted[2:end];
dfORA_compare

pval_fishoil_db_cat_unadjusted = ccdf.(TDist(sampleN-1), abs.(TstatZdb[2,:])).*2;
pval_fishoil_tc_cat_unadjusted = ccdf.(TDist(sampleN-1), abs.(TstatZtc[2,:])).*2;


dfORA_compare.Pval_MLM_unadjusted = vcat(
    pval_fishoil_tc_cat_unadjusted[2:end],
    pval_fishoil_db_cat_unadjusted[2:end]
    )  

dfORA_compare







