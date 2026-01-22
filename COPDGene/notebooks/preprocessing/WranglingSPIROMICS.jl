# -*- coding: utf-8 -*-
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

# # Wrangling SPIROMICS
# ---

# This notebook carrieds out the wrangling process for the SPIROMICS metabolomics data.

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

using DataFrames, CSV, Missings
using FreqTables #, CategoricalArrays
using Statistics

# ### Ext. Functions

include(joinpath(@__DIR__,"..","..","src","wrangle_utils.jl" ));

# ### Load data ST001639: SPIROMICS

# #### Participants

fileIndividuals = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","SPIROMICS_ClinicalCovariates.csv"))
dfIndividuals = CSV.read(fileIndividuals, DataFrame;  delim = ',', missingstring = "NA");
first(dfIndividuals, 3)

unique(dfIndividuals.GOLD_STAGE_COPD_SEVERITY)

# #### Metabolites References

# The metabolites references are split in to 4 files:
# - neg_metadata.txt
# - polar_metadata.txt
# - pos_early_metadata.txt
# - pos_late_metadata.txt

# +
# Reference metabolomics
vRefMetaFiles = repeat(["path"],4)

vRefMetaFiles[1] = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS",
                                  "neg_metadata.txt"))
vRefMetaFiles[2] = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS",
                                  "polar_metadata.txt"))
vRefMetaFiles[3] = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS",
                                  "pos_early_metadata.txt"))
# Missing 1st column name: metabolite_name, in pos_polar_metadata.txt   
vRefMetaFiles[4] = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS",
                                  "pos_late_metadata_revised20220829.txt"))#"pos_polar_metadata.txt"))
prepend(vRefMetaFiles[4], "metabolite_name	") # prepend "metabolite_name  "
vRefMetaFiles[4] = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS",
                                  "new_pos_late_metadata_revised20220829.txt"))#"new_pos_polar_metadata.txt"))

# Generate dataframe
dfRefMetabo = CSV.read(vRefMetaFiles[1], DataFrame;  delim = '	');

for i in 2:length(vRefMetaFiles)
    dfRefMetabo = vcat(dfRefMetabo, CSV.read(vRefMetaFiles[i], DataFrame;  delim = '	'));
end
# add code to remove prepend new file
# rm(vRefMetaFiles[4])
first(dfRefMetabo, 1)
# -

names(dfRefMetabo)

# #### Negative

# **Notes:** *The metabolics files do not contains a column header for the name of the metabolites; we need to add a column name.*

fileNegMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","neg.txt"))
prepend(fileNegMetabo, "metabolite_name	")
fileNegMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","new_neg.txt"))
dfNegMetabo = CSV.read(fileNegMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(fileNegMetabo);
first(dfNegMetabo, 3)

# #### Polar

filePolarMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","polar.txt"))
prepend(filePolarMetabo, "metabolite_name	")
filePolarMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","new_polar.txt"))
dfPolarMetabo = CSV.read(filePolarMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePolarMetabo);
first(dfPolarMetabo, 3)

# #### Positive early

# + tags=[]
filePosEarlyMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","pos_early.txt"))
prepend(filePosEarlyMetabo, "metabolite_name	")
filePosEarlyMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","new_pos_early.txt"))
dfPosEarlyMetabo = CSV.read(filePosEarlyMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePosEarlyMetabo)
first(dfPosEarlyMetabo, 3)
# -

# #### Positive late

# + tags=[]
filePosLateMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","pos_late.txt"))
prepend(filePosLateMetabo, "metabolite_name	")
filePosLateMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","SPIROMICS","new_pos_late.txt"))
dfPosLateMetabo = CSV.read(filePosLateMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePosLateMetabo)
first(dfPosLateMetabo, 3)
# -

# #### Check that each dataframe contains the same column names(sample names). 

names(dfNegMetabo) == names(dfPolarMetabo) == names(dfPosEarlyMetabo) == names(dfPosLateMetabo)

# ## ST001639: SPIROMICS exploration

# + [markdown] tags=[]
# ### Individuals
# -

println("The participants dataset contains $(size(dfIndividuals, 1)) individuals and $(size(dfIndividuals, 2)) covariates.")

describe(dfIndividuals)

# Notes:`ATS_PackYears`, i.e. smoking number of pack-years, and `Insp_LAA950_total_Thirona`, i.e. percent emphysema, contain `NA`.

# Check how many `NA`:

vMissing = map(eachcol(dfIndividuals)) do col
               sum(ismissing.(col))
           end
idxColMiss = findall(vMissing .!= 0)
for i in idxColMiss
    println("$(names(dfIndividuals)[i]) contains $(vMissing[i]) missing values.")
end

# #### Change variable names: 

rename!(dfIndividuals, Dict(:SUBJID => "SampleName", :SAMPLE_NAME => "SampleID",
                            :SITE => "Site", :GOLD_STAGE_COPD_SEVERITY => "FinalGold",
                            :GENDER => "Sex", :AGE_DERV_01 => "Age", :RACE => "race",
                            :SMOKING_PACK_YEARS01 => "SmokingPackYears", :BMI_CM01 => "BMI",
                            :POST_FEV1FVC_DERV => "FEV1_FVC", :PCT_POST_FEV1_V1 => "FEV1pp",
                            :CURRENT_SMOKER_V1 => "CurrentSmoker",
                            :V1_PERCENT_EMPHYSEMA_TOTAL => "PercentEmphysema"));

# Verify how many factors per categorical variables, *i.e.* gender, race, smoking status, COPD case status.

vCovariateNames = [:Sex, :race, :CurrentSmoker, :COPD]
vUniqueCat = map(eachcol(dfIndividuals[:, vCovariateNames])) do col
                 join(unique(col), ",", " and ")
             end
for i in 1:length(vUniqueCat)
    println("$(string(vCovariateNames[i])) variable contains: $(vUniqueCat[i]) values.")
end


# Convert:
# - `race` into `NHW`, where the value *1* corresponds to non-Hispanic White and *0* otherwise.
# - `smoking_status` into `CurrentSmoker`, , where the value *1* corresponds to Current Smoker and *0* to Former Smoker.

#  non-Hipanic White
vNHW = zeros(Int, size(dfIndividuals, 1));
idxNHW = findall(dfIndividuals.race .== "Non-Hispanic, White");
vNHW[idxNHW] .= 1
dfIndividuals.NHW = (vNHW);

# Drop `race` and `smoking_status`
select!(dfIndividuals, Not([:race]));

# Get demographics SPIROMICS cohort by sex.

# +
# Group by sex
gdf = groupby(dfIndividuals, :Sex);

# Get mean values
mymean(X) = mean(skipmissing(X))
df1a = combine(gdf, [:Age, :BMI, :SmokingPackYears, :PercentEmphysema] .=> mymean)
df1a[:,2:end] = round.((df1a[:,2:end]); digits = 1)
rename!(df1a, Dict(:Age_mymean => "Age", :BMI_mymean => "BMI",
                   :SmokingPackYears_mymean => "SmokingPackYears", :PercentEmphysema_mymean => "PercentEmphysema"));

# Get standard deviation values
mystd(X) = std(skipmissing(X))
df1b = combine(gdf, [:Age, :BMI, :SmokingPackYears, :PercentEmphysema] .=> mystd)
df1b[:,2:end] = round.((df1b[:,2:end]); digits = 1)
rename!(df1b, Dict(:Age_mystd => "Age", :BMI_mystd => "BMI",
                   :SmokingPackYears_mystd => "SmokingPackYears", :PercentEmphysema_mystd => "PercentEmphysema"));

# Join mean and standard deviation values
dfDem1 = string.(df1a[:,2:end]).*repeat(["("], size(df1a,1),size(df1a,2)-1).* 
         string.(df1b[:,2:end]).*repeat([")"], size(df1a,1),size(df1a,2)-1);
insertcols!(dfDem1, 1, :Sex => df1a.Sex, :Participants => combine(gdf, nrow)[:,2])

# Get sum values
sum_skipmissing(x)= sum(skipmissing(x));
df2a = combine(gdf, [:NHW, :CurrentSmoker, :COPD] .=> sum_skipmissing)

# Get number of participant excluding those with missing values
nrow_skipmissing(x)= length(collect(skipmissing(x)));
dfParticipantSkipmissing = combine(gdf, [:NHW, :CurrentSmoker, :COPD] .=> nrow_skipmissing)


# Get percentage values
df2b = round.((df2a[:,2:end]./ Matrix(dfParticipantSkipmissing[:, 2:end])).*100, digits = 1)

# Join sum and percentage values
dfDem2 = string.(df2a[:,2:end]).*repeat(["("], size(df2a,1),size(df2a,2)-1).* 
         string.(df2b[:,1:end]).*repeat([")"], size(df2a,1),size(df2a,2)-1)
insertcols!(dfDem2, 1, :Sex => df2a.Sex)

rename!(dfDem2, Dict(:NHW_sum_skipmissing => "NHW", :CurrentSmoker_sum_skipmissing => "CurrentSmoker",
                   :COPD_sum_skipmissing => "COPD"))
# Join demographics dataframes
dfDem = leftjoin(dfDem1, dfDem2, on = :Sex )

# Pivot table
dfDem = permutedims(dfDem, 1, "Variable")
# -

# The demographic table is identical to the article "*Metabolomic Profiling Reveals Sex Specific Associations with
# Chronic Obstructive Pulmonary Disease and Emphysema*"(2021).

# #### Save processed individuals dataset:

first(dfIndividuals)

fileIndividuals = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","SPIROMICS_ClinicalCovariates.csv");
dfIndividuals |> CSV.write(fileIndividuals);

# + [markdown] tags=[]
# ### Metabolomics References
# -

# #### Create dataframe whith pathways

names(dfRefMetabo)

# Since the [PUBCHEM ID list](https://www.metabolomicsworkbench.org/rest/study/study_id/ST001639/metabolites/text) is incomplete, we need to generate a pseudo CompID. First sort by metabolite name then generate the IDs.

sort!(dfRefMetabo, :metabolite_name)
insertcols!(dfRefMetabo, 2, :CompID => replace.(string.(collect(1:size(dfRefMetabo, 1)).+99990000), "9999"=>"comp"));

first(dfRefMetabo, 3)

# Keep metabolite_name, comp ID, super pathways and sub pathways.

first(dfRefMetabo[:, Symbol.(["metabolite_name", "CompID", "SUPER.PATHWAY", "SUB.PATHWAY"])], 3)

last(dfRefMetabo[:, Symbol.(["metabolite_name", "CompID", "SUPER.PATHWAY", "SUB.PATHWAY"])], 3)

# Select variables of interest and rename accordingly
rename!(dfRefMetabo, Dict(Symbol("SUB.PATHWAY") => "SubPathway", Symbol("SUPER.PATHWAY") => "SuperPathway")) 
select!(dfRefMetabo, [:metabolite_name, :CompID, :SubPathway, :SuperPathway]);

# +
# Create 2 new variables name SubClassID and SuperClassID that 
# contain a codification of pathways

# Group by Super Pathway
gdf = groupby(dfRefMetabo, :SuperPathway);

nTotalSub = length(unique(dfRefMetabo.SubPathway))
vInit = repeat(["NA"], nTotalSub);
dfNewRef = DataFrame(SubPathway = vInit, SubClassID = vInit,
                     SuperPathway = vInit, SuperClassID = vInit);

# +
# Generate pathway ID references for the metabolites
idxStart = 1

for i in 1:(length(gdf)-1)
    vSub = unique(gdf[i].SubPathway)
    nSub = length(vSub)
    
    idxEnd = idxStart + nSub - 1
    
    dfNewRef.SubPathway[idxStart:idxEnd] = vSub;
    dfNewRef.SubClassID[idxStart:idxEnd] = uppercase(gdf[i].SuperPathway[1][1:3]).*string.(collect(1:nSub));
    dfNewRef.SuperPathway[idxStart:idxEnd] .= gdf[i].SuperPathway[1];
    dfNewRef.SuperClassID[idxStart:idxEnd] .= uppercase(gdf[i].SuperPathway[1][1:3]);
    
    idxStart = idxEnd + 1
end

# +
# Initiatlize vector
nMeta = size(dfRefMetabo, 1);
vSubClass = repeat(["NA"], nMeta);
vSupClass = repeat(["NA"], nMeta);

for i in 1:length(dfNewRef.SubPathway)
    idx = findall(dfRefMetabo.SubPathway.== dfNewRef.SubPathway[i])
    vSubClass[idx] .= dfNewRef.SubClassID[i]
    vSupClass[idx] .= dfNewRef.SuperClassID[i]
end
dfRefMetabo.SubClassID = vSubClass; 
dfRefMetabo.SuperClassID = vSupClass;


# Insert 0 in SubID when ID number less than 10. It helps for sorting.
idxSub2Change = findall(length.(dfRefMetabo.SubClassID) .== 4)
for i in 1:length(idxSub2Change) 
    dfRefMetabo.SubClassID[idxSub2Change[i]] = dfRefMetabo.SubClassID[idxSub2Change[i]][1:3]*"0"*dfRefMetabo.SubClassID[idxSub2Change[i]][4]
end
# -

# #### Check cotinine bio chemical

# The cotinine levels will be imputed differently if missing is more than 20%.

# check for cotinine
idxCotinine = findall(occursin.(r"(?i)cotinine", dfRefMetabo.metabolite_name))
dfRefMetabo[idxCotinine, :]

# #### Explore frequency table

freqtable(dfRefMetabo.SuperPathway)

# + tags=[]
idxLipid = findall(dfRefMetabo.SuperPathway .== "Lipid")
freqtable(dfRefMetabo[idxLipid, :SubPathway]); #  |> show;
# -

# #### Save processed individuals dataset:

fileRef = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","refMeta.csv");
dfRefMetabo |> CSV.write(fileRef);

# + [markdown] tags=[]
# ### Negative
# -

# **Notes:** *The metabolics files use the metabolite name to identify each metabolite. We need to use an ID similar to the compound ID from the COPDGenegene cohort. Since compound IDs are not available we use the workbench metabolite IDs.*

# Filter `dfNegMetabo` sample according to the individuals dataframe `dfIndividuals`:

# #### Keep complete cases

# + tags=[]
dfNegMetabo = keepComplete(dfNegMetabo, dfIndividuals, dfRefMetabo; sampleCol = :SampleID, metaCol = :metabolite_name );
# -

# #### Save filtered sample negative metabolites levels dataset:

fileNeg = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","negMeta.csv");
dfNegMetabo |> CSV.write(fileNeg);

println("The negative metabolite dataset contains $(size(dfNegMetabo, 2)-1) samples and $(size(dfNegMetabo, 1)) metabolites.")

# + [markdown] tags=[]
# ### Polar
# -

# Filter `dfPolarMetabo` sample according to the individuals dataframe `dfIndividuals`:

# #### Keep complete cases

# + tags=[]
dfPolarMetabo = keepComplete(dfPolarMetabo, dfIndividuals, dfRefMetabo; sampleCol = :SampleID, metaCol = :metabolite_name );
# -

# #### Save filtered sample polar metabolites levels dataset:

filePolar = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","polarMeta.csv");
dfPolarMetabo |> CSV.write(filePolar);

println("The polar metabolite dataset contains $(size(dfPolarMetabo, 2)-1) samples and $(size(dfPolarMetabo, 1)) metabolites.")

# + [markdown] tags=[]
# ### Positive Early
# -

# Filter `dfPosEarlyMetabo` sample according to the individuals dataframe `dfIndividuals`:
#
# [Polar molecules elute earlier and nonpolar molecules later.](https://www.sciencedirect.com/topics/immunology-and-microbiology/metabolome-analysis)

# #### Keep complete cases

# + tags=[]
dfPosEarlyMetabo = keepComplete(dfPosEarlyMetabo, dfIndividuals, dfRefMetabo; sampleCol = :SampleID, metaCol = :metabolite_name );
# -

# #### Save filtered sample positive early metabolites levels dataset:

filePosEarly = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","posEarlyMeta.csv");
dfPosEarlyMetabo |> CSV.write(filePosEarly);

println("The positive early metabolite dataset contains $(size(dfPosEarlyMetabo, 2)-1) samples and $(size(dfPosEarlyMetabo, 1)) metabolites.")

# + [markdown] tags=[]
# ### Positive Late
# -

# Filter `dfPosLateMetabo` sample according to the individuals dataframe `dfIndividuals`:
#
# [Polar molecules elute earlier and nonpolar molecules later.](https://www.sciencedirect.com/topics/immunology-and-microbiology/metabolome-analysis)

# #### Keep complete cases

# + tags=[]
dfPosLateMetabo = keepComplete(dfPosLateMetabo, dfIndividuals, dfRefMetabo; sampleCol = :SampleID, metaCol = :metabolite_name );
# -

# #### Save filtered sample positive late metabolites levels dataset:

filePosLate = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","posLateMeta.csv");
dfPosLateMetabo |> CSV.write(filePosLate);

println("The positive late metabolite dataset contains $(size(dfPosLateMetabo, 2)-1) samples and $(size(dfPosLateMetabo, 1)) metabolites.")

names(dfRefMetabo)

# +
# innerjoin(dfRefMetabo[:, [:metabolite_name, :CompID]], dfPosLateMetabo, on = :metabolite_name);

# +
# dfPosLateMetabo.metabolite_name ⊆ dfRefMetabo.metabolite_name

# +
# fileList = joinpath(@__DIR__,"..","..","data","processed","SPIROMICS","list.csv");
# DataFrame(Name = dfPosLateMetabo.metabolite_name) |> CSV.write(fileList);
# -


