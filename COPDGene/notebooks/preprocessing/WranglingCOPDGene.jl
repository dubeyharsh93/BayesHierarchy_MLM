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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# # Wrangling COPDGene
# ---

# This notebook carrieds out the wrangling process for the COPDGene metabolomics data.

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

# ### Load data ST001443: COPDGene

# #### Participants

fileIndividuals = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","COPDGene_ClinicalCovariates.csv"))
dfIndividuals = CSV.read(fileIndividuals, DataFrame;  delim = ',', missingstring = "NA");
first(dfIndividuals, 3)

# #### Metabolites References

# Reference metabolomics
fileRefMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene",
                                  "COPDGene_2018_Metabolon_metabolite_metadata.txt"))
prepend(fileRefMetabo, "MetaID	")
fileRefMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene",
                                  "new_COPDGene_2018_Metabolon_metabolite_metadata.txt"))
dfRefMetabo = CSV.read(fileRefMetabo, DataFrame;  delim = '	');
rm(fileRefMetabo);
names(dfRefMetabo)

# Notes: we added `MetaID` variable name.  

# #### Negative

fileNegMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","neg.txt"))
prepend(fileNegMetabo, "MetaID	")
fileNegMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","new_neg.txt"))
dfNegMetabo = CSV.read(fileNegMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(fileNegMetabo);
first(dfNegMetabo, 3)

# #### Polar

filePolarMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","polar.txt"))
prepend(filePolarMetabo, "MetaID	")
filePolarMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","new_polar.txt"))
dfPolarMetabo = CSV.read(filePolarMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePolarMetabo)
first(dfPolarMetabo, 3)

# #### Positive early

filePosEarlyMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","pos.early.txt"))
prepend(filePosEarlyMetabo, "MetaID	")
filePosEarlyMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","new_pos.early.txt"))
dfPosEarlyMetabo = CSV.read(filePosEarlyMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePosEarlyMetabo)
first(dfPosEarlyMetabo, 3)

# #### Positive late

filePosLateMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","pos.late.txt"))
prepend(filePosLateMetabo, "MetaID	")
filePosLateMetabo = realpath(joinpath(@__DIR__,"..","..","data","raw","COPDGene","new_pos.late.txt"))
dfPosLateMetabo = CSV.read(filePosLateMetabo, DataFrame;  delim = '	', missingstring = "NA");
rm(filePosLateMetabo)
first(dfPosLateMetabo, 3)

# ## ST001443: COPDGene exploration

# ### Individuals

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

rename!(dfIndividuals, Dict(:sid => "SampleName", :sample_name => "SampleID",
                            :ccenter => "Site", :finalgold_visit => "FinalGold",
                            :FEV1_FVC_utah => "FEV1_FVC", :FEV1pp_utah => "FEV1pp",
                            :gender => "Sex", :age_visit => "Age", 
                            :ATS_PackYears => "SmokingPackYears",
                            :Insp_LAA950_total_Thirona => "PercentEmphysema"));

# Verify how many factors per categorical variables, *i.e.* gender, race, smoking status, GOLD index, and COPD case status.

vCovariateNames = [:Sex, :race, :smoking_status, :FinalGold, :COPD]
vUniqueCat = map(eachcol(dfIndividuals[:, vCovariateNames])) do col
                 join(unique(col), ",", " and ")
             end
for i in 1:length(vUniqueCat)
    println("$(string(vCovariateNames[i])) variable contains: $(vUniqueCat[i]) values.")
end


# Convert:
# - `race` into `NHW`, where the value *1* corresponds to non-Hispanic White and *0* otherwise.
# - `smoking_status` into `CurrentSmoker`, where the value *1* corresponds to Current Smoker and *0* to Former Smoker.
# - `FinalGold` values are changed to  "0", "2", "3", and "4".

# +
#  non-Hipanic White
vNHW = zeros(Int, size(dfIndividuals, 1));
idxNHW = findall(dfIndividuals.race .== "White");
vNHW[idxNHW] .= 1
dfIndividuals.NHW = (vNHW);

# GOLD index
dfIndividuals.FinalGold = [match.(r"\d+", s).match for s in dfIndividuals.FinalGold];

# Current Smokers
vCurrentSmoker = zeros(Int, size(dfIndividuals, 1));
idxCurrentSmoker = findall(dfIndividuals.smoking_status .== "Current smoker");
vCurrentSmoker[idxCurrentSmoker] .= 1
dfIndividuals.CurrentSmoker = vCurrentSmoker;

# Drop `race` and `smoking_status`
select!(dfIndividuals, Not([:race, :smoking_status]));
# -

# Get demographics COPDGene cohort by sex.

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
df2a = combine(gdf, [:NHW, :CurrentSmoker, :COPD] .=> sum)

# Get percentage values
df2b = round.((df2a[:,2:end]./ dfDem1.Participants).*100, digits = 1)

# Join sum and percentage values
dfDem2 = string.(df2a[:,2:end]).*repeat(["("], size(df2a,1),size(df2a,2)-1).* 
         string.(df2b[:,1:end]).*repeat([")"], size(df2a,1),size(df2a,2)-1)
insertcols!(dfDem2, 1, :Sex => df2a.Sex)

rename!(dfDem2, Dict(:NHW_sum => "NHW", :CurrentSmoker_sum => "CurrentSmoker",
                   :COPD_sum => "COPD"))
# Join demographics dataframes
dfDem = leftjoin(dfDem1, dfDem2, on = :Sex )

# Pivot table
dfDem = permutedims(dfDem, 1, "Variable")
# -

# The demographic table is identical to the article "*Metabolomic Profiling Reveals Sex Specific Associations with
# Chronic Obstructive Pulmonary Disease and Emphysema*"(2021).

# #### Save processed individuals dataset:

first(dfIndividuals)

fileIndividuals = joinpath(@__DIR__,"..","..","data","processed","COPDGene","COPDGene_ClinicalCovariates.csv");
dfIndividuals |> CSV.write(fileIndividuals);

# ### Metabolomics References

# #### Create dataframe whith pathways

# Keep ID, biochemical, comp ID, super pathways and sub pathways.

first(dfRefMetabo[:, Symbol.(["MetaID", "BIOCHEMICAL", "COMP.ID"])], 3)

last(dfRefMetabo[:, Symbol.(["MetaID", "BIOCHEMICAL", "COMP.ID"])], 3)

# Select variables of interest and rename accordingly
rename!(dfRefMetabo, Dict(:BIOCHEMICAL => "Biochemical", Symbol("COMP.ID") => "CompID", 
                         Symbol("SUB.PATHWAY") => "SubPathway", Symbol("SUPER.PATHWAY") => "SuperPathway")) 
select!(dfRefMetabo, [:MetaID, :Biochemical, :CompID, :SubPathway, :SuperPathway]);

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
    vSub = sort(unique(gdf[i].SubPathway))
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
vClass = repeat(["NA"], nMeta);
vSupClass = repeat(["NA"], nMeta);

for i in 1:length(dfNewRef.SubPathway)
    idx = findall(dfRefMetabo.SubPathway.== dfNewRef.SubPathway[i])
    vClass[idx] .= dfNewRef.SubClassID[i]
    vSupClass[idx] .= dfNewRef.SuperClassID[i]
end
dfRefMetabo.SubClassID = vClass; 
dfRefMetabo.SuperClassID = vSupClass;
dfRefMetabo.CompID = "comp".*string.(dfRefMetabo.CompID);

# Insert 0 in SubID when ID number less than 10. It helps for sorting.
idxSub2Change = findall(length.(dfRefMetabo.SubClassID) .== 4)
for i in 1:length(idxSub2Change) 
    dfRefMetabo.SubClassID[idxSub2Change[i]] = dfRefMetabo.SubClassID[idxSub2Change[i]][1:3]*"0"*dfRefMetabo.SubClassID[idxSub2Change[i]][4]
end
# -

dfRefMetabo[2,Not([4,5])]

# #### Check cotinine bio chemical

# The cotinine levels will be imputed differently if missing is more than 20%.

# check for cotinine
idxCotinine = findall(occursin.(r"(?i)cotinine", dfRefMetabo.Biochemical))
dfRefMetabo[idxCotinine, :]

# #### Explore frequency table

freqtable(dfRefMetabo.SuperPathway)

idxLipid = findall(dfRefMetabo.SuperPathway .== "Lipid")
freqtable(dfRefMetabo[idxLipid, :SubPathway]); #  |> show;

# #### Save processed individuals dataset:

fileRef = joinpath(@__DIR__,"..","..","data","processed","COPDGene","refMeta.csv");
dfRefMetabo |> CSV.write(fileRef);

# ### Negative

# Filter `dfNegMetabo` sample according to the individuals dataframe `dfIndividuals`:

# #### Keep complete cases

dfNegMetabo = keepComplete(dfNegMetabo, dfIndividuals, dfRefMetabo; sampleCol=  :SampleID);

# #### Save filtered sample negative metabolites levels dataset:

fileNeg = joinpath(@__DIR__,"..","..","data","processed","COPDGene","negMeta.csv");
dfNegMetabo |> CSV.write(fileNeg);

println("The negative metabolite dataset contains $(size(dfNegMetabo, 2)-1) samples and $(size(dfNegMetabo, 1)) metabolites.")

# ### Polar

# Filter `dfPolarMetabo` sample according to the individuals dataframe `dfIndividuals`:

# #### Keep complete cases

dfPolarMetabo = keepComplete(dfPolarMetabo, dfIndividuals, dfRefMetabo; sampleCol=  :SampleID);

# #### Save filtered sample polar metabolites levels dataset:

filePolar = joinpath(@__DIR__,"..","..","data","processed","COPDGene","polarMeta.csv");
dfPolarMetabo |> CSV.write(filePolar);

println("The polar metabolite dataset contains $(size(dfPolarMetabo, 2)-1) samples and $(size(dfPolarMetabo, 1)) metabolites.")

# check name of Polar chemicals, they seem to be negative polar metablites
idxPol = findall(x -> [x] ⊆ dfPolarMetabo.CompID, dfRefMetabo.CompID);
dfRefMetabo.Biochemical[idxPol];

# Polar dataset contains negative polar metabolites.

# ### Positive Early

# Filter `dfPosEarlyMetabo` sample according to the individuals dataframe `dfIndividuals`:
#
# [Polar molecules elute earlier and nonpolar molecules later.](https://www.sciencedirect.com/topics/immunology-and-microbiology/metabolome-analysis)

# #### Keep complete cases

dfPosEarlyMetabo = keepComplete(dfPosEarlyMetabo, dfIndividuals, dfRefMetabo; sampleCol=  :SampleID);

# #### Save filtered sample positive early metabolites levels dataset:

filePosEarly = joinpath(@__DIR__,"..","..","data","processed","COPDGene","posEarlyMeta.csv");
dfPosEarlyMetabo |> CSV.write(filePosEarly);

println("The positive early metabolite dataset contains $(size(dfPosEarlyMetabo, 2)-1) samples and $(size(dfPosEarlyMetabo, 1)) metabolites.")

# ### Positive Late

# Filter `dfPosLateMetabo` sample according to the individuals dataframe `dfIndividuals`:
#
# [Polar molecules elute earlier and nonpolar molecules later.](https://www.sciencedirect.com/topics/immunology-and-microbiology/metabolome-analysis)

# #### Keep complete cases

dfPosLateMetabo = keepComplete(dfPosLateMetabo, dfIndividuals, dfRefMetabo; sampleCol=  :SampleID);

# #### Save filtered sample positive late metabolites levels dataset:

filePosLate = joinpath(@__DIR__,"..","..","data","processed","COPDGene","posLateMeta.csv");
dfPosLateMetabo |> CSV.write(filePosLate);

println("The positive late metabolite dataset contains $(size(dfPosLateMetabo, 2)-1) samples and $(size(dfPosLateMetabo, 1)) metabolites.")


