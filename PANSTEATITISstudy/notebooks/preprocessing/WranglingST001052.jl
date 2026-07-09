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

# # Wrangling PANSTEATITIS study ST001052
# ---

# This notebook carries out the wrangling process for the [LIVER study ST001052 lipidomics data](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001052&StudyType=MS&ResultType=1) [1].

# ## Libraries

using Pkg 

Pkg.activate(joinpath(@__DIR__, "..", ".."))

Pkg.instantiate()

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

using DataFrames, CSV
using FreqTables #, CategoricalArrays
using StatsBase
using RCall# Conda, PyCall
using MetabolomicsWorkbenchAPI

# ## Ext. Functions

include(joinpath(@__DIR__,"..","..","src","wrangling_utils.jl" ));
include(joinpath(@__DIR__,"..","..","src","demog.jl" ));

# ## Load data ST001052

ST  = "ST001052";

# ## Extract clinical covariates

# Use the Julia's API to get the samples data from the [metabolomics workbench](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001052).

# get clinical covariates
dfIndividuals =  fetch_samples(ST);
print_df_size(dfIndividuals)

# List of the covariate names: 

names(dfIndividuals)

println("From the study description, $(ST) has $(fetch_total_subjects(ST)) subjects.")

# The clinical covariates dataframe contains 2 extra rows. We need to indicate what values corresponds to the `missing` data. In our case, all "-" will be replaced by `missing`.    

# assign missing value to "-"
dfIndividuals = ifelse.(dfIndividuals .== "-", missing, dfIndividuals);

# Check number of missing per columns.

print_variables_missing(dfIndividuals)

# ### Clinical dictionary

fileClinicalDict = joinpath(@__DIR__,"..","..","data","processed", "ClinicalDataDictionary.csv");
open(fileClinicalDict,"w") do io
   println(io,
        "Variable name, Variable description\n",
        "Gender,Sex\n",
        "Age, Years\n",
        "Weight (KG), Kilogram\n",
        "Length (CM), Centimeter\n",
        "Annuli,Number of opaque zones on fish scales\n",
        # "TG (CM),???\n",
        "VET SCORE,Veterinarian score where vet score < 1 indicates healthy tilapia and  vet score ≥ 1 indicates pansteatitis-affected tilapia.\n" ,
        "PCV Color,Pigmentation visually observed\n",
        "PCV,Pigmentation concentration volume\n",
        "Histology Adipose,Histological examination score of the adipose tissue\n",
        "Histology Liver,Histological examination score of the liver tissue\n",
        "Histology Swim Bladder,Histological examination score of the swim bladder tissue"
    )
end

# ### Independent variables

# Select variables of interest:

select!(dfIndividuals, Symbol.(["Sample ID",
                                "Group",
                                "Gender",
                                "Annuli",
                                "Age",
                                "WEIGHT (KG)",
                                "LENGTH (CM)",
                                "Histology Adipose",
                                # "Histology Liver",
                                # "Histology Swim Bladder",
]));

# Rename variables if needed:

rename!(dfIndividuals, Dict(
        :Group => "Status",
        :Gender => "Sex",
        Symbol("Sample ID") => "SampleID",
        Symbol("WEIGHT (KG)") => "Weight",
        Symbol("LENGTH (CM)") => "Length",
        Symbol("Histology Adipose") => "Histological_Score",
        # Symbol("Histology Liver") => "Histology_Liver",
        # Symbol("Histology Swim Bladder") => "Histology_Swim_Bladder",
));

# #### Extract histology score   
#
# The histology score spans between 0 and 5. Let write the histology dictionnary:   

fileHystologyDict = joinpath(@__DIR__,"..","..","data","processed", "HistologyDictionary.csv");
open(fileHystologyDict,"w") do io
   println(io,
        "Score, Histological Score Progression\n",
        "0, no signs\n",
        "1, minimal\n",
        "2, mild/few\n",
        "3, moderate",
        "4, moderate/severe\n",
        "5, severe\n",
    )
end

# Create a function to extract each score:

function xtrcHistScore(vecHist) 
    idxnotmissing = findall(.!(ismissing.(vecHist)));
    vScore = Vector{Union{Missing, Int}}(undef, length(vecHist));
    for i in idxnotmissing
        if isdigit(vecHist[i][1])
           vScore[i] = parse(Int, vecHist[i][1]) 
        end
    end
    
    return vScore
end;

# Replace original histology text by extracted score:

dfIndividuals.Histological_Score .= xtrcHistScore(dfIndividuals.Histological_Score);
# dfIndividuals.Histology_Liver .= xtrcHistScore(dfIndividuals.Histology_Liver);
# dfIndividuals.Histology_Swim_Bladder .= xtrcHistScore(dfIndividuals.Histology_Swim_Bladder);

# Filter incomplete cases:

# +
# filter complete cases
idxComplete = findall(completecases(dfIndividuals))
dfIndividuals = dfIndividuals[idxComplete, :]

# add a prefix to the ID samples
dfIndividuals.SampleID = "ID_".*string.(dfIndividuals.SampleID);

first(dfIndividuals, 5)
# -

# Insert a `GroupStatus` variable. The `Group` variable includes the diseases status and gender:

unique(dfIndividuals.Status)

# Let redefine the `Status` variable:

# +
# insertcols!(dfIndividuals, 3, :GroupStatus => occursin.("D", dfIndividuals.Group));
idxDiseased = findall(occursin.("D", dfIndividuals.Status)) ;
idxHealthy = findall(occursin.("H", dfIndividuals.Status));

dfIndividuals.Status[idxDiseased] .= "Diseased";
dfIndividuals.Status[idxHealthy] .= "Healthy";

idxMale = findall(occursin.("M", dfIndividuals.Sex)) ;
idxFemale = findall(occursin.("F", dfIndividuals.Sex));

dfIndividuals.Sex[idxMale] .= "Male";
dfIndividuals.Sex[idxFemale] .= "Female";
# -

first(dfIndividuals, 5)

# #### Save processed individuals dataset:

fileIndividuals = joinpath(@__DIR__,"..","..","data","processed","ST001052_ClinicalCovariates.csv");
dfIndividuals |> CSV.write(fileIndividuals);

# ### Demography

dfDemog = getDemographicST001052()

fileDemog = joinpath(@__DIR__,"..","..","data","processed","Demog.csv");
dfDemog |> CSV.write(fileDemog);

# ## Extract Metabolite references

# get clinical covariates
dfRef =  fetch_metabolites(ST);
print_df_size(dfRef)

# List the name of available properties:

names(dfRef)

first(dfRef, 5)

# Create a metabolite ID and keep only name and ID:

dfRef.MetaboliteID = "MT" .* string.(10000 .+ collect(1:size(dfRef, 1)));
select!(dfRef, [:Metabolite, :MetaboliteID]);

# ### Get Classification information

# To get the classifciation information, we use the package `MetabolomicsWorkbenchAPI.jl`.

dfClassification = fetch_properties(dfRef.Metabolite);
# insertcols!(dfClassification, 1, :metabolite_name => dfRef.metabolite_name);
first(dfClassification, 3)

dfRef |> names

dfClassification |> names

# replace "-" by missing
for c ∈ eachcol(dfClassification)
           replace!(c, "-" => missing)
end

idxmissing = findall(ismissing.(dfClassification.main_class))
dfRef.Metabolite[idxmissing]

size(dfRef)

# To be able to use all the lipid, especially the Triglycerides, we need to adjust the name in a more standardized way to be able to extract their properties information.  
# - The *Ox-* prefix mean oxidized, such as in *OxTG(16:0_20:5_20:3(OH))*. The *OH* indicates that it is a [TG hydroperoxide](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6550225/) [2].   
# - The *P-* and the *O-* indicate respectively if the lipid is a plasmalogen or a plasmanyl, such as *plasmanyl-TG(O-20:0_18:0_18:4)* and *plasmenyl-TG(P-20:1_15:0_16:0)*[3].
# - The non-oxidized lipids that contains Oxygen are described in [4].
#

dfRefOriginal = copy(dfRef);

dfRef = copy(dfRefOriginal);

dfRef.StandardizedName = copy(dfRef.Metabolite);

function standardizename(df::DataFrame, colname::String, matchstring)
    if !([colname] ⊆ names(df))
        df[:, colname] = repeat([false], size(dfRef, 1));
    end
    idx = findall(occursin.(matchstring, dfRef.Metabolite));
    df[idx, colname] .= true;
    # standardize name 
    df.StandardizedName[idx] .= replace.(dfRef.StandardizedName[idx], matchstring=>""); 
    return df
end

newcolname = ["Oxidized", "OH", "O", "O₂", "O₃", "O₄", "Plasmanyl", "Plasmalogen",
              "O₂", "O₃", "Plasmanyl", "Plasmalogen", "CHO", "Ke", "O" ] 
rmvstring = ["Ox", "(OH)", "+O", "(OO)", "(OOO)", "(OOOO)", "O-", "P-",
             "+OO", "+3O", r"(?i)plasmanyl-", r"(?i)plasmenyl-", "(CHO)", "(Ke)", "O" ]
for i in 1:length(newcolname)
    dfRef = standardizename(dfRef, newcolname[i], rmvstring[i]);
end    

dfRef

dfClassification = fetch_properties(dfRef.StandardizedName);
insertcols!(dfClassification, 1, :Metabolite => dfRef.Metabolite);
first(dfClassification, 3)

# replace "-" by missing
for c ∈ eachcol(dfClassification)
           replace!(c, "-" => missing)
end

idxmissing = findall(ismissing.(dfClassification.main_class))
dfRef.StandardizedName[idxmissing]

# At this stage, only *DMPE(16:0_22:6)* can be processed by [goslin](https://apps.lifs-tools.org/goslin/).   
# We will filter *PMe(16:0/18:1)* (phosphatidylmethanol) [5] and *ZyE(22:5)*.
#

dfRef |> names

dfClassification |> names

dfRef = leftjoin(dfRef, dfClassification, on = :Metabolite); size(dfRef)

# filter
# deleteat!(dfRef, idxmissing[[2,3]]);
deleteat!(dfRef, idxmissing[[5,11]]);

dfRef.sub_class |> unique

?skipmissing

findall(occursin.("PE-NMe2", skipmissing(dfRef.sub_class)))

skipmissing(dfRef.sub_class)[29]

# ### Use GOSLIN

R"""
suppressMessages(library('rgoslin'))
suppressMessages(library('tidyverse'));
"""

@rput dfRef;

R"""
# check validity
dfRef$Valid <- suppressWarnings(sapply(dfRef$StandardizedName, isValidLipidName))
if (sum(dfRef$Valid) == dim(dfRef)[1]) {
    cat("All valid.")
} else {
    print("Check invalid names.")
}
""";


dfRef.Total_C = zeros(Int, size(dfRef,1));
dfRef.Total_DB = zeros(Int, size(dfRef,1));
dfRef.Class = repeat(["NA"], size(dfRef,1));

for i in 1:size(dfRef, 1)
    @rput i;
    R"""
    #rsltGoslin <- as_tibble(parseLipidNames(dfRef$StandardizedName[i]))[, c("Original Name", "Total C", "Total DB", "Lipid Maps Main Class")];
    rsltGoslin <- as_tibble(parseLipidNames(dfRef$StandardizedName[i]))[, c("Original.Name", "Total.C", "Total.DB", "Lipid.Maps.Main.Class")];
    """
    @rget rsltGoslin
    # dfRef.Total_C[i] = parse(Int, rsltGoslin."Total C"[1])
    # dfRef.Total_DB[i] = parse(Int, rsltGoslin."Total DB"[1])
    # dfRef.Class[i] = rsltGoslin."Lipid Maps Main Class"[1]
        
    dfRef.Total_C[i] = rsltGoslin."Total_C"[1]
    dfRef.Total_DB[i] = rsltGoslin."Total_DB"[1]
    dfRef.Class[i] = rsltGoslin."Lipid_Maps_Main_Class"[1]
end

names(dfRef)

# #### Save metabolites reference dataset:

fileMetaboRef = joinpath(@__DIR__,"..","..","data","processed","refMeta.csv");
dfRef |> CSV.write(fileMetaboRef);

# ### Sub class dictionary

fileClinicalDict = joinpath(@__DIR__,"..","..","data","processed", "SubClassDictionary.csv");
open(fileClinicalDict,"w") do io
   println(io,
        "SubClass,Name\n",
        "Cer,Ceramides\n",
        "Chol. esters,Cholesteryl esters\n",
        "DAG,Diglycerides\n",
        "LPC,Lysophosphatidylcholines\n",
        "LPE,Lysocephalins\n",
        "PC,Phosphatidylcholines\n",
        "PE,Phosphatidylethanolamines\n" ,
        "PI,Phosphatidylinositols\n",
        "PS,Phosphatidylserines\n",
        "SM,Sphingomyelins\n",
        "TAG,Triglycerides\n"
    )
end

# ## Extract Metabolites dataset 

dfMetabo = fetch_data(ST);

# rename sample ID with suffix
vHeader = names(dfMetabo);
vHeader[2:end] .= "ID_".*vHeader[2:end];
rename!(dfMetabo, Symbol.(vHeader));

first(dfMetabo, 5)

# Replace `Metabolite` name information with `MetaboliteID` values:

dfMetaboAll = leftjoin(select(dfRefOriginal, [:Metabolite, :MetaboliteID]), dfMetabo, on = [:Metabolite]);
select!(dfMetaboAll, Not([:Metabolite]));

# Select the samples that only present in the filtered clinical dataset, `dfIndividuals`:   

select!(dfMetaboAll, vcat([:MetaboliteID], Symbol.(dfIndividuals.SampleID)));
size(dfMetaboAll)

# #### Save metabolites levels dataset:

fileMetabo = joinpath(@__DIR__,"..","..","data","processed","Metabo.csv");
dfMetaboAll = permutedims(dfMetaboAll, 1, :SampleID);
dfMetaboAll |> CSV.write(fileMetabo);

# ## References
#
# [1] Koelmel, J. P., Ulmer, C. Z., Fogelson, S., Jones, C. M., Botha, H., Bangma, J. T., Guillette, T. C., Luus-Powell, W. J., Sara, J. R., Smit, W. J., Albert, K., Miller, H. A., Guillette, M. P., Olsen, B. C., Cochran, J. A., Garrett, T. J., Yost, R. A., & Bowden, J. A. (2019). Lipidomics for wildlife disease etiology and biomarker discovery: a case study of pansteatitis outbreak in South Africa. Metabolomics : Official journal of the Metabolomic Society, 15(3), 38. https://doi.org/10.1007/s11306-019-1490-9    
#
# [2] Kato, S., Shimizu, N., Hanzawa, Y., Otoki, Y., Ito, J., Kimura, F., Takekoshi, S., Sakaino, M., Sano, T., Eitsuka, T., Miyazawa, T., & Nakagawa, K. (2018). Determination of triacylglycerol oxidation mechanisms in canola oil using liquid chromatography-tandem mass spectrometry. NPJ science of food, 2, 1. https://doi.org/10.1038/s41538-017-0009-x    
#
# [3] Koelmel, J. P., Ulmer, C. Z., Jones, C. M., Yost, R. A., & Bowden, J. A. (2017). Common cases of improper lipid annotation using high-resolution tandem mass spectrometry data and corresponding limitations in biological interpretation. Biochimica et biophysica acta. Molecular and cell biology of lipids, 1862(8), 766–770. https://doi.org/10.1016/j.bbalip.2017.02.016    
#
# [4] Riewe, D., Wiebach, J., & Altmann, T. (2017). Structure Annotation and Quantification of Wheat Seed Oxidized Lipids by High-Resolution LC-MS/MS. Plant physiology, 175(2), 600–618. https://doi.org/10.1104/pp.17.00470    
#
# [5] Koelmel, J. P., Jones, C. M., Ulmer, C. Z., Garrett, T. J., Yost, R. A., Schock, T. B., & Bowden, J. A. (2018). Examining heat treatment for stabilization of the lipidome. Bioanalysis, 10(5), 291–305. https://doi.org/10.4155/bio-2017-0209    
#
#

