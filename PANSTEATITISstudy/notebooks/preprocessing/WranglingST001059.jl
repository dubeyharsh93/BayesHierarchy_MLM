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

# # Wrangling PANSTEATITIS study ST001059
# ---

# This notebook carries out the wrangling process for the [LIVER study ST001059 lipidomics data](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001059&StudyType=MS&ResultType=1) [1].

# ## Libraries

# To use RCall for the first time, one needs to 
# the location of the R home directory.
firstTimeRCall = false
if firstTimeRCall 
    using Pkg
    ENV["R_HOME"] = "C:/PROGRA~1/R/R-42~1.1" # from R.home() in R
    Pkg.build("RCall")
end     

using DataFrames, CSV
using FreqTables #, CategoricalArrays
using StatsBase
using Conda, RCall, PyCall
using MetabolomicsWorkbenchAPI

# ## Ext. Functions

include(joinpath(@__DIR__,"..","..","src","wrangling_utils.jl" ));
include(joinpath(@__DIR__,"..","..","src","demog.jl" ));

# ## Load data ST001059

ST  = "ST001059";

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

# **Notes:** In VET SCORE the character "-" actually indicates 0 or absence of pansteatitis.

# ### Clinical dictionary

fileClinicalDict = joinpath(@__DIR__,"..","..","data","processed", "ClinicalDataDictionary_ST001059.csv");
open(fileClinicalDict,"w") do io
   println(io,"Group,Indicates sex and disease status.\n",
        "Gender,Sex.\n",
        "Age, Years\n",
        "Weight (KG), Kilogram\n",
        "Length (CM), Centimeter\n",
        "Annuli,Number of opaque zones on fish scales.\n",
        # "TG (CM),???\n",
        "VET SCORE  (Adipose),Veterinarian score where vet score < 1 indicates healthy tilapia and  vet score ≥ 1 indicates pansteatitis-affected tilapia.\n" ,
        "TOTAL PROTEIN (g/100mL),Total protein measurement in grams per deciliter",
        "PCV Color,Pigmentation visually observed.\n",
        "PCV,Pigmentation concentration volume.\n",
        "Histology Adipose,Histological examination score of the adipose tissue.\n",
        "Histology Liver,Histological examination score of the liver tissue.\n",
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
                                "VET SCORE (Adipose)",
]));

# Rename variables if needed:

rename!(dfIndividuals, Dict(Symbol("Sample ID") => "SampleID",
                            Symbol("WEIGHT (KG)") => "Weight",
                            Symbol("LENGTH (CM)") => "Length",
                            Symbol("VET SCORE (Adipose)",) => "VetScore",
));

# Filter incomplete cases:

# Replace the missings in VetScore by 0
dfIndividuals.VetScore = coalesce.(dfIndividuals.VetScore, 0);

# +
# filter complete cases
idxComplete = findall(completecases(dfIndividuals))
dfIndividuals = dfIndividuals[idxComplete, :]

# add a prefix to the ID samples
dfIndividuals.SampleID = "ID_".*string.(dfIndividuals.SampleID);

first(dfIndividuals, 5)
# -

# Insert a `GroupStatus` variable. The `Group` variable includes the diseases status and gender:

unique(dfIndividuals.Group)

# Let redefine the `Group` variable:

# +
# insertcols!(dfIndividuals, 3, :GroupStatus => occursin.("D", dfIndividuals.Group));
idxDiseased = findall(occursin.("D", dfIndividuals.Group)) ;
idxHealthy = findall(occursin.("H", dfIndividuals.Group));

dfIndividuals.Group[idxDiseased] .= "Diseased";
dfIndividuals.Group[idxHealthy] .= "Healthy";

idxMale = findall(occursin.("M", dfIndividuals.Gender)) ;
idxFemale = findall(occursin.("F", dfIndividuals.Gender));

dfIndividuals.Gender[idxMale] .= "Male";
dfIndividuals.Gender[idxFemale] .= "Female";
# -

first(dfIndividuals, 5)

# #### Save processed individuals dataset:

fileIndividuals = joinpath(@__DIR__,"..","..","data","processed","ST001059_ClinicalCovariates.csv");
dfIndividuals |> CSV.write(fileIndividuals);

# ### Demography

dfDemog = getDemographicST001059("ST001059")

fileDemog = joinpath(@__DIR__,"..","..","data","processed","Demog_ST001059.csv");
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

"""
standardizename(df::DataFrame, colname::String, matchstring, newstring="")

Takes a dataframe that contains a column of standardized name, and 
create a new column name according to `colname` argument. The type
of this new column is boolean, where each entry is true if the 
standardize name value contains the matching string argument.
In addition, it replaces the matching string by a new string, default
is "".
"""
function standardizename(df::DataFrame, colname::String, matchstring, newstring="")
    if [colname] ⊈ names(df) #!([colname] ⊆ names(df))
        df[:, colname] = repeat([false], size(dfRef, 1));
    end
    idx = findall(occursin.(matchstring, dfRef.Metabolite));
    df[idx, colname] .= true;
    # standardize name 
    df.StandardizedName[idx] .= replace.(dfRef.StandardizedName[idx], matchstring=>newstring); 
    return df
end

newcolname = ["ProteinBound", "Hexosyl", 
              "CeramideAP", "CeramideAS", "CeramideNS", "CeramideNP",   "CeramideNDS",
              "Oxidized", "OH", "O", "O₂", "O₃",
              "O₄", "Plasmanyl", "Plasmalogen", "O₂", "O₂",
              "O₃", "O₄", "O₆", "Plasmanyl", "Plasmalogen", 
              "CHO", "Ke", "Ke_OH", "O", 
              "IgnoreCol",  
             ] 
rmvstring = ["+pO_", "CerG1", 
             "_AP", "_AS", "_NS", "_NP", "_NDS",
             "Ox", "(OH)", "+O", "(OO)", "(OOO)",
             "(OOOO)", "O-", "P-", "+OO", "+2O", 
             "+3O", "+4O", "+6O", r"(?i)plasmanyl-", r"(?i)plasmenyl-", 
             "(CHO)", "(Ke)", "(Ke,OH)", "O", 
             "_",
            ]
newstring = vcat( "/", "HexCer", repeat([""], 24), "/")
for i in 1:length(newcolname)
    dfRef = standardizename(dfRef, newcolname[i], rmvstring[i], newstring[i]);
end    

dfClassification = fetch_properties(dfRef.StandardizedName);
insertcols!(dfClassification, 1, :Metabolite => dfRef.Metabolite);
first(dfClassification, 3)
idxmissing = findall(ismissing.(dfClassification.main_class))
dfRef.StandardizedName[idxmissing]



# We need to replace "+pO_" or "+p_" by "/", the PO type ceramides are protein bound ceramide [7].
# The ceramides "CerG1" are hexosylceramides and need to be replaced by "HexCer"[6]   
# The ceramides "GlcCer" are glucosylceramide.
# The "OAHFA"s are (O-acyl) ω-hydroxy fatty acids.
#
# At this stage, only *DMPE(16:0_22:6)* can be processed by [goslin](https://apps.lifs-tools.org/goslin/).   
# We will filter *PMe(16:0/18:1)* (phosphatidylmethanol) [5] and *ZyE(22:5)*.
#

dfRef = leftjoin(dfRef, dfClassification, on = :Metabolite); size(dfRef)
# filter
deleteat!(dfRef, idxmissing[[2,3]]);

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
# [6] Muhammad Z. Chauhan, Paul H. Phillips, Joseph G. Chacko, David B. Warner, Daniel Pelaez, Sanjoy K. Bhattacharya,
# Temporal Alterations of Sphingolipids in Optic Nerves After Indirect Traumatic Optic Neuropathy, Ophthalmology Science, Volume 3, Issue 1, 2023, 100217, ISSN 2666-9145, https://doi.org/10.1016/j.xops.2022.100217.
#
# [7] Madoka Suzuki, Yusuke Ohno, Akio Kihara, Whole picture of human stratum corneum ceramides, including the chain-length diversity of long-chain bases, Journal of Lipid Research, Volume 63, Issue 7, 2022, 100235, ISSN 0022-2275, https://doi.org/10.1016/j.jlr.2022.100235.
#
 dfIndividuals =  fetch_samples(ST);

    # dfIndividuals = ifelse.(dfIndividuals .== "-", "0", dfIndividuals)
    dfIndividuals = ifelse.(dfIndividuals .== "-", missing, dfIndividuals);
    
    select!(dfIndividuals, Symbol.(["Sample ID",
                                "Group",
                                "Gender",
                                "Annuli",
                                "Age",
                                "WEIGHT (KG)",
                                "LENGTH (CM)",
    ]));
    
    # filter complete cases
    idxComplete = findall(completecases(dfIndividuals))
    dfIndividuals = dfIndividuals[idxComplete, :]
    
    insertcols!(dfIndividuals, 3, :GroupStatus => occursin.("D", dfIndividuals.Group));
    
    # Categorical Variables
    catVar = ["GroupStatus", "Gender"]
    vDFcat = Vector{DataFrame}(undef, length(catVar));
    
    for i in 1:length(catVar)
        dfCat = combine(groupby(dfIndividuals, [Symbol(catVar[i])]), nrow => :Count)
        sort!(dfCat, [Symbol(catVar[i])]) 
        dfCat = DataFrame(vcat([names(dfCat)[1] " "], Matrix(dfCat)), ["Clinical Features", "Count/ mean(SD)"])
        vDFcat[i] = dfCat 
    end
    


    # Continuous Variables
    contVar = ["Annuli", "Age", "WEIGHT (KG)", "LENGTH (CM)"]
    vDFcont = Vector{DataFrame}(undef, length(contVar));
   


 
    for i in 1:length(contVar)
    # i = 1
        # calculate mean
        vVar = dfIndividuals[:,Symbol(contVar[i])];
        idxNotMiss = findall(.!ismissing.(vVar));

        vVar = parse.(Float32, string.(vVar[idxNotMiss]))

        myMean = vVar |> mean |> x->round(x, digits = 2)
        # calculate SD
        myStd = vVar |> std |> x->round(x, digits = 2)
    
        dfCont  = DataFrame([contVar[i] string(myMean,"(", myStd, ")")], ["Clinical Features", "Count/ mean(SD)"])
        vDFcont[i] = dfCont 
    end
    vDF = vcat(vDFcat, vDFcont)
    
    dfDem = reduce(vcat, vDF)


findall(occursin.("D", dfIndividuals.Group)) 


