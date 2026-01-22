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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# # Preprocessing Step
# ---

# This notebook carrieds out the preprocessing steps for the metabolomics data.

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

using CSV, DataFrames, Missings, CategoricalArrays
using StatsBase, Statistics, MultivariateStats, RCall
using FreqTables, Plots, StatsPlots

# Change display options
Base.displaysize() = (280, 88)

# ### Ext. Functions

include(joinpath(@__DIR__,"..","..","src","preprocessing.jl" ));
include(joinpath(@__DIR__,"..","..","src","utils.jl" ));

# ### Load data

# +
# Load Negative
negMeta = realpath((@__DIR__)*"/../../data/data_primary/Elam_metabolites_2020_NEG_PH_cleaned_blkfil-batch-effect-analysis.csv")
dfNegMeta = DataFrame(CSV.File(negMeta));

negMetaKln = realpath((@__DIR__)*"/../../data/data_primary/Elam_metabolites-negative-2020_batch-corrected-all-data-elements.csv")
dfNegMetaKln = DataFrame(CSV.File(negMetaKln));

# Load Positive
posMeta = realpath(joinpath(@__DIR__,"..", "..","data","data_primary","Elam_metabolites_2020_POS_PH_cleaned_blkfil-batch-effect-analysis.csv" ))
dfPosMeta = DataFrame(CSV.File(posMeta));

posMetaKln = realpath(joinpath(@__DIR__,"..", "..","data","data_primary","Elam_metabolites-positive-2020_batch-corrected-all-data-elements.csv" ))
dfPosMetaKln = DataFrame(CSV.File(posMetaKln));
# -

# ## Chain of preprocessing

diGraph(["Wrangle", "Imputation", "Normalization", "Log2 Transformation", "Batch Effect Correction"])

# ## Wrangle data
# ----

# ### Metabolites reference dictionnary

# We create a reference dictionnary for the metabolites name.

first(dfNegMeta)

first(dfPosMeta)

# +
# isCossRef = false
dfNegMeta, dfNegCrossRef = crossRefLipids(dfNegMeta, "negMeta", startCol = 5);
dfPosMeta, dfPosCrossRef = crossRefLipids(dfPosMeta, "posMeta", startCol = 4);

# # save cross reference
# dfNegCrossRef |> CSV.write("../data/dataprocessed/inl2b_NegMeta_Xref.csv");
# dfPosCrossRef |> CSV.write("../data/dataprocessed/inl2b_PosMeta_Xref.csv");

# merge neg and pos look up table
dfMetaCrossRef = deepcopy(dfNegCrossRef);
append!(dfMetaCrossRef, dfPosCrossRef)
# dfMetaCrossRef |> CSV.write("../data/dataprocessed/inl2b_Meta_Xref.csv");
# -

dfPosCrossRef[949,:]

first(dfPosMeta)

# ### Create a Group variable: CN, CS, CSbaseline

# Adjust for column name according to primary csv file.

# drop CLASS column 
dfNegMetaKln = select(dfNegMetaKln, Not(:CLASS));
rename!(dfNegMetaKln, :NAME=> :Sample);
# rename CLASS to Batch and copy batch info from dfNegMeta
rename!(dfPosMetaKln, :CLASS=> :Batch);
rename!(dfPosMetaKln, :NAME=> :Sample);
dfPosMetaKln.Batch = dfNegMetaKln.Batch;

dfNegMetaKln = addGroupCatMeta(dfNegMetaKln, false);
dfPosMetaKln = addGroupCatMeta(dfPosMetaKln, false);

# Check group population:

countmap(dfNegMetaKln.Group; alg = :dict)

first(dfNegMetaKln)

# ### Convert to categorical: Batch, Statin and Oil

dfNegMetaKln = catBGSF(dfNegMetaKln, false);
dfPosMetaKln = catBGSF(dfPosMetaKln, false);

dfNegMetaKln[3:7,:]

dfPosMetaKln[3:7,:]

# ### Check tables

freqtable(dfNegMetaKln, :Group, :Batch)

freqtable(dfNegMetaKln, :Group, :FishOil)

# ## Get original values with batch effect 

# ### Keep only IDs' sample and values

# Keep only Sample IDs and metabolites'values of the original data an combine with `Group`, `Batch` and `Fish Oil` variables the corrected files.

# Negative:

first(dfNegMeta)

# Remove the columns: *Group, Batch and Column4*

select!(dfNegMeta, Not([:Group, :Batch, :Column4]));
first(dfNegMeta)

# Positive:

first(dfPosMeta)

select!(dfPosMeta, Not([:Group, :Batch]));
first(dfPosMeta)

# ### Check that metabolites column names are similar

# +
vNameNegMeta = names(dfNegMeta)[2:end];
vNamePosMeta = names(dfPosMeta)[2:end];

vNameNegMetaKln = names(dfNegMetaKln)[6:end];
vNamePosMetaKln = names(dfPosMetaKln)[6:end];
# -

# Check length:

DataFrame(Meta = ["Neg", "Pos"], 
      NotClean = [length(vNameNegMeta), length(vNamePosMeta)],
      Clean = [length(vNameNegMetaKln), length(vNamePosMetaKln)])

# ### Join info fom corrected data set to original values

first(dfNegMetaKln)

select!(dfNegMetaKln, [:Sample, :Batch, :Group, :Statin, :FishOil]);
select!(dfPosMetaKln, [:Sample, :Batch, :Group, :Statin, :FishOil]);

# Left-join to the original values:

dfNegMeta = leftjoin(dfNegMetaKln, dfNegMeta, on= :Sample)
first(dfNegMeta, 5)

dfPosMeta = leftjoin(dfPosMetaKln, dfPosMeta, on= :Sample);
first(dfPosMeta,5)

# ## Impute missing data  
# ---
# HM (Half of the Minimum): This method replaces missing elements with half of the minimum of non-missing elements in the corresponding variable.

dfNegMeta = imputeHM(dfNegMeta, startCol = 6);
dfPosMeta = imputeHM(dfPosMeta, startCol = 6);

findall(Matrix(dfNegMeta[:,6:end]) .=== missing)

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

dfNegMeta = pqnorm(dfNegMeta, startCol = 6);
dfPosMeta = pqnorm(dfPosMeta, startCol = 6);

# ## Transformation
# ---
#
# A simple and widely used transformation to make data more symmetric and homoscedastic is the log-transformation.

dfNegMeta = log2tx(dfNegMeta, startCol = 6);
dfPosMeta = log2tx(dfPosMeta, startCol = 6);

first(dfNegMeta)

# ## Adjusting for batch effects
# ---

# ### Check for batch effects

R"""
suppressMessages(library(mixOmics))
suppressMessages(library(tidyverse));
"""

dfNeg = catBGSF(dfNegMeta);

dfPos = catBGSF(dfPosMeta);

@rput dfNeg;
@rput dfPos;

R"summary(dfNeg[,c(1:5)])"

# +
# get matrix data
Xneg= copy(transpose(Matrix(dfNeg[:,6:end])));
Xpos= copy(transpose(Matrix(dfPos[:,6:end])));

# train a PCA model
Mneg = fit(PCA, Xneg; maxoutdim=10)
Mpos = fit(PCA, Xpos; maxoutdim=10)
# -


size(Xneg)

# get explained variance
explainedVarPCAneg = principalvars(Mneg)./tvar(Mneg)
explainedVarPCApos = principalvars(Mpos)./tvar(Mpos);

# +
ticklabel = string.(collect(1:10))
pNeg =bar(explainedVarPCAneg, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Negative Metabolites")
xlabel!("Principal Components")
ylabel!("Explained Variance");

pPos =bar(explainedVarPCApos, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Positive Metabolites")
xlabel!("Principal Components")
ylabel!("Explained Variance");

plot(pNeg, pPos, size = (800, 400))

# +
# get batch group labels
XbatchNeg = Vector(dfNeg[:,2]);
XbatchPos = Vector(dfPos[:,2]);

# apply PCA model 
Yneg = MultivariateStats.transform(Mneg, Xneg)
Ypos = MultivariateStats.transform(Mpos, Xpos)

# group results by testing set labels for color coding
B1neg = Yneg[:,XbatchNeg.=="B1"]; B1pos = Ypos[:,XbatchPos.=="B1"] 
B2neg = Yneg[:,XbatchNeg.=="B2"]; B2pos = Ypos[:,XbatchPos.=="B2"]
B3neg = Yneg[:,XbatchNeg.=="B3"]; B3pos = Ypos[:,XbatchPos.=="B3"]
B4neg = Yneg[:,XbatchNeg.=="B4"]; B4pos = Ypos[:,XbatchPos.=="B4"];

# visualize first 2 principal components
pNeg = scatter(B1neg[1,:],B1neg[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2neg[1,:],B2neg[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3neg[1,:],B3neg[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4neg[1,:],B4neg[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pNeg,xlabel="PC1",ylabel="PC2");

pPos = scatter(B1pos[1,:],B1pos[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2pos[1,:],B2pos[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3pos[1,:],B3pos[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4pos[1,:],B4pos[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pPos,xlabel="PC1",ylabel="PC2");

plot(pNeg, pPos, legend = :outertopright, title = ["Negative Metabolites" "Positive Metabolites"], size = (700, 400))
# -


plotattr("size")

# ### Lipids most influenced by batches

# Get variance explained
dfVarExplNeg = getVarExpl(Xneg, XbatchNeg, names(dfNeg)[6:end]);
dfVarExplPos = getVarExpl(Xpos, XbatchPos, names(dfPos)[6:end]);

first(dfVarExplNeg, 5)

first(dfVarExplPos, 5)

# +
nTop = 25# sum(dfVarExpl.VarExpl>0.1)

ticklabel = dfVarExplNeg.Lipids[1:nTop]
pNeg =bar(dfVarExplNeg.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Negative Metabolites Influenced by Batch", ylims = (0, 1))
xlabel!("Negative Metabolites")
ylabel!("Explained Variance")
# -

ticklabel = dfVarExplPos.Lipids[1:nTop]
pPos =bar(dfVarExplPos.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Positive Metabolites Influenced by Batch", ylims = (0, 1))
xlabel!("Positive Metabolites")
ylabel!("Explained Variance")

# ### Most influential batch 

adjRsquaredPerLipidsPerBatch = getVarExplPerMetaPerBatch(Xneg, XbatchNeg, dfVarExplNeg.Lipids); 

ticklabel = dfVarExplNeg.Lipids[1:nTop]
groupedbar(ticklabel, adjRsquaredPerLipidsPerBatch[1:nTop, :], 
    bar_position = :dodge, bar_width=0.7, alpha=0.5,
    xticks=(1:4:nTop, ticklabel[1:4:nTop]), 
    legend = :outertopright, label = ["B1" "B2" "B3" "B4"]
)

size(dfNeg)

# ## Correct batch effect with combat 

R"""
suppressMessages(library(sva))
fCombat <- function(myDf){
mLipids <- as.matrix(myDf[,c(-1,-2,-3,-4,-5)])

modcombat <- model.matrix(~1, data = myDf[,c(2,3,5)])

combatLipids <- ComBat(dat=t(mLipids), batch = myDf$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# modGroupFishOil <- model.matrix(~Group*FishOil, data = myDf[,c(2,3,5)])
# combatFit = lm.fit(modGroupFishOil, t(combatLipids))

return(combatLipids)
}

mLipidsBatchAdjNeg <- t(fCombat(dfNeg));
mLipidsBatchAdjPos <- t(fCombat(dfPos));

"""
@rget mLipidsBatchAdjNeg;
@rget mLipidsBatchAdjPos;

# +
# get matrix data
Xneg= copy(transpose(mLipidsBatchAdjNeg));
Xpos= copy(transpose(mLipidsBatchAdjPos));

# train a PCA model
Mneg = fit(PCA, Xneg; maxoutdim=10)
Mpos = fit(PCA, Xpos; maxoutdim=10)
# -


size(Xneg)

# get explained variance
explainedVarPCAneg = principalvars(Mneg)./tvar(Mneg)
explainedVarPCApos = principalvars(Mpos)./tvar(Mpos);

# +
ticklabel = string.(collect(1:10))
pNeg =bar(explainedVarPCAneg, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Negative Metabolites")
xlabel!("Principal Components")
ylabel!("Explained Variance");

pPos =bar(explainedVarPCApos, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Positive Metabolites")
xlabel!("Principal Components")
ylabel!("Explained Variance");

plot(pNeg, pPos, size = (800, 400))

# +
# get batch group labels
XbatchNeg = Vector(dfNeg[:,2]);
XbatchPos = Vector(dfPos[:,2]);

# apply PCA model 
Yneg = MultivariateStats.transform(Mneg, Xneg)
Ypos = MultivariateStats.transform(Mpos, Xpos)

# group results by testing set labels for color coding
B1neg = Yneg[:,XbatchNeg.=="B1"]; B1pos = Ypos[:,XbatchPos.=="B1"] 
B2neg = Yneg[:,XbatchNeg.=="B2"]; B2pos = Ypos[:,XbatchPos.=="B2"]
B3neg = Yneg[:,XbatchNeg.=="B3"]; B3pos = Ypos[:,XbatchPos.=="B3"]
B4neg = Yneg[:,XbatchNeg.=="B4"]; B4pos = Ypos[:,XbatchPos.=="B4"];

# visualize first 2 principal components
pNeg = scatter(B1neg[1,:],B1neg[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2neg[1,:],B2neg[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3neg[1,:],B3neg[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4neg[1,:],B4neg[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pNeg,xlabel="PC1",ylabel="PC2");

pPos = scatter(B1pos[1,:],B1pos[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2pos[1,:],B2pos[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3pos[1,:],B3pos[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4pos[1,:],B4pos[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pPos,xlabel="PC1",ylabel="PC2");

plot(pNeg, pPos, legend = :outertopright, title = ["Negative Lipids" "Positive Lipids"], size = (700, 400))
# -


# ### Metabolites most influenced by batches after correction

# Get variance explained
dfVarExplNeg = getVarExpl(Xneg, XbatchNeg, names(dfNeg)[6:end]);
dfVarExplPos = getVarExpl(Xpos, XbatchPos, names(dfPos)[6:end]);

first(dfVarExplNeg, 5)

first(dfVarExplPos, 5)

# +
nTop = 25# sum(dfVarExpl.VarExpl>0.1)

ticklabel = dfVarExplNeg.Lipids[1:nTop]
pNeg =bar(dfVarExplNeg.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Negative Metabolites Influenced by Batch", ylims = (0, 1))
xlabel!("Negative Metabolites")
ylabel!("Explained Variance")
# -

ticklabel = dfVarExplPos.Lipids[1:nTop]
pPos =bar(dfVarExplPos.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Positive Metabolites Influenced by Batch", ylims = (0, 1))
xlabel!("Positive Metabolites")
ylabel!("Explained Variance")

# + [markdown] kernel="SoS"
# ## Save pretreatments
# -

dfNegMeta[:, 6:end] = mLipidsBatchAdjNeg;
dfPosMeta[:, 6:end] = mLipidsBatchAdjPos;

println(string("Number of Missing for Neg: " , 
        length(findall(Matrix(dfNegMeta[:,6:end]) .=== missing))))
println(string("Number of Missing for Pos: " , 
        length(findall(Matrix(dfPosMeta[:,6:end]) .=== missing))))   

# + kernel="Julia 1.5.3"
dfNegMeta |> CSV.write("../../data/data_processed/inl2b_NegMeta.csv")
# -

dfPosMeta |> CSV.write("../../data/data_processed/inl2b_PosMeta.csv")

dfNegMeta[80:100,:]

first(dfPosMeta,3)

# +
# Join negative and positive lipids data frames

# Keep only "[ ]" and ID inside bracket:
for i in 1:size(dfNegMeta)[1]
    dfNegMeta.Sample[i]=  dfNegMeta.Sample[i][end-6:end-1]
    dfPosMeta.Sample[i]=  dfPosMeta.Sample[i][end-6:end-1]
end

# join
# difference is due to CN04 and CN05 inversion 
# dfMeta = leftjoin(dfNegMeta, dfPosMeta[:, [1;4; collect(6:end)]], on = [:Sample, :Statin]);
dfMeta = leftjoin(dfPosMeta, dfNegMeta[:, [1;4; collect(6:end)]], on = [:Sample, :Statin]);
# -

println(string("Number of Missing for Meta: " , 
        length(findall(Matrix(dfMeta[:,6:end]) .=== missing))))

dfMeta |> CSV.write("../../data/data_processed/inl2b_Meta.csv")

versioninfo()

R"""
sessionInfo()
"""
