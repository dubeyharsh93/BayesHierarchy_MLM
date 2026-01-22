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

# This notebook carrieds out the preprocessing steps for the lipidomics data.

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

# ### Ext. Functions

include(joinpath(@__DIR__,"..","..","src","preprocessing.jl" ));
include(joinpath(@__DIR__,"..","..","src","utils.jl" ));

# ### Load data

negLipids = realpath(joinpath(@__DIR__,"..","..","data","data_primary","Elam_NEG_LIPIDS_all-variables.csv" ))
dfNegLipids = DataFrame(CSV.File(negLipids));

posLipids = realpath(joinpath(@__DIR__,"..","..","data","data_primary","Elam_POS_LIPIDS-all-variables.csv" ))
dfPosLipids = DataFrame(CSV.File(posLipids));

# ## Chain of preprocessing

# +
# diGraph(["Wrangle", "Imputation", "Normalization", "Log2 Transformation", "Batch Effect Correction"])
# -

# Wrangle => Imputation => Normalization => Log2 Transformation => Batch Effect Correction

# ## Wrangle data
# ----

# ### Lipids reference dictionnary

# * We create a reference dictionnary for the lipids names.

# +
# isCossRef = true
# Generate data frames based on cross reference naming for positive and negative lipids
dfNegLipids, dfNegCrossRef = crossRefLipids(dfNegLipids, "negLip");
dfPosLipids, dfPosCrossRef = crossRefLipids(dfPosLipids, "posLip", startCol = 6);

# save cross reference
dfNegCrossRef |> CSV.write("../../data/data_processed/inl2b_NegLipids_Xref.csv");
dfPosCrossRef |> CSV.write("../../data/data_processed/inl2b_PosLipids_Xref.csv");

# merge neg and pos look up table
dfLipCrossRef = deepcopy(dfNegCrossRef);
append!(dfLipCrossRef, dfPosCrossRef)
dfLipCrossRef |> CSV.write("../../data/data_processed/inl2b_Lipids_Xref.csv");
# -

# * Display the first row of negative lipids data frame:

first(dfNegLipids)

# * Display the first row of positive lipids data frame:

first(dfPosLipids)

# ### Create a Group variable: CN, CS, CSbaseline

# * Create a categorical variable named `Group` that indicates if a patient belonged to the control group (CN), the symptomatic group not on statin (CSbaseline) or the rechallenge symptomatic group on statin (CS).

dfNegLipids = addGroupCat(copy(dfNegLipids), false);
dfPosLipids = addGroupCat(copy(dfPosLipids), false);

# Check group population:

countmap(dfNegLipids.Group; alg = :dict)

# ### Convert to categorical: Batch, Statin and Oil

dfNegLipids = catBGSF(dfNegLipids, false);
dfPosLipids = catBGSF(dfPosLipids, false);

# ### Check tables

freqtable(dfNegLipids, :Group, :Batch)

freqtable(dfNegLipids, :Group, :FishOil)

# ## Impute missing data  
# ___

# HM (Half of the Minimum): This method replaces missing elements with half of the minimum of non-missing elements in the corresponding variable.

# QRILC ( quantile regression imputation of left-censored): This method imputes the left-censored data by randomly drawing values from a truncated normal distribution.

# +
# dfNegLipids = imputeHM(dfNegLipids, startCol = 6);
# dfPosLipids = imputeHM(dfPosLipids, startCol = 6);

dfNegLipids = imputeQRILC(dfNegLipids; startCol = 6, logTX = true);
dfPosLipids = imputeQRILC(dfPosLipids; startCol = 6, logTX = true);
# -

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

dfNegLipids = pqnorm(dfNegLipids, startCol = 6);
dfPosLipids = pqnorm(dfPosLipids, startCol = 6);

# ## Transformation
# ---
#
# A simple and widely used transformation to make data more symmetric and homoscedastic is the log-transformation.

dfNegLipids = log2tx(dfNegLipids, startCol = 6);
dfPosLipids = log2tx(dfPosLipids, startCol = 6);

first(dfPosLipids)

# ## Adjusting for batch effects
# ---

# ### Check for batch effects

R"""
suppressMessages(library(mixOmics))
suppressMessages(library(tidyverse));
"""

dfNeg = catBGSF(dfNegLipids);
dfPos = catBGSF(dfPosLipids);

@rput dfNeg;
@rput dfPos;

R"summary(dfNeg[,c(1:5)])"

# +
# get matrix data
Xneg= copy(transpose(Matrix(dfNeg[:,6:end])));
Xpos= copy(transpose(Matrix(dfPos[:,6:end])));

# train a PCA model
Mneg = fit(PCA, Xneg; maxoutdim=10)
Mpos = fit(PCA, Xpos; maxoutdim=10);
Mpos;
# -


size(Xneg)

# get explained variance
explainedVarPCAneg = principalvars(Mneg)./tvar(Mneg)
explainedVarPCApos = principalvars(Mpos)./tvar(Mpos);

# +
ticklabel = string.(collect(1:10))
pNegPCA =bar(explainedVarPCAneg, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Negative Lipids")
xlabel!("Principal Components")
ylabel!("Explained Variance");

pPosPCA =bar(explainedVarPCApos, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Positive Lipids")
xlabel!("Principal Components")
ylabel!("Explained Variance");

plot(pNegPCA, pPosPCA, size = (800, 400))

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
pNegScat = scatter(B1neg[1,:],B1neg[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2neg[1,:],B2neg[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3neg[1,:],B3neg[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4neg[1,:],B4neg[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pNegScat,xlabel="PC1",ylabel="PC2");

pPosScat = scatter(B1pos[1,:],B1pos[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2pos[1,:],B2pos[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3pos[1,:],B3pos[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4pos[1,:],B4pos[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pPosScat,xlabel="PC1",ylabel="PC2");

plot(pNegScat, pPosScat, legend = :outertopright, title = ["Negative Lipids" "Positive Lipids"], size = (700, 400))
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
    yflip=false, legend = false, title = "Top Negative Lipids Influenced by Batch", ylims = (0, 1))
xlabel!("Negative Lipids")
ylabel!("Explained Variance")
# -

ticklabel = dfVarExplPos.Lipids[1:nTop]
pPos =bar(dfVarExplPos.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Positive Lipids Influenced by Batch", ylims = (0, 1))
xlabel!("Positive Lipids")
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
Mpos = fit(PCA, Xpos; maxoutdim=10);
# -


size(Xneg)

# get explained variance
explainedVarPCAneg = principalvars(Mneg)./tvar(Mneg)
explainedVarPCApos = principalvars(Mpos)./tvar(Mpos);

# +
ticklabel = string.(collect(1:10))
pNegPCAAdj =bar(explainedVarPCAneg, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Negative Lipids")
xlabel!("Principal Components")
ylabel!("Explained Variance");

pPosPCAAdj =bar(explainedVarPCApos, orientation=:v, xticks=(1:10, ticklabel),
    yflip=false, legend = false, title = "Positive Lipids")
xlabel!("Principal Components")
ylabel!("Explained Variance");

plot(pNegPCAAdj, pPosPCAAdj, size = (800, 400))

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
pNegScatAdj = scatter(B1neg[1,:],B1neg[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2neg[1,:],B2neg[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3neg[1,:],B3neg[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4neg[1,:],B4neg[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pNegScatAdj,xlabel="PC1",ylabel="PC2");

pPosScatAdj = scatter(B1pos[1,:],B1pos[2,:], marker=:auto, markersize=4, linewidth=0, label = "B1")
scatter!(B2pos[1,:],B2pos[2,:], marker=:utriangle,linewidth=0, label = "B2")
scatter!(B3pos[1,:],B3pos[2,:], marker=:+,linewidth=0, label = "B3")
scatter!(B4pos[1,:],B4pos[2,:], marker=:x,linewidth=0, label = "B4")
plot!(pPosScatAdj,xlabel="PC1",ylabel="PC2");

plot(pNegScatAdj, pPosScatAdj, legend = :outertopright, title = ["Negative Lipids" "Positive Lipids"], size = (700, 400))
# -


# ### Lipids most influenced by batches after correction

# Get variance explained
dfVarExplNeg = getVarExpl(Xneg, XbatchNeg, names(dfNeg)[6:end]);
dfVarExplPos = getVarExpl(Xpos, XbatchPos, names(dfPos)[6:end]);

first(dfVarExplNeg, 5)

first(dfVarExplPos, 5)

# +
nTop = 25# sum(dfVarExpl.VarExpl>0.1)

ticklabel = dfVarExplNeg.Lipids[1:nTop]
pNegAdj =bar(dfVarExplNeg.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Negative Lipids Influenced by Batch", ylims = (0, 1))
xlabel!("Negative Lipids")
ylabel!("Explained Variance")
# -

ticklabel = dfVarExplPos.Lipids[1:nTop]
pPosAdj =bar(dfVarExplPos.VarExpl[1:nTop], orientation=:v, xticks=(1:4:nTop, ticklabel[1:4:nTop]),
    yflip=false, legend = false, title = "Top Positive Lipids Influenced by Batch", ylims = (0, 1))
xlabel!("Positive Lipids")
ylabel!("Explained Variance")

plot(pNeg, pNegAdj,
    legend = :false, 
    grid = false, 
    title = ["Not Corrected" "Batch Corrected"],
    ylabel = ["Explained Variance" ""],
    xaxis = false,
    left_margin = (10, :mm),
    size = (1200, 400)
)

plot(pPos, pPosAdj,
    legend = :false, 
    grid = false, 
    title = ["Not Corrected" "Batch Corrected"],
    ylabel = ["Explained Variance" ""],
    xaxis = false,
    left_margin = (10, :mm),
    size = (1200, 400)
)

# +
# Xneg, XbatchNeg, names(dfNeg)[6:end]
# coeftable(lm_fit) |> c -> c.cols[c.pvalcol][c.rownms .== "x"]
n = size(Xneg)[1]

adjPval = zeros(n)

for i in 1:n
    dftest = DataFrame(X = CategoricalArray(string.(XbatchNeg)), Y = Xneg[i,:])
    out = lm(@formula(Y ~ X), dftest);
    out0 = lm(@formula(Y ~ 1), dftest);
    my_ftest=ftest(out0.model, out.model)
    adjPval[i] = my_ftest.pval[2]
end

R"""
suppressMessages(library(stats))
suppressMessages(library(qvalue));
"""
@rput adjPval;

R"""
qobj <- qvalue(p = adjPval)
qVals <- qobj$qvalues;
"""
@rget qVals;

describe(qVals)

# +
n = size(Xpos)[1]

adjPval = zeros(n)

for i in 1:n
    dftest = DataFrame(X = CategoricalArray(string.(XbatchPos)), Y = Xpos[i,:])
    out = lm(@formula(Y ~ X), dftest);
    out0 = lm(@formula(Y ~ 1), dftest);
    my_ftest=ftest(out0.model, out.model)
    adjPval[i] = my_ftest.pval[2]
end

R"""
suppressMessages(library(stats))
suppressMessages(library(qvalue));
"""
@rput adjPval;

R"""
qobj <- qvalue(p = adjPval)
qVals <- qobj$qvalues;
"""
@rget qVals;

describe(qVals)

# + [markdown] kernel="SoS"
# ## Save pretreatments
# -

dfNegLipids[:, 6:end] = mLipidsBatchAdjNeg;
dfPosLipids[:, 6:end] = mLipidsBatchAdjPos;

# + kernel="Julia 1.5.3"
dfNegLipids |> CSV.write("../../data/data_processed/inl2b_NegLipids.csv")
# -

dfPosLipids |> CSV.write("../../data/data_processed/inl2b_PosLipids.csv")

# Join negative and positive lipids data frames
dfLipids = leftjoin(dfNegLipids, dfPosLipids[:, [1; collect(6:end)]], on = :Sample);

dfLipids |> CSV.write("../../data/data_processed/inl2b_Lipids.csv")

versioninfo()

R"""
sessionInfo()
"""


