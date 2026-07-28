################################################################################
# construct_pans_matrices.jl
#
# Purpose:
#   Construct the matrices for the PANSTEATITIS study using the processed data.
################################################################################

using CSV
using DataFrames
using DataFramesMeta
using Missings
using CategoricalArrays
using StatsBase
using Statistics
using StatsModels
using LinearAlgebra
using MatrixLM
using FreqTables
using RecipesBase

include(joinpath(@__DIR__, "..", "..", "src", "wrangling_utils.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "mLinearModel.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "recipe_plots.jl"))

function fix_covar_name(s::String)
    s = replace(s, "(" => "", ")" => "", ": " => "_")
    s = replace(s, " & " => "Ξ")
    return s
end

function make_coef_lookup(coef_names::Vector{String})
    coef_lookup = Dict{String, Int}()

    for (i, nm) in enumerate(coef_names)
        clean_nm = fix_covar_name(nm)
        coef_lookup[nm] = i
        coef_lookup[clean_nm] = i

        base_nm = split(clean_nm, "_")[1]
        coef_lookup[base_nm] = i
    end

    return coef_lookup
end

function construct_pans_matrices(;
    xCovariates = ["Status", "Sex", "Age", "Weight", "Length"],
    center_y::Bool = true
)

    ############################################################
    # Individuals
    ############################################################

    fileIndividuals = joinpath(
        @__DIR__, "..", "..", "data", "processed",
        "ST001052_ClinicalCovariates.csv"
    )

    dfInd = CSV.read(fileIndividuals, DataFrame)
    dfInd = dfInd[findall(completecases(dfInd)), :]

    ############################################################
    # Metabolite responses
    ############################################################

    fileMeta = joinpath(
        @__DIR__, "..", "..", "data", "processed",
        "nl2_Meta.csv"
    )

    dfMet = CSV.read(fileMeta, DataFrame)

    ############################################################
    # Metabolite annotations
    ############################################################

    fileRef = joinpath(
        @__DIR__, "..", "..", "data", "processed",
        "refMeta.csv"
    )

    dfRef = CSV.read(fileRef, DataFrame)

    sort!(dfRef, [:MetaboliteID])

    dfMetID = DataFrame(MetaboliteID = names(dfMet)[2:end])
    dfRef = rightjoin(dfRef, dfMetID, on = :MetaboliteID)

    dfRef.super_class = String.(replace(dfRef.super_class, missing => "NA"))
    dfRef.sub_class   = String.(replace(dfRef.sub_class, missing => "NA"))
    dfRef.Class       = String.(replace(dfRef.Class, missing => "NA"))
    dfRef.main_class  = String.(replace(dfRef.main_class, missing => "NA"))

    dfRef.sub_class = replace.(dfRef.sub_class, "O-PC" => "PC", "O-PS" => "PS")

    ############################################################
    # Standardize and order metabolite matrix
    ############################################################

    if center_y
        funStandardize!(dfMet, tocenter = true)
    end

    dfMet = dfMet[!, vcat([:SampleID], Symbol.(sort(names(dfMet)[2:end])))]

    Y = Matrix{Float64}(dfMet[:, 2:end])
    m = size(Y, 2)

    ############################################################
    # Build X matrix
    ############################################################

    vPredictorNames = copy(xCovariates)

    if "Intercept" in xCovariates
        vPredictorNames[findall(vPredictorNames .== "Intercept")] .= "1"
    end

    frml = join(vPredictorNames, " + ")

    formulaX = eval(Meta.parse(string("@formula(0 ~ ", frml, ").rhs")))

    contrasts = Dict(
        :Status => EffectsCoding(base = sort(unique(dfInd.Status))[2]),
        :Sex    => EffectsCoding(base = sort(unique(dfInd.Sex))[1]),
    )

    X = modelmatrix(formulaX, dfInd, hints = contrasts)

    if frml == "1"
        coef_names = ["(Intercept)"]
    else
        sch = schema(formulaX, dfInd, contrasts)
        coef_names = coefnames(apply_schema(formulaX, sch))
    end

    pseudo_coef_names = fix_covar_name.(coef_names)
    coef_lookup = make_coef_lookup(coef_names)

    ############################################################
    # Z matrix
    ############################################################

    Z = Matrix{Float64}(I, m, m)

    ############################################################
    # Superclass and subclass grouping
    ############################################################

    metabolite_ids = names(dfMet)[2:end]

    ref_super = Dict(dfRef.MetaboliteID .=> dfRef.super_class)
    ref_sub   = Dict(dfRef.MetaboliteID .=> dfRef.sub_class)

    vSuperClass = [get(ref_super, id, missing) for id in metabolite_ids]
    vSubClass   = [get(ref_sub, id, missing) for id in metabolite_ids]

    vSuperClass = categorical(vSuperClass)
    vSubClass   = categorical(vSubClass)

    subclass_levels = unique(vSubClass)
    H = length(subclass_levels)

    subclass_index = Dict(subclass_levels[i] => i for i in 1:H)
    subclass_of_met = [subclass_index[x] for x in vSubClass]

    superclass_levels = unique(vSuperClass)
    G = length(superclass_levels)

    superclass_index = Dict(superclass_levels[i] => i for i in 1:G)
    superclass_of_met = [superclass_index[x] for x in vSuperClass]

    superclass_of_sub = Vector{Int}(undef, H)

    for h in 1:H
        J = findall(==(h), subclass_of_met)
        @assert !isempty(J) "Subclass $h has no metabolites."

        gvals = unique(superclass_of_met[J])
        @assert length(gvals) == 1 "Subclass $h appears in multiple superclasses: $gvals"

        superclass_of_sub[h] = gvals[1]
    end

    @assert maximum(subclass_of_met) == H
    @assert maximum(superclass_of_sub) == G
    @assert size(X, 1) == size(Y, 1) "X and Y must have the same number of subjects"
    @assert size(X, 2) == length(coef_names) "Coefficient names must align with X columns"
    @assert size(Z) == (m, m) "Z must be an m × m identity matrix"
    @assert length(subclass_of_met) == m "Each metabolite must have a subclass assignment"

    return (
        dfInd = dfInd,
        dfMet = dfMet,
        dfRef = dfRef,
        X = Matrix{Float64}(X),
        Y = Y,
        Z = Z,
        formulaX = formulaX,
        contrasts = contrasts,
        coef_names = coef_names,
        pseudo_coef_names = pseudo_coef_names,
        coef_lookup = coef_lookup,
        idx_covar = collect(1:length(coef_names)),
        subclass_of_met = Int.(subclass_of_met),
        superclass_of_met = Int.(superclass_of_met),
        superclass_of_sub = Int.(superclass_of_sub),
        H = H,
        G = G,
        subclass_levels = subclass_levels,
        superclass_levels = superclass_levels,
        metabolite_ids = metabolite_ids
    )
end
