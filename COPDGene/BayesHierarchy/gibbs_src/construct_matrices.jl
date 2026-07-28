############################################################
# construct_matrices.jl
#
# Purpose:
#   Construct the COPDGene analysis objects needed downstream:
#       X matrix
#       Y matrix
#       Z matrix
#       coefficient names
#       covariate indices
#
############################################################

using CSV
using DataFrames
using CategoricalArrays
using StatsBase
using Statistics
using StatsModels
using LinearAlgebra
using MatrixLM
using RecipesBase

############################################################
# COPDGene data-wrangling functions
########################################################

include(joinpath(@__DIR__, "..", "..", "src", "wrangle_utils.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "utils_copd_spiro.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "demog.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "mLinearModel.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "myPlots.jl"))

############################################################
# Helper function
############################################################

function fix_covar_name(s::String)
    return replace(
        replace(
            replace(s, "(" => "", ")" => ""),
            ": " => "_"
        ),
        " & " => "Ξ"
    )
end

############################################################
# Construct matrices for COPDGene analysis
############################################################

function construct_copdgene_matrices(;
    xCovariates = [
        "Sex",
        "Age",
        "BMI",
        "SmokingPackYears",
        "PercentEmphysema",
        "COPD",
        "NHW",
        "CurrentSmoker"
    ],
    site_adjusted::Bool = true,
    center_y::Bool = true
)

    ########################################################
    # Load COPDGene data
    ########################################################

    copd = get_data("COPDGene")

    ########################################################
    # Build formula for X
    ########################################################

    vPredictorNames = copy(xCovariates)

    for i in eachindex(vPredictorNames)
        if vPredictorNames[i] == "Intercept"
            vPredictorNames[i] = "1"
        end
    end

    frml = join(vPredictorNames, " + ")

    if site_adjusted
        frml = frml * " + Site"
    end

    formulaX = eval(Meta.parse(string("@formula(0 ~ ", frml, ").rhs")))

    ########################################################
    # Contrasts for categorical variables
    ########################################################

    contrasts_copd = Dict(
        :Sex           => EffectsCoding(base = sort(unique(copd.dfInd.Sex))[2]),
        :NHW           => EffectsCoding(base = sort(unique(copd.dfInd.NHW))[1]),
        :Site          => EffectsCoding(base = sort(unique(copd.dfInd.Site))[1]),
        :CurrentSmoker => EffectsCoding(base = sort(unique(copd.dfInd.CurrentSmoker))[1]),
        :COPD          => EffectsCoding(base = sort(unique(copd.dfInd.COPD))[1]),
    )

    ########################################################
    # Construct X
    ########################################################

    X = modelmatrix(formulaX, copd.dfInd; hints = contrasts_copd)

    ########################################################
    # Coefficient names
    ########################################################

    sch_copd = schema(formulaX, copd.dfInd, contrasts_copd)
    coef_names = coefnames(apply_schema(formulaX, sch_copd))
    pseudo_coef_names = fix_covar_name.(coef_names)

    ########################################################
    # Indices of non-site covariates
    ########################################################

    idx_covar = findall(.!occursin.("Site", coef_names))

    ########################################################
    # Construct Z
    # For this model, Z is the identity matrix
    ########################################################

    Z = Matrix{Float64}(I, size(copd.mY, 2), size(copd.mY, 2))

    ########################################################
    # Construct Y
    ########################################################

    Y = center_y ? copd.mY .- mean(copd.mY; dims = 1) : copd.mY

    ########################################################
    # Metabolite superclass and subclass grouping
    ########################################################

    copd.dfRef.SuperClassID = categorical(copd.dfRef.SuperClassID)
    copd.dfRef.SubClassID   = categorical(copd.dfRef.SubClassID)

    superclass_of_met = levelcode.(copd.dfRef.SuperClassID)
    subclass_of_met   = levelcode.(copd.dfRef.SubClassID)

    G = length(levels(copd.dfRef.SuperClassID))
    H = length(levels(copd.dfRef.SubClassID))

    super_of_sub = similar(1:H)

    for (h, sublev) in enumerate(levels(copd.dfRef.SubClassID))
        idx = findfirst(==(sublev), copd.dfRef.SubClassID)
        super_of_sub[h] = superclass_of_met[idx]
    end

    ########################################################
    # Validate matrix and annotation alignment
    ########################################################

    @assert size(X, 1) == size(Y, 1) "X and Y must have the same number of subjects"
    @assert size(Z, 1) == size(Y, 2) "Z rows must match the number of metabolites"
    @assert length(subclass_of_met) == size(Y, 2) "Subclass annotations must align with Y"
    @assert length(coef_names) == size(X, 2) "Coefficient names must align with X columns"

    # Create a lookup dictionary for covariates
    coef_lookup = Dict{String,Int}()
    for (i, nm) in enumerate(coef_names)
        base = split(nm, ":")[1]
        coef_lookup[base] = i
    end

    ########################################################
    # Return all objects needed downstream
    ########################################################

    return (
        copd = copd,
        X = X,
        Y = Y,
        Z = Z,
        formulaX = formulaX,
        contrasts = contrasts_copd,
        coef_names = coef_names,
        pseudo_coef_names = pseudo_coef_names,
        idx_covar = idx_covar,
        xCovariates = xCovariates,
        coef_lookup = coef_lookup,
        superclass_of_met = superclass_of_met,
        subclass_of_met = subclass_of_met,
        super_of_sub = super_of_sub,
        G = G,
        H = H
    )

end
