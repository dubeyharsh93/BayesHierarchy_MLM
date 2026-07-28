######################################################################
# construct_sams_matrices.jl
#
# Purpose:
#   Construct the matrices for the SAMS study using the processed data.
######################################################################

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using LinearAlgebra
using MatrixLM
using FreqTables
using StatsBase
using RecipesBase

include(joinpath(@__DIR__, "..", "..", "src", "mLinearModel.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "myPlots.jl"))

function construct_sams_matrices(;
    slctFishOil::String = "true",
    rdZ::String = "Triglycerides",
    rdX::String = "false",
    isunpairFlag::Bool = true,
)

    # Load data
    fileLipids = realpath(joinpath(@__DIR__,"..", "..", "data", "data_processed", "inl2b_Lipids.csv"))
    dfLipids = DataFrame(CSV.File(fileLipids))

    df = getCases(dfLipids, isunpaired = isunpairFlag)

    hasFishOil = eval(Meta.parse(slctFishOil))

    if !hasFishOil
        filter!(row -> row.FishOil == "no", df)
    end

    funStandardize!(df, isunpaired = isunpairFlag)

    # Select Z/raw lipid reference file
    fZraw = Dict(
        "All" => joinpath(@__DIR__, "..", "..", "data", "data_processed", "ZmatRawAll.csv"),
        "Triglycerides" => joinpath(@__DIR__, "..", "..", "data", "data_processed", "ZmatRawTG.csv"),
        "Phospholipids" => joinpath(@__DIR__, "..", "..", "data", "data_processed", "ZmatRawPhos.csv"),
        "Phospholipids: PC" => joinpath(@__DIR__, "..", "..", "data", "data_processed", "ZmatRawPhosPC.csv"),
    )

    dfZraw = DataFrame(CSV.File(fZraw[rdZ]))
    dfZraw = dfZraw[findall(x -> [x] ⊆ names(df), dfZraw.lipID), :]

    m = length(dfZraw.lipID)
    Z = Matrix{Float64}(I, m, m)

    is4ways = eval(Meta.parse(rdX))

    # Construct X and Y
    if hasFishOil
        if is4ways
            X, Y = getXY4ways(
                df;
                responseSelection = dfZraw.lipID,
                isunpaired = isunpairFlag
            )

            coef_names = [
                "CN-No Fish Oil",
                "CN-Fish Oil",
                "CS-No Fish Oil",
                "CS-Fish Oil"
            ]
        else
            X, Y = getXY(
                df;
                responseSelection = dfZraw.lipID,
                isunpaired = isunpairFlag
            )

            coef_names = [
                "Intercept",
                "SAMS status",
                "Fish Oil",
                "Interaction SAMS-Fish Oil"
            ]
        end
    else
        X, Y = getXYnoFishOil(
            df;
            responseSelection = dfZraw.lipID,
            isunpaired = isunpairFlag
        )

        coef_names = [
            "Intercept",
            "SAMS status"
        ]
    end

    X = Matrix{Float64}(X)
    Y = Matrix{Float64}(Y)

    # The manuscript analysis excludes the intercept.
    X = X[:, 2:end]
    coef_names = coef_names[2:end]

    # Total DB grouping
    db_vals = dfZraw.Total_DB

    cat_DB_idx = map(db_vals) do x
        if x < 3
            1
        elseif x < 6
            2
        elseif x < 9
            3
        else
            4
        end
    end

    subclass_of_met = Int.(cat_DB_idx)
    H = length(unique(subclass_of_met))

    group_labels = ["DB <3", "3–6", "6–9", "≥9"]

    coef_lookup = Dict{String, Int}()
    for (i, nm) in enumerate(coef_names)
        coef_lookup[nm] = i
    end

    @assert size(X, 1) == size(Y, 1) "X and Y must have the same number of subjects"
    @assert size(Z, 1) == size(Y, 2) "Z rows must match the number of lipids"
    @assert length(subclass_of_met) == size(Y, 2) "Subclass annotations must align with Y"
    @assert length(coef_names) == size(X, 2) "Coefficient names must align with X columns"

    return (
        df = df,
        dfZraw = dfZraw,
        X = X,
        Y = Y,
        Z = Z,
        coef_names = coef_names,
        coef_lookup = coef_lookup,
        pseudo_coef_names = coef_names,
        idx_covar = collect(1:length(coef_names)),
        subclass_of_met = subclass_of_met,
        H = H,
        group_labels = group_labels,
        lipid_names = dfZraw.lipID,
        hasFishOil = hasFishOil,
        is4ways = is4ways,
        rdZ = rdZ
    )
end
