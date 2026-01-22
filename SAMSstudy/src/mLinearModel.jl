# -*- coding: utf-8 -*-
#############
# LIBRARIES #
#############
# using CSV, DataFrames, DataFramesMeta, Missings #, CategoricalArrays
# using StatsBase, Statistics, MatrixLM, Random, Distributions, StatsModels#,  MultivariateStats,
# using LinearAlgebra
# using FreqTables, Plots


# #################
# FUNCTIONS LIST #
# #################

#=
Functions list 

- getCases: 
    Returns a data frame with only paired (rechallenge cases) or unpaired cases, default `isunpaired` is true.

- funStandardize: 
    Standardize data, but no centering in the paired cases.

- getXY: 
    Returns X covariates and Y response matrix.

- getXYnoFishOil:
    Returns X covariates and Y response matrix, for individuals without fish oil.

- getXY4ways:
    Returns X covariates and Y response matrix, where X covariates contains 4 categories
    {CN-Fish Oil, CN-No Fish Oil, CS-Fish Oil, CS-No Fish Oil}.

- getCoefs:
    Returns coefficients, confidence of interval, T-stats.

- getPermPvals:
    Returns T-statistics and corresponding P-value based on permutations tests.

=#


# ############
# FUNCTIONS #
# ############

"""
**getCases** -*Function*.

    getCases(df::DataFrame, isunpaired::Bool = false) => DataFrame

Returns a data frame with only paired (rechallenge cases) or unpaired cases, default `isunpaired` is true. 
"""
function getCases(dfInput::DataFrame; isunpaired::Bool = true)
    
    if isunpaired
        # Keep CN and CS group
        dfOutput = filter(row -> row.Group .!= "CSbaseline", dfInput);
#         saveFileNameRslts = "../data/dataprocessed/rslts_CNvsCS_Lipids_Z.csv";
    else    
        # Rechallenge: Keep CSbaseline and CS group
        filter!(row -> row.Group .!= "CN", dfInput);
#         saveFileNameRslts = "../data/dataprocessed/rslts_CSbaselinevsCS_Lipids_Z.csv";
        
        # compare baseline and rechallenge cases
        dfSmpl1 = @linq filter(row -> row.Group .== "CSbaseline", dfInput) |> select!(Not(:Batch))
        dfSmpl2 = @linq filter(row -> row.Group .== "CS", dfInput) |> select!(Not(:Batch));

        # Keep only "[ ]" and ID inside bracket:
        for i in 1:size(dfSmpl1)[1]
            dfSmpl1.Sample[i]=  dfSmpl1.Sample[i][end-5:end]
            dfSmpl2.Sample[i]=  dfSmpl2.Sample[i][end-5:end]
        end

        # Check if both samples ID are identical:
        dfDiffSmpl = append!(antijoin(dfSmpl1, dfSmpl2, on = :Sample),antijoin(dfSmpl2, dfSmpl1, on = :Sample)) 

        # Remove not common IDs
        for i in 1:size(dfDiffSmpl)[1]
            if dfDiffSmpl.Group[i] == "CSbaseline"
                dfSmpl1 = @linq filter!(row -> row.Sample .!= dfDiffSmpl.Sample[i], dfSmpl1) |> DataFrames.sort(:Sample)
            else
                dfSmpl2 = @linq filter!(row -> row.Sample .!= dfDiffSmpl.Sample[i], dfSmpl2) |> DataFrames.sort(:Sample);
            end    
        end

        # Check Fish Oil is identical
        if dfSmpl1.FishOil != dfSmpl2.FishOil
            @warn "Fish Oil variable different into the paired sample"
        end

#         dfSmpl1 |> CSV.write("../../data/dataprocessed/dfBaseline.csv")
#         dfSmpl2 |> CSV.write("../../data/dataprocessed/dfCS.csv")


        # Get difference:
        mCSbaseline = Matrix(dfSmpl1[:, 5:end]);
        mCS = Matrix(dfSmpl2[:, 5:end]);
        dfPaired = deepcopy(dfSmpl1);

        # compute difference
        dfPaired[:, 5:end] = mCS-mCSbaseline;
        select!(dfPaired, Not([:Group,:Statin])); # remove Group and Statin column since we get the difference
        dfOutput = deepcopy(dfPaired);
            
    end
    
    return dfOutput
end 








function getCases_new(dfInput::DataFrame; isunpaired::Bool = true)
    
    if isunpaired
        # Keep CN and CS group
        dfOutput = filter(row -> row.Group .!= "CSbaseline", dfInput);
#         saveFileNameRslts = "../data/dataprocessed/rslts_CNvsCS_Lipids_Z.csv";
    else    
        # Rechallenge: Keep CSbaseline and CS group
        filter!(row -> row.Group .!= "CN", dfInput);
#         saveFileNameRslts = "../data/dataprocessed/rslts_CSbaselinevsCS_Lipids_Z.csv";
        
        # compare baseline and rechallenge cases
        dfSmpl1 = @linq filter(row -> row.Group .== "CSbaseline", dfInput) |> select!(Not(:Batch))
        dfSmpl2 = @linq filter(row -> row.Group .== "CS", dfInput) |> select!(Not(:Batch));

        # Keep only "[ ]" and ID inside bracket:
        for i in 1:size(dfSmpl1)[1]
            dfSmpl1.Sample[i]=  dfSmpl1.Sample[i][end-5:end]
            dfSmpl2.Sample[i]=  dfSmpl2.Sample[i][end-5:end]
        end

        # Check if both samples ID are identical:
        dfDiffSmpl = append!(antijoin(dfSmpl1, dfSmpl2, on = :Sample),antijoin(dfSmpl2, dfSmpl1, on = :Sample)) 

        # Remove not common IDs
        for i in 1:size(dfDiffSmpl)[1]
            if dfDiffSmpl.Group[i] == "CSbaseline"
                dfSmpl1 = @linq filter!(row -> row.Sample .!= dfDiffSmpl.Sample[i], dfSmpl1) |> DataFrames.sort(:Sample)
            else
                dfSmpl2 = @linq filter!(row -> row.Sample .!= dfDiffSmpl.Sample[i], dfSmpl2) |> DataFrames.sort(:Sample);
            end    
        end

        # Check Fish Oil is identical
        if dfSmpl1.FishOil != dfSmpl2.FishOil
            @warn "Fish Oil variable different into the paired sample"
        end

#         dfSmpl1 |> CSV.write("../../data/dataprocessed/dfBaseline.csv")
#         dfSmpl2 |> CSV.write("../../data/dataprocessed/dfCS.csv")


        # Get difference:
        mCSbaseline = Matrix(dfSmpl1[:, 5:end]);
        mCS = Matrix(dfSmpl2[:, 5:end]);
        dfPaired = deepcopy(dfSmpl1);

        # compute difference
        dfPaired[:, 5:end] = mCS-mCSbaseline;
        select!(dfPaired, Not([:Group,:Statin])); # remove Group and Statin column since we get the difference
        dfOutput = deepcopy(dfPaired);
            
    end
    
    lipidnames = names(dfOutput)
    idxPredictors = findall(
                        .!((occursin.("neg",lipidnames)) .|| (occursin.("pos",lipidnames))) 
                        )
    
    
    return select(dfOutput, idxPredictors), select(dfOutput, Not(idxPredictors))
end 


"""
**funStandardize** -*Function*.

    funStandardize(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, but no centering in the paired cases. 
"""
function funStandardize(dfInput::DataFrame; unpaired::Bool = true)
    dfOutput = deepcopy(dfInput);
    if unpaired
        Xlip = Matrix(dfOutput[:,6:end]);
        dt  = fit(ZScoreTransform, Xlip, dims=1)
        mLipids = StatsBase.transform(dt, Xlip)
        dfOutput[:,6:end] = mLipids
    else 
        # No centering for paired model
        Xlip = Matrix(dfOutput[:,3:end]);
        dt  = fit(ZScoreTransform, Xlip, dims=1, center = false)
        mLipids = StatsBase.transform(dt, Xlip)
        dfOutput[:,3:end] = mLipids
    end
    
    return dfOutput
end

"""
**funStandardize!** -*Function*.

    funStandardize!(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, but no centering in the paired cases. 
"""
function funStandardize!(dfInput::DataFrame; isunpaired::Bool = true)
    if isunpaired
        Xlip = Matrix(dfInput[:,6:end]);
        dt  = fit(ZScoreTransform, Xlip, dims=1)
        mLipids = StatsBase.transform(dt, Xlip)
        dfInput[:,6:end] = mLipids
    else 
        # No centering for paired model
        Xlip = Matrix(dfInput[:,3:end]);
        dt  = fit(ZScoreTransform, Xlip, dims=1, center = false)
        mLipids = StatsBase.transform(dt, Xlip)
        dfInput[:,3:end] = mLipids
    end
    
    return dfInput
end



"""
**funStandardize** -*Function*.

    funStandardize(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, but no centering in the paired cases. 
"""
function funStandardize_new(dfInput::DataFrame; unpaired::Bool = true)
    dfOutput = deepcopy(dfInput);
    if unpaired
        Xlip = Matrix(dfOutput);
        dt  = fit(ZScoreTransform, Xlip, dims=1)
        mLipids = StatsBase.transform(dt, Xlip)
        dfOutput[:,:] = mLipids
    else 
        # No centering for paired model
        Xlip = Matrix(dfOutput[:,:]);
        dt  = fit(ZScoreTransform, Xlip, dims=1, center = false)
        mLipids = StatsBase.transform(dt, Xlip)
        dfOutput[:,:] = mLipids
    end
    
    return dfOutput
end

"""
**funStandardize!** -*Function*.

    funStandardize!(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, but no centering in the paired cases. 
"""
function funStandardize_new!(dfInput::DataFrame; isunpaired::Bool = true)
    if isunpaired
        Xlip = Matrix(dfInput);
        dt  = fit(ZScoreTransform, Xlip, dims=1)
        mLipids = StatsBase.transform(dt, Xlip)
        dfInput[:,:] = mLipids
    else 
        # No centering for paired model
        Xlip = Matrix(dfInput);
        dt  = fit(ZScoreTransform, Xlip, dims=1, center = false)
        mLipids = StatsBase.transform(dt, Xlip)
        dfInput[:,:] = mLipids
    end
    
    return dfInput
end


"""
**getCoefs** -*Function*.

getCoefsgetCoefs(
        Y::AbstractMatrix, 
        X::AbstractMatrix, 
        Zin::AbstractMatrix; 
        zcritical = 1.96,
        hasXIntercept = true, hasZIntercept = true) => Array{Float,2}, Array{Float,2}, Array{Float,2}, Array{Float,2}

Returns coefficients, confidence of interval, T-stats and variance.

"""
function getCoefs(Y::AbstractMatrix, X::AbstractMatrix, Zin::AbstractMatrix; 
    zcritical = 1.96, hasXIntercept = true, hasZIntercept = true)
       
    
    #########
    # MODEL #
    #########
    
    # Construct a RawData object
    dat = RawData(Response(Y), Predictors(X,Zin, hasXIntercept, hasZIntercept));

    # Estimate coefficients
    # est = mlm(dat, hasXIntercept = false, hasZIntercept = false); # v. 0.1.3
    est = mlm(dat, addXIntercept = hasXIntercept, addZIntercept = hasZIntercept);
   
    # get estimate
    estCoefOut = MatrixLM.coef(est);

    # get variance estimate
    σ² = est.varB
    # get SE
    SE = sqrt.(σ²)
    CIOut = SE.*zcritical
    
    # get T-statistics
    tStatsOut = t_stat(est, true)
    
    return estCoefOut, CIOut, tStatsOut, σ² 
end

"""
**getPermPvals** -*Function*.

getPermPvals(Y::AbstractMatrix, X::AbstractMatrix, Zin::AbstractMatrix; 
                nPerms = 10) => Array{Float,2}, Array{Float,2}

Returns T-statistics and corresponding P-value based on permutations tests.
"""
function getPermPvals(Y::AbstractMatrix, X::AbstractMatrix, Zin::AbstractMatrix; nPerms = 10)    
    #########
    # MODEL #
    #########
    
    # Construct a RawData object
    dat = RawData(Response(Y), Predictors(X,Zin));

    # Estimate t-stats and pVals
    tStatsOut, pValsOut = mlm_perms(dat, nPerms, hasXIntercept=false, hasZIntercept=false);
    
    return tStatsOut, pValsOut 
end

##############################################################################################













"""
**getXY** -*Function*.

    getXY(dfInput::DataFrame; 
                    responseSelection = [""],
                    isunpaired::Bool = true) => DataFrame, DataFrame

Returns X covariates and Y response matrix.
"""
function getXY(dfInput::DataFrame; 
                    responseSelection = [""],
                    isunpaired::Bool = true)
    
    responseSelection = intersect(names(dfInput), responseSelection) 
    # Select responses
    if isunpaired
        colSelection = vcat(["Sample", "Batch", "Group", "Statin", "FishOil"], responseSelection)
    else
        colSelection = vcat(["Sample", "FishOil"], responseSelection)        
    end;
    dfInput= dfInput[:, colSelection]
    
    # get dimensions of sample and covariates
    n = size(dfInput)[1]
    p = ifelse(isunpaired, 4, 2); # Intercept | Group | Fish Oil | Group:Fish Oil 
    
    #####
    # X #
    ##### 
    # Generate X matrix
    X = zeros(n,p)
    # Intercept
    X[:,1] = ones(n,1)
    if isunpaired
        # Group
        X[:,2] = (dfInput.Group .== "CS").*2 .- 1
        # Fish Oil
        X[:,3] = (dfInput.FishOil .== "yes").*2 .- 1
        # Interaction Group:Fish Oil
        X[:, 4] = X[:, 2] .* X[:, 3]
        
        # Create data frame containing X 
        dfX = DataFrame(Intercept = X[:,1], Group = X[:,2], FishOil = X[:,3], InteractionGroupFishOil = X[:,4])
                
        # SAVE X MATRIX
        # dfX |> CSV.write("../data/dataprocessed/Xmat_CN_CS.csv")
    else
        # Fish Oil
        X[:,2] = (dfInput.FishOil .== "yes").*2 .- 1
        
        # Create data frame containing X
        dfX = DataFrame(Intercept = X[:,1], InteractionGroupFishOil = X[:,2])
        
        # SAVE X MATRIX
        # dfX |> CSV.write("../../data/dataprocessed/Xmat_CSbaseline_CS.csv")
    end
    
    #####
    # Y #
    ##### 
    # Generate Y matrix
    
    if isunpaired
        dfY = dfInput[:, 6:end];
        # SAVE Y MATRIX
#         df[:, 6:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv");
    else
        dfY = dfInput[:, 3:end];
        # SAVE Y MATRIX
#         df[:, 3:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv")
    end;   
    
    return dfX, dfY
end


"""
**getXYnoFishOil** -*Function*.

    getXYnoFishOil(dfInput::DataFrame; 
                        responseSelection = [""],
                        isunpaired::Bool = true) => DataFrame, DataFrame

Returns X covariates and Y response matrix.

"""
function getXYnoFishOil(dfInput::DataFrame; 
                        responseSelection = [""],
                        isunpaired::Bool = true)
    
    responseSelection = intersect(names(dfInput), responseSelection) 
    # Select responses
    if isunpaired
        colSelection = vcat(["Sample", "Batch", "Group", "Statin", "FishOil"], responseSelection)
    else
        colSelection = vcat(["Sample", "FishOil"], responseSelection)        
    end;
    dfInput= dfInput[:, colSelection]
    
    # get dimensions of sample and covariates
    n = size(dfInput)[1]
    p = ifelse(isunpaired, 2, 1); # Intercept | Group  
    
    #####
    # X #
    ##### 
    # Generate X matrix
    X = zeros(n,p)
    # Intercept
    X[:,1] = ones(n,1)
    if isunpaired
        # Group
#         X[:,2] = (dfInput.Group .== "CS").*2 .- 1 
        X[:,2] = (dfInput.Group .== "CS").*1 
        # SAVE X MATRIX
        dfX = DataFrame(Intercept = X[:,1], Group = X[:,2])
        # dfX |> CSV.write("../data/dataprocessed/Xmat_CN_CS.csv")
    else
        # SAVE X MATRIX
        # dfX = DataFrame(Intercept = X[:,1], InteractionGroupFishOil = X[:,2])
        # dfX |> CSV.write("../../data/dataprocessed/Xmat_CSbaseline_CS.csv")
    end
    
    #####
    # Y #
    ##### 
    # Generate Y matrix
    
    if isunpaired
        dfY = dfInput[:, 6:end];
        # SAVE Y MATRIX
        # df[:, 6:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv");
    else
        dfY = dfInput[:, 3:end];
        # SAVE Y MATRIX
        # df[:, 3:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv")
    end;   
    
    return dfX, dfY
end


"""
**getXY4ways** -*Function*.

    getXY4ways(dfInput::DataFrame; 
                        responseSelection = [""],
                        isunpaired::Bool = true) => DataFrame, DataFrame

Returns X covariates and Y response matrix.

"""
function getXY4ways(dfInput::DataFrame; 
                        responseSelection = [""],
                        isunpaired::Bool = true)
    
    responseSelection = intersect(names(dfInput), responseSelection) 
    # Select responses
    if isunpaired
        colSelection = vcat(["Sample", "Batch", "Group", "Statin", "FishOil"], responseSelection)
    else
        colSelection = vcat(["Sample", "FishOil"], responseSelection)        
    end;
    dfInput= dfInput[:, colSelection]
    
    # get dimensions of sample and covariates
    n = size(dfInput)[1]
    p = ifelse(isunpaired, 2, 1); # Intercept | Group  
    
    #####
    # X #
    ##### 
    # Generate X matrix
    # Formula for model matrix
    frmlStr = string("@formula(", names(dfInput)[6], " ~ 1 + Group + FishOil)");
    frml = eval(Meta.parse(frmlStr));
    # println(frml)
    # Generate model matrix
    mf = ModelFrame(frml, @view(dfInput[:,1:6]))
    mX = modelmatrix(mf)
       
    
    
    if isunpaired
        # SAVE X MATRIX
        dfX = DataFrame(CN_NoFishOil = ones(n),
                        CN_FishOil = iszero.(isone.(mX[:,2]) .| .~isone.(mX[:,3]))*1.0,
                        CS_NoFishOil = iszero.(.~isone.(mX[:,2]) .| isone.(mX[:,3]))*1.0,
                        CS_FishOil = (isone.(mX[:,2]) .& isone.(mX[:,3]))*1.0)

        

#         dfX = DataFrame(CN_NoFishOil = iszero.(isone.(mX[:,2]) .| isone.(mX[:,3])).*1,
#                         CN_FishOil = iszero.(isone.(mX[:,2]) .| .~isone.(mX[:,3])).*1,
#                         CS_NoFishOil = iszero.(.~isone.(mX[:,2]) .| isone.(mX[:,3])).*1,
#                         CS_FishOil = (isone.(mX[:,2]) .& isone.(mX[:,3])).*1)

        # dfX |> CSV.write("../data/dataprocessed/Xmat_CN_CS.csv")
    else
        # SAVE X MATRIX
        # dfX = DataFrame(Intercept = X[:,1], InteractionGroupFishOil = X[:,2])
        # dfX |> CSV.write("../../data/dataprocessed/Xmat_CSbaseline_CS.csv")
    end
    
    #####
    # Y #
    ##### 
    # Generate Y matrix
    
    if isunpaired
        dfY = dfInput[:, 6:end];
        # SAVE Y MATRIX
        # df[:, 6:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv");
    else
        dfY = dfInput[:, 3:end];
        # SAVE Y MATRIX
        # df[:, 3:end] |> CSV.write("../data/dataprocessed/Ymat_TG.csv")
    end;   
    
    return dfX, dfY
end



"""
**getCoefs** -*Function*.

    getCoefs(dfInput::DataFrame, Zin::Array{Float64,2}; 
                    responseSelection = [""],
                    isunpaired::Bool = true,
                    hasFishOil::Bool = true) => Array{Float,2}, Array{Float,2}, Array{Float,2}

Returns coefficients, confidence of interval, T-stats.
"""
function getCoefs(dfInput::DataFrame, Zin::Array{Float64,2}; 
                    responseSelection = [""],
                    isunpaired::Bool = true, 
                    hasFishOil::Bool = true, is4ways = false)
       
    ###########
    # X and Y #
    ###########
    # responseSelection = intersect(names(dfInput), responseSelection) 
    if hasFishOil
        if is4ways
            dfX, dfY =  getXY4ways(dfInput; responseSelection = responseSelection, isunpaired = isunpaired)
        else
            dfX, dfY =  getXY(dfInput, responseSelection = responseSelection, isunpaired = isunpaired)
        end
    else
        dfX, dfY =  getXYnoFishOil(dfInput, responseSelection = responseSelection, isunpaired = isunpaired)
    end
    X = Matrix(dfX);
    Y = Matrix(dfY);
    
    #########
    # MODEL #
    #########
    
    # Construct a RawData object
    dat = RawData(Response(Y), Predictors(X,Zin));

    # Estimate coefficients
    est = mlm(dat, addXIntercept=false, addZIntercept=false);
   
    # get estimate
    estCoefOut = MatrixLM.coef(est);

    # get variance estimate
    σ² = est.varB
    # get SE
    SE = sqrt.(σ²)
    CIOut = SE.*1.96
    
    # get T-statistics
    tStatsOut = t_stat(est)
    
    return estCoefOut, CIOut, tStatsOut, σ² 
end

"""
**getPermPvals** -*Function*.

    getPermPvals(dfInput::DataFrame, Zin::Array{Float64,2}; 
                    responseSelection = [""],
                    isunpaired::Bool = true,
                    hasFishOil::Bool = true) => Array{Float,2}, Array{Float,2}

Returns T-statistics and corresponding P-value based on permutations tests.
"""
function getPermPvals(dfInput::DataFrame, Zin::Array{Float64,2};
                    nPerms = 100,
                    responseSelection = [""],
                    isunpaired::Bool = true, 
                    hasFishOil::Bool = true, is4ways = false)
       
    ###########
    # X and Y #
    ########### 
    if hasFishOil
        if is4ways
            dfX, dfY =  getXY4ways(dfInput; responseSelection = responseSelection, isunpaired = isunpaired)
        else
            dfX, dfY =  getXY(dfInput, responseSelection = responseSelection, isunpaired = isunpaired)
        end
    else
        dfX, dfY =  getXYnoFishOil(dfInput, responseSelection = responseSelection, isunpaired = isunpaired)
    end
    X = Matrix(dfX);
    Y = Matrix(dfY);
    
    #########
    # MODEL #
    #########
    
    # Construct a RawData object
    dat = RawData(Response(Y), Predictors(X,Zin));

    # Estimate t-stats and pVals
    tStatsOut, pValsOut = mlm_perms(dat, nPerms, addXIntercept=false, addZIntercept=false);
    
    return tStatsOut, pValsOut 
end




