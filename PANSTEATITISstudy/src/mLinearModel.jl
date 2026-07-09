
#############
# LIBRARIES #
#############
# using MatrixLM
# using CSV, DataFrames, DataFramesMeta, Missings, CategoricalArrays
# using StatsBase, Statistics, MultivariateStats, MatrixLM, Random, Distributions, StatsModels
# using LinearAlgebra
# using FreqTables, Plots


# #################
# FUNCTIONS LIST #
# #################

#=
Functions list 

- getCoefs:
    Returns coefficients, confidence of interval, T-stats.

- getPermPvals:
    Returns T-statistics and corresponding P-value based on permutations tests.

- funStandardize: 
    Standardize data, with the option of centering or not.

=#


# ############
# FUNCTIONS #
# ############




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
    tStatsOut, pValsOut = mlm_perms(dat, nPerms, hasXIntercept=false, hasZIntercept=false);
    
    return tStatsOut, pValsOut 
end


"""
**funStandardize** -*Function*.

    funStandardize(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, with the option of centering or not. 
"""
function funStandardize(dfInput::DataFrame; tocenter::Bool = true)
    dfOutput = deepcopy(dfInput);
    if tocenter
        A = Matrix(dfOutput[:,2:end]);
        dt  = fit(ZScoreTransform, A, dims=1)
        Astndz = StatsBase.transform(dt, A)
        dfOutput[:,2:end] = Astndz
    else 
        # No centering for paired model
        A = Matrix(dfOutput[:,2:end]);
        dt  = fit(ZScoreTransform, A, dims=1, center = false)
        Astndz = StatsBase.transform(dt, A)
        dfOutput[:,2:end] = Astndz
    end
    
    return dfOutput
end

"""
**funStandardize!** -*Function*.

    funStandardize!(df::DataFrame, unpaired::Bool = false) => DataFrame

Standardize data, with the option of centering or not.
"""
function funStandardize!(dfInput::DataFrame; tocenter::Bool = true)
    if tocenter
        A = Matrix(dfInput[:,2:end]);
        dt  = fit(ZScoreTransform, A, dims=1)
        Astndz = StatsBase.transform(dt, A)
        dfInput[:,2:end] = Astndz
    else 
        # No centering for paired model
        A = Matrix(dfInput[:,2:end]);
        dt  = fit(ZScoreTransform, A, dims=1, center = false)
        Astndz = StatsBase.transform(dt, A)
        dfInput[:,2:end] = Astndz
    end
    
    return dfInput
end