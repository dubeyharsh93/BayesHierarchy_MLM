using GLM

"""
**imputeHM** -*Function*.

    imputeHM(df::DataFrame; startCol::Int64 = 1) => DataFrame

Replaces missing elements with half of the minimum of non-missing elements in the corresponding variable.

**Example:**  
```
julia> df = DataFrame(A = [1, 2, 3], 
          B = [missing, missing, missing],
          C = [missing, 4, 5],
          D = [6, missing, 7],
          E = [missing, missing, 10])

3×5 DataFrame
│ Row │ A     │ B       │ C       │ D       │ E       │
│     │ Int64 │ Missing │ Int64?  │ Int64?  │ Int64?  │
├─────┼───────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ 1     │ missing │ missing │ 6       │ missing │
│ 2   │ 2     │ missing │ 4       │ missing │ missing │
│ 3   │ 3     │ missing │ 5       │ 7       │ 10      │

julia> imputeHM(tt)

3×5 DataFrame
│ Row │ A     │ B       │ C       │ D       │ E       │
│     │ Int64 │ Missing │ Int64?  │ Int64?  │ Int64?  │
├─────┼───────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ 1.0   │ 0.5     │ 0.5     │ 6.0     │ 0.5     │
│ 2   │ 2.0   │ 1.0     │ 4.0     │ 1.0     │ 1.0     │
│ 3   │ 3.0   │ 1.5     │ 5.0     │ 7.0     │ 1.5     │

```

"""
function imputeHM(df::DataFrame; startCol::Int64 = 1)
    
    mLip = Matrix(df[:,startCol:end])
    mLip = Array{Union{Missing, Float64},2}(mLip);

    # find rows with missing 
    idxRowMissing = findall(x -> x!=0, (sum(ismissing.(mLip), dims = 2)[:]))
    
#     # replace missing values
#     if length(idxRowMissing) != 0
#         for i in 1:length(idxRowMissing)
#             mLip[i,:] = collect(Missings.replace(mLip[i,:], minimum(skipmissing(mLip[i,:]))/2))
#         end
#     end
    
        # replace missing values
    if length(idxRowMissing) != 0
        for i in 1:length(idxRowMissing)
            idx = idxRowMissing[i]
            mLip[idx,:] = collect(Missings.replace(mLip[idx,:], minimum(skipmissing(mLip[idx,:]))/2))
        end
    end
    

    mLip = Array{Float64,2}(mLip)

    dfLip = DataFrame(mLip, :auto); # convert(DataFrame, mLip);
    rename!(dfLip, Symbol.(names(df)[startCol:end]))
    
    if startCol > 1
        startCol = startCol-1 
        dfLip = hcat(df[:,1:startCol], dfLip)
    end
    
    return dfLip
end

"""
**imputeQRILC** -*Function*.
    
    imputeQRILC(df::DataFrame; tune_sigma=1.0, startCol::Int64 = 1, logTX::Bool = false) => DataFrame

Returns imputated dataframe based on quantile regression. 
The function is based on the function `impute.QRILC` from the `imputeLCMD.R`
package.

`matdata`: matrix of omics value, e.g., metabolomics matrix, where the rows are the features 
           and the columns are the samples.
`tune_sigma`: coefficient that controls the sd of the MNAR distribution
                        - tune_sigma = 1 if the complete data distribution is supposed 
                          to be gaussian
                        - 0 < tune_sigma < 1 if the complete data distribution is  
                          supposed to be left-censored
`startCol`: indicate starting column in the dataframe
`logTX`: if `true`, apply a logarithm transformation before to impute then 
         apply exponential transformation, default is `false`.
`threshold`: maximum proportion of missing value of sample accpeted, if there exist more than 
             `threshold` of samples, then features may be removed from the data. The defaults is 
             `0.2` which is called "80% rule".
"""
function imputeQRILC(df::DataFrame; 
        tune_sigma=1.0, 
        startCol::Int64 = 1, 
        logTX::Bool = false, 
        threshold = 0.2)
 
    df[!,startCol:end] = convert.(Union{Missing, Float64}, df[!,startCol:end])
    
    mLip = Matrix(df[:,startCol:end])
    mLip = Array{Union{Missing, Float64},2}(mLip);

    # Select variables with missingness of less or equal to threshold percentage value
    # Calculate missingness per column
    n = size(mLip, 1)
    vMissing = map(eachcol(mLip)) do col
                   sum(ismissing.(col))/n
               end
    
    # Get indices of the metabolites with a missingness less or equal to threshold %
    idxThresh = findall(vMissing .<= threshold);
    
    # Remove metabolites with more than threshold % of missing data
    mLip = mLip[:, idxThresh]

    # log transform if true
    if logTX
        mLip = log.(mLip.+1);
    end
    
    # Impute using quantile regression algorithm    
    data = permutedims(mLip);   
    imputed_data = impute_QRILC_R(data; tune_sigma=tune_sigma) |> permutedims

    # exponential transform if true
    if logTX
        imputed_data = exp.(imputed_data).-1;
    end
    
    # Remove metabolites with more than threshold % of missing data
    Δx = startCol - 1
    if Δx > 0
        df = select(df, vcat(collect(1:Δx), idxThresh .+ Δx));
        df[:,startCol:end] = imputed_data;
    else
        df = select(df, idxThresh);
        df[:,:] = imputed_data;
    end

    df[!,startCol:end] = convert.(Float64, df[!,startCol:end])
    
    return df
end

"""
**impute_QRILC** -*Function*.
    
    impute_QRILC(matdata::Matrix{Union{Missing, Float64}}; tune_sigma=1.0) => Matrix

Returns imputated matrix based on quantile regression. 
The function is based on the function `impute.QRILC` from the `imputeLCMD.R`
package.

# Arguments

`matdata`: matrix of omics value, e.g., metabolomics matrix, where the rows are the features 
           and the columns are the samples.
`tune.sigma`: coefficient that controls the sd of the MNAR distribution
                        - tune.sigma = 1 if the complete data distribution is supposed 
                          to be gaussian
                        - 0 < tune.sigma < 1 if the complete data distribution is  
                          supposed to be left-censored
"""
function impute_QRILC(matdata::Matrix{Union{Missing, Float64}}; tune_sigma=1.0)
    # Get dimensions of the data
    nFeatures, nSamples = size(matdata)

    # Initialize the matrix of complete abundances
    data_imputed = copy(matdata)

    for i in 1:nSamples
        curr_sample = matdata[:, i]
    
        # Calculate the percentage of missing values
        pNAs = count(ismissing, curr_sample) / length(curr_sample)
        
        # Estimate the mean and standard deviation of the original distribution using quantile regression
        upper_q = 0.99
        
        q_normal = quantile(
            Normal(0, 1),
            LinRange(pNAs + 0.001, upper_q + 0.001, 100)
        )
        
        q_curr_sample = quantile(skipmissing(curr_sample), LinRange(0.001,upper_q + 0.001, 100))
        temp_QR = lm(hcat(ones(length(q_normal),1),reshape(q_normal,:, 1)), q_curr_sample)
    
        mean_CDD = coef(temp_QR)[1]
        sd_CDD = abs(coef(temp_QR)[2])
    
       # Generate data from a truncated normal distribution with the estimated parameters
        truncated_dist = truncated(
            Normal(mean_CDD, sd_CDD * tune_sigma); 
            upper = quantile(Normal(mean_CDD, sd_CDD), pNAs + 0.001)
        )
    
        curr_sample_imputed = copy(curr_sample)
        missing_idx = findall(ismissing, curr_sample)
        curr_sample_imputed[missing_idx] .= rand(truncated_dist, nFeatures)[missing_idx]
        data_imputed[:, i] = curr_sample_imputed
     end

    return data_imputed
end        

"""
**impute_QRILC_R** -*Function*.
    
    impute_QRILC_R(matdata::Matrix{Union{Missing, Float64}}; tune_sigma=1.0) => Matrix

Returns imputated matrix based on quantile regression. 
Imputation is performed in R using `RCall` and package `imputeLCMD.R`.

# Arguments

`matdata`: matrix of omics value, e.g., metabolomics matrix, where the rows are the features 
           and the columns are the samples.
`tune_sigma`: coefficient that controls the sd of the MNAR distribution
                        - tune_sigma = 1 if the complete data distribution is supposed 
                          to be gaussian
                        - 0 < tune_sigma < 1 if the complete data distribution is  
                          supposed to be left-censored
"""
function impute_QRILC_R(matdata; tune_sigma=1.0) 
    @rput matdata;
    @rput tune_sigma;
    
    R"""
    suppressMessages(library(imputeLCMD))
    
    # imputation
    imputed_data <- impute.QRILC(matdata, tune.sigma = tune_sigma)[[1]];
    
    """;
    @rget imputed_data;

    return imputed_data
end

"""
**log2tx** -*Function*.

    log2tx(df::DataFrame; startCol::Int64 = 1) => DataFrame

Computes logarithm base 2 on dataframe.
"""
function log2tx(df::DataFrame; startCol::Int64 = 1)   
    
    if any(any.(eachcol(df[:,startCol:end].==0)))
        df[:,startCol:end] = log2.(df[:,startCol:end].+1)
    else
        df[:,startCol:end] = log2.(df[:,startCol:end])
    end
        
    return df
end


"""
**center** -*Function*.

    center(df::DataFrame) => DataFrame

Mean centering for each row.

Example:
# Test
tt = [0.5 1 2 3 3.5;
      7 3 5 1.5 3.5;
      8 2 5 6 9]
ttN = center(tt)
display(ttN)
mean(ttN, dims = 1)

"""
function center(df::DataFrame; startCol = 1)
    mLip = Matrix(df[:,startCol:end])
    
    vMean = mean(mLip, dims = 1)
    
    df[:,startCol:end] = mLip - repeat(vMean, size(mLip)[1])
    
    return df
end




"""
**intnorm** -*Function*.

    intnorm(mat::Array{Float64,2}) => Array{Float64,2}

Total Area Normalization for each row.

Example:
# Test
tt = [0.5 1 2 3 3.5;
      7 3 5 1.5 3.5;
      8 2 5 6 9]
ttN = intNorm(tt)
display(ttN)
sum(ttN, dims = 2)

"""
function intnorm(mat::Array{Float64,2})
    # total area normalization
    cstIntegral = 0.01
    mat = mat./(cstIntegral.*sum(mat, dims = 2))    
    
    return mat
end



"""
**`pqnorm`** -*Function*.

    pqnorm(df::DataFrame; startCol = 1) => DataFrame

Performs a probabilistic quotient normalization (PQN) for sample intensities.

Example:

# Test
tt = [0.5 1 2 3 3.5;
      7 3 5 1.5 3.5;
      8 2 5 6 9]
ttN = pqnorm(tt)
display(ttN)
sum(ttN, dims = 2)


"""
function pqnorm(df::DataFrame; startCol = 1)
    mat = Matrix(df[:,startCol:end])
    
    # Integral normalization
    mat = intnorm(mat)
    
    # Calculate the reference spectrum (default: median) of all samples
    refSpec = median(mat, dims = 1)
    
    # Calculate the quotients of all variables of interest of the test spectrum 
    # with those of the reference spectrum.
    qtnts = mat ./ refSpec
    
    # Calculate the median of these quotients.
    medQ = median(qtnts, dims= 2)
    
    # Divide all variables of the test spectrum by this median
    matNorm = mat./medQ
    
    df[:,startCol:end] = matNorm
    
    return df
end



"""
**`getVarExpl`** -*Function*.

    getVarExpl(mMtblmcs::Array{Float64, 2}, Xbatch, nameCol::Array{String,1}) => DataFrame  

Assess variance explained of every metabolites according to the Xbatch variable.
"""
function getVarExpl(mMtblmcs::Array{Float64, 2}, Xbatch, nameCol) 
    n = size(mMtblmcs)[1]
    Xbatch = CategoricalArray(string.(Xbatch))
    
    # Initialize adjusted coefficient of determination for each lipids 
    adjRsquared = zeros(n)

    # Get adjusted R^2 for every lipids
    for i in 1:n
        data = DataFrame(X = Xbatch, Y = mMtblmcs[i,:])
        out = lm(@formula(Y ~ X), data);
        adjRsquared[i] = adjr2(out)
    end

    # If there are negative values force to 0
    adjRsquared[findall(x -> x < 0, adjRsquared)] .= 0
    
    # Create a tibble with variation explained results 
    dfVarExpl = DataFrame(Lipids = nameCol, VarExpl = adjRsquared)

    # Order results
    sort!(dfVarExpl, :VarExpl, rev=true)
    
    return dfVarExpl
end


"""
**`getVarExplPerMetaPerBatch`** -*Function*.

    getVarExplPerMetaPerBatch(mMtblmcs::Array{Float64, 2}, Xbatch, nameCol::Array{String,1}) => DataFrame  

Assess variance explained of every metabolites per batch.
"""
function getVarExplPerMetaPerBatch(mMtblmcs::Array{Float64, 2}, Xbatch, nameMeta) 
    
    # Get level batches
    lvlBatches = levels(Xbatch)
    n =  length(nameMeta)   

    idxSorted = zeros(Int, n);
    for i in 1:n
       idxSorted[i] = parse(Int, match.(r"\d+", nameMeta[i]).match) 
    end

    batchLevel = fill(string("Not_B"), size(mMtblmcs)[2], length(lvlBatches))


    for i in 1:length(lvlBatches)    
        batchLevel[Xbatch .== string(lvlBatches[i]), i] .= lvlBatches[i]
    end

    # Initialize
    adjRsquaredPerLipidsPerBatch = zeros(n, 4)

    for j in 1:n

        # Initialize
        adjRsquaredBatch = zeros(length(lvlBatches))

        for i in 1:length(lvlBatches)
            data = DataFrame(X = CategoricalArray(batchLevel[:, i]), Y = mMtblmcs[idxSorted[j],:])
            out = lm(@formula(Y ~ X), data);
            adjRsquaredBatch[i]= adjr2(out)
        end

        # If there are negative values force to 0
        adjRsquaredBatch[findall(x -> x < 0, adjRsquaredBatch)] .= 0


        # Insert results in list
        adjRsquaredPerLipidsPerBatch[j, :] = adjRsquaredBatch
    end 
    
    return adjRsquaredPerLipidsPerBatch
    
    
end