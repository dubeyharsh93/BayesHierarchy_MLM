#=
Author: Gregory Farage
Date: 2022-04-28
Synopsosis: `pretreatment.jl` contains function to preprocess metabolomics data.

Functions list:

- imputeKNN: 
    Replaces missing elements based on k-nearest neighbors (KNN) imputation.

- imputeCOPD:
    Replaces missing elements based on k-nearest neighbors (KNN) imputation. 
    But the cotinine missing values are imputed with 0s. Variable with a
    missingness more the threshold value will be dropped.

- imputeHM:
    Replaces missing elements with half of the minimum of non-missing elements in the corresponding variable.

- getTransposedMetaboData:
    Returns the transposed dataframe containing metabolites values for each individual
    sample, where the columns are metabolites and the rows are individual samples.

=#


# using GLM, StatsBase, RCall


"""
**imputeKNN** -*Function*.

    imputeKNN(df::DataFrame; threshold = 20) => DataFrame

Replaces missing elements based on k-nearest neighbors (KNN) imputation. 
The missing values are imputed based on similarity among non-missing rows.

"""
function imputeKNN(df::DataFrame; threshold = 20)
    
    # Select variables with missingness of less or equal to threshold percentage value
    vSubjectID = df.SampleID;
    
    # Calculate missingness per column
    n = size(df, 1)
    mat = Matrix(df[:, 2:end])
    vMissing = map(eachcol(mat)) do col
                   sum(ismissing.(col))*100/n
               end
    
    # Get indices of the metabolites with a missingness less or equal to threshold %
    idxThresh = findall(vMissing .<= threshold);
    
    # Remove metabolites with more than threshold % of missing data
    # select!(dfNegMetabo, idxTresh20)
    mat = mat[:, idxThresh]
    
    # Impute using kNN algorithm
    data = permutedims(mat);
    @rput data;
    
    R"""
    suppressMessages(library(impute))

    # scale before imputation
    matScaledT <- scale(t(data))
    matScaled <- t(matScaledT)

    # imputation
    matImputed <- impute.knn(matScaled, k = 10)$data
    imputedData <- (matImputed * attr(matScaledT, 'scaled:scale') + attr(matScaledT, 'scaled:center'));

    """
    @rget imputedData;
    imputedData = permutedims(imputedData);
    
    # Remove metabolites with more than threshold % of missing data
    df = select(df, vcat([1],idxThresh .+ 1));
    df[:,2:end] = imputedData;
    
    return df

end   


"""
**imputeCOPD** -*Function*.

    imputeCOPD(df::DataFrame; threshold = 20) => DataFrame

Replaces missing elements based on k-nearest neighbors (KNN) imputation. 
But the cotinine missing values are imputed with 0s. Variable with a 
missingness more the threshold value will be dropped.

"""
function imputeCOPD(df::DataFrame; threshold = 20)
    
    # Preserve Cotinine variables
    # Get column indices for cotinine variables
    dfLookUp = DataFrame(CompID = names(df)[2:end])
    dfLookUp = leftjoin(dfLookUp, select(dfRef, [:CompID, :Biochemical]), on = :CompID);
    idxCotinine = findall(dfLookUp.Biochemical .== "cotinine")

    if !isempty(idxCotinine)
        # Create a dataframe that will be joined to the kNN imputed data
        idxColCotinine = findall(names(df) .== dfLookUp.CompID[idxCotinine])
        dfCotinine = df[:, vcat([1], idxColCotinine)];
        for i in 1:length(idxCotinine)
            metaName = dfLookUp.Biochemical[idxCotinine[i]]
            percMissing = round(100*sum(ismissing.(dfCotinine[:,i+1]))/size(dfCotinine, 1), digits = 2)
            println("The metabolite $(metaName) contains $(percMissing)% missing samples.")
            dfCotinine[:, i+1] = coalesce.(dfCotinine[:, i+1], 0);
        end
    end
    
    numMeta = size(df, 2)-1
    
    if !isempty(idxCotinine)
        df = select(df, Not(Symbol.(dfLookUp.CompID[idxCotinine])))
        df = imputeKNN(df; threshold = threshold)
        df = leftjoin(dfCotinine, df, on = :SampleID);
    else
        df = imputeKNN(df; threshold = threshold)
    end
    
    numMetaImputed = size(df, 2)-1
    
        println("We dropped $(numMeta - numMetaImputed) metabolites",
                " due to a missingness greater than $(threshold)%.")
        println("We preserved $(numMetaImputed) metabolites.")
    
    
    return df

end


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

    dfLip = convert(DataFrame, mLip);
    rename!(dfLip, Symbol.(names(df)[startCol:end]))
    
    if startCol > 1
        startCol = startCol-1 
        dfLip = hcat(df[:,1:startCol], dfLip)
    end
    
    return dfLip
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
# function getVarExpl(mMtblmcs::Array{Float64, 2}, Xbatch, nameCol) 
#     n = size(mMtblmcs)[1]
#     Xbatch = CategoricalArray(Xbatch)
    
#     # Initialize adjusted coefficient of determination for each lipids 
#     adjRsquared = zeros(n)

#     # Get adjusted R^2 for every lipids
#     for i in 1:n
#         data = DataFrame(X = Xbatch, Y = mMtblmcs[i,:])
#         out = lm(@formula(Y ~ X), data);
#         adjRsquared[i] = adjr2(out)
#     end

#     # If there are negative values force to 0
#     adjRsquared[findall(x -> x < 0, adjRsquared)] .= 0
    
#     # Create a tibble with varaition explained results 
#     dfVarExpl = DataFrame(Lipids = nameCol, VarExpl = adjRsquared)

#     # Order results
#     sort!(dfVarExpl, :VarExpl, rev=true)
    
#     return dfVarExpl
# end


"""
**`getVarExplPerMetaPerBatch`** -*Function*.

    getVarExplPerMetaPerBatch(mMtblmcs::Array{Float64, 2}, Xbatch, nameCol::Array{String,1}) => DataFrame  

Assess variance explained of every metabolites per batch.
"""
# function getVarExplPerMetaPerBatch(mMtblmcs::Array{Float64, 2}, Xbatch, nameMeta) 
    
#     # Get level batches
#     lvlBatches = levels(Xbatch)
#     n =  length(nameMeta)   

#     idxSorted = zeros(Int, n);
#     for i in 1:n
#        idxSorted[i] = parse(Int, match.(r"\d+", nameMeta[i]).match) 
#     end

#     batchLevel = fill(string("Not_B"), size(mMtblmcs)[2], length(lvlBatches))


#     for i in 1:length(lvlBatches)    
#         batchLevel[Xbatch .== string(lvlBatches[i]), i] .= lvlBatches[i]
#     end

#     # Initialize
#     adjRsquaredPerLipidsPerBatch = zeros(n, 4)

#     for j in 1:n

#         # Initialize
#         adjRsquaredBatch = zeros(length(lvlBatches))

#         for i in 1:length(lvlBatches)
#             data = DataFrame(X = CategoricalArray(batchLevel[:, i]), Y = mMtblmcs[idxSorted[j],:])
#             out = lm(@formula(Y ~ X), data);
#             adjRsquaredBatch[i]= adjr2(out)
#         end

#         # If there are negative values force to 0
#         adjRsquaredBatch[findall(x -> x < 0, adjRsquaredBatch)] .= 0


#         # Insert results in list
#         adjRsquaredPerLipidsPerBatch[j, :] = adjRsquaredBatch
#     end 
    
#     return adjRsquaredPerLipidsPerBatch
    
    
# end