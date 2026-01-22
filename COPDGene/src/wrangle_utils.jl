#=
Synopsosis: `wrangle_utils.jl` contains function to wrangle data set.

Functions list:

- keepComplete: 
    Returns a metabolites level dataframe with only individual's complete cases (i.e. no 
    missing in independent variables).

- readCOPDdata:
    Returns a metabolites COPD dataframe.

- moveIndices:
    Returns a vector where orginal vector elements were moved to a new index inside the orginal vector.
    
=#



"""
**keepComplete** -*Function*

    keepComplete(dfMet::DataFrame, dfInd::DataFrame, dfRef::DataFrame) => DataFrame

Returns a metabolites level dataframe with only individual's complete cases (i.e. no 
missing in independent variables).

***Arguments***


- `dfMet` dataframe contains original metabolites level, where columns are samples and rows
    are metabolites.
- `dfInd` dataframe contains original individual covariates, where columns are covraiates
    and rows are samples.
- `dfRef` dataframe contains metabolites references.
"""
function keepComplete(dfMet::DataFrame, dfInd::DataFrame, dfRef::DataFrame; 
                      sampleCol::Symbol = :SampleName, metaCol::Symbol = :MetaID)
    # select complete cases
    idxComplete = findall(completecases(dfInd))
    vComplete = dfInd[idxComplete, sampleCol];
    dfMet = select(dfMet, Symbol.(vcat([string(metaCol)], vComplete)));

    # dfNegMetabo CompID variable
    dfMet = innerjoin(dfRef[:, [metaCol, :CompID]], dfMet, on = metaCol) 
    dfMet = select(dfMet,Not([metaCol]));
    
    return dfMet
end


"""
**readCOPDdata** -*Function*

    readCOPDdata(fileName::String) => DataFrame

Returns a metabolites COPD dataframe.

***Arguments***

- `fileName` path of the file.
"""
function readCOPDdata(fileName::String)
    
    dfMetabo = CSV.read(fileName, DataFrame; missingstring = ["NA", ""]);
    dfMetabo = permutedims(dfMetabo, 1, :SampleID);
    
    return dfMetabo
end


"""
**moveIndices** -*Function*

    moveIndices(v::AbstractVector, rng::AbstractVector, loc::Int64 ) => Vector

Returns a vector where orginal vector elements were moved to a new index inside the orginal vector.

***Arguments***

- `v` vector.
- `rng` contains the indices of the element to be moved.
- `loc` new index.
"""
function moveIndices(v::AbstractVector, rng::AbstractVector, loc::Int64 )
    if rng[1] < loc
        vCirc = circshift(v, -loc -length(rng)+1)
        idxNew = findall(circshift(collect(1:length(v)),-loc -length(rng)+1).== rng[1])
        rngNew = collect((1:length(rng)).+(idxNew[1] -1))
        vNew =  circshift(vcat(vCirc[rngNew], vCirc[Not(rngNew)]), loc-1)
    else
        vCirc = circshift(v, -loc+1);
        idxNew = findall(circshift(collect(1:length(v)), -loc+1).== rng[1])
        rngNew = collect((1:length(rng)).+(idxNew[1] -1))
        rngNew = mod.(rngNew, length(v)+1)
        rngNew[rngNew.==0] .= 1
        vNew = circshift(vcat(vCirc[rngNew], vCirc[Not(rngNew)]), loc-1)
    end
    
    return vNew
end


"""
**prepend** -*Function*

    prepend(fileName::String, mystring::String)

Prepend a string into a text file; it keeps the original and generate a new file.

***Arguments***

- `fileName` contains the path name of the text file.
- `mystring` contains the string to prepend.
"""
function prepend(fileName::String, mystring::String)
    
    # Keep original file and make a new file
    splitFileName = splitdir(realpath(fileName));
    tempFilePath = joinpath(splitFileName[1], "new_"*splitFileName[2])
    
    # Prepend
    open(tempFilePath, "w") do file1
        open(fileName) do file2
            write(file1, mystring, read(file2, String))
        end
    end;
    
    
end


"""
**prepend!** -*Function*

    prepend!(fileName::String, mystring::String)

Prepend a string into a text file.

***Arguments***

- `fileName` contains the path name of the text file.
- `mystring` contains the string to prepend.
"""
function prepend!(fileName::String, mystring::String)
    
    # Keep original file and make a new file
    splitFileName = splitdir(realpath(fileName));
    tempFilePath = joinpath(splitFileName[1], "new_"*splitFileName[2])
    
    # Prepend
    open(tempFilePath, "w") do file1
        open(fileName) do file2
            write(file1, mystring, read(file2, String))
        end
    end;
    
    mv(tempFilePath, fileName, force = true)
end