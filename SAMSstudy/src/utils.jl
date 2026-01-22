# using Catlab.WiringDiagrams, Catlab.Graphics
# using Catlab.Theories


"""
**crossRefLipids** -*Function*.

    crossRefLipids(df::DataFrame, baseName::String = "lip", startCol::Int64 = 5) => Dataframe, DataFrame

Create a cross reference for the lipids name; returns a data frame with new names
and a look up table.
"""
function crossRefLipids(df::DataFrame, 
                        baseName::String = "lip"; 
                        startCol::Int64 = 5)
    
    numLipids = size(df)[2] - startCol + 1 # first columns are not lipids
    newNames = names(df)
    newNames[startCol:end] = repeat([baseName], numLipids).*string.(collect(1:numLipids))
    
    dfCrossRef = DataFrame(lipID = newNames[startCol:end], 
                            OriginalNames = names(df)[startCol:end] )
    
    DataFrames.rename!(df, Symbol.(newNames))
    
    return df, dfCrossRef    
end 




"""
**addGroupCat** -*Function*.

    addGroupCat(df::DataFrame; categorical::Bool = true) => DataFrame

Create a Group variable: CN, CS, CSbaseline, based on the lipidomic data.
"""
function addGroupCat(df::DataFrame, categorical::Bool = true)
    # initialize a group vector
    vGroup = Array{String}(undef, size(df)[1])
    
    # create CN category
    idxCN = findall(x -> occursin("CN", x), df.Sample[:])
    vGroup[idxCN] .= "CN"
    
    # create CS category
    idxCS = findall(occursin.("CS", df.Sample[:]) .& (df.Statin[:].=="Statin"))
    vGroup[idxCS] .= "CS"

    # create CSbaseline category
    idxCSbaseline = findall(occursin.("CS", df.Sample[:]) .& (df.Statin[:].!="Statin"))
    vGroup[idxCSbaseline] .= "CSbaseline"

    # convert string to categorical array
    if categorical
        vGroup = CategoricalArray(vGroup)
    end


    if any(occursin.("Group", names(df)))
        df[!, :Group] = vGroup
        else 
#         insert!(df, 3, vGroup, :Group);
        insertcols!(df, 3, :Group => vGroup)
    end
    
    return df
end

"""
**addGroupCatMeta** -*Function*.

    addGroupCatMeta(df::DataFrame; categorical::Bool = true) => DataFrame

Create a Group variable: CN, CS, CSbaseline, based on the metabolomic data.
"""
function addGroupCatMeta(df::DataFrame, categorical::Bool = true)
    # initialize a group vector
    vGroup = Array{String}(undef, size(df)[1])
    
    # create CN category
    idxCN = findall(x -> occursin("CN", x), df.Sample[:])
    vGroup[idxCN] .= "CN"
    
    # create CS category
    idxCS = findall(occursin.("CS", df.Sample[:]) .& (df.Group[:].!="Baseline"))
    vGroup[idxCS] .= "CS"

    # create CSbaseline category
    idxCSbaseline = findall(occursin.("CS", df.Sample[:]) .& (df.Group[:].=="Baseline"))
    vGroup[idxCSbaseline] .= "CSbaseline"

    # convert string to categorical array
    if categorical
        vGroup = CategoricalArray(vGroup)
    end


    if any(occursin.("Group", names(df)))
        df[!, :Group] = vGroup
        else 
#         insert!(df, 3, vGroup, :Group);
        insertcols!(df, 3, :Group => vGroup)
    end
    
    return df
end





"""
**catBGSF** -*Function*.

    catBGSF(df::DataFrame) => DataFrame

Converts Batch, Statin and Oil variable to categorical array type variable.
"""
function catBGSF(df::DataFrame, con2cat::Bool = true)
    # Batch
    if con2cat
        df.Batch = CategoricalArray(df.Batch)
    end
    
    # Group
    if con2cat
        df.Group = CategoricalArray(df.Group)
    end
    
    # Statin
    if any(occursin.("Statin", names(df)))    
        idxYes = findall(x -> occursin("Statin", x), df.Statin[:])
        idxNo = findall(x -> occursin("No Statin", x), df.Statin[:])
        df.Statin[idxYes] .= "yes"
        df.Statin[idxNo] .= "no"
        
        if con2cat
            df.Statin = CategoricalArray(df.Statin)
        end
    else
        # initialize a statin vector
        vStatin = Array{String}(undef, size(df)[1])

        # create Statin: yes category
        idxStatin = findall((df.Group[:].!="CS") .| (df.Group[:].!="CN"))
        vStatin[idxStatin] .= "no"

        # create Statin: no category
        idxNoStatin = findall(df.Group[:].!="CSbaseline")
        vStatin[idxNoStatin] .= "yes"
        
        # convert to categorical
        if con2cat
            vStatin = CategoricalArray(vStatin)
        end 
        
#         insert!(df, 4, vStatin, :Statin);
        insertcols!(df, 4, :Statin => vStatin)
    end
    
    # Fish Oil
    vName = names(df)[1:6];
    idxFishOil = findall(x-> (x=="Fish Oil") | (x =="Oil") | (x == "Fish oil"), vName)
    
    if !isempty(idxFishOil)
        rename!(df, idxFishOil[1]=> :FishOil)
    end
    idxYes = findall(x -> occursin("Fish Oil", x), df.FishOil[:])
    idxNo = findall(x -> occursin("No Fish Oil", x), df.FishOil[:])
    df.FishOil[idxYes] .= "yes"
    df.FishOil[idxNo] .= "no"
    
    if con2cat
        df.FishOil = CategoricalArray(df.FishOil)
    end
    
    return df
end


"""
**findDuplicates** -*Function*.

    findDuplicates(mat::Array{Float64,2}) => Array{Float64,2}

Find duplicates in array of string.
"""
function findDuplicates(vS::Array{String,1})

    buffS = deepcopy(vS)
    dupIndices = []
    dupNames = []

    # Find potential duplicates

    while length(buffS) > 1
        idxDup = findall(occursin.(buffS[1], vS))
        if length(idxDup) < 2
            # No duplicate
            deleteat!(buffS, [1]) 
        else
            # Check distance for duplicates, it should be 0
            lngthDup = length.(vS[idxDup])
            dist = abs.(lngthDup .- minimum(lngthDup))
            idxDupFilt  = idxDup[findall(x -> x==0 || x == 2, dist)]
            if length(idxDupFilt) > 1
                push!(dupIndices, Tuple(idxDupFilt))
                push!(dupNames, vS[idxDupFilt[1]])
                deleteat!(buffS, findall(occursin.(buffS[1], buffS)))
            else
                # No duplicate
                deleteat!(buffS, [1]) 
            end
        end
        #println("There exist at least ", length(dupIndices), " duplicates.")
    end
    
    if isempty(dupIndices)
        println("No duplicates to be found!", length(dupIndices))
    else    
        println("There exist at least ", length(dupIndices), " duplicates.")
    end
    
    return dupIndices, dupNames

end


"""
**diGraph** -*Function*.

    diGraph(steps::Array{String,1}) => GraphViz diagram

Create a directional diagram.
"""
function diGraph(steps::Array{String,1})
    steps = Symbol.(steps)

    n = length(steps)

    ptrtmnt = Array{Any,1}(undef, n)
    A, B = Ob(FreeSymmetricMonoidalCategory, :A, :B)

    for i in 1:n    
        ptrtmnt[i] = Hom(steps[i],A,A)
    end  
 
    if n>= 2 
        composite = compose(ptrtmnt[1], ptrtmnt[2])
    else
        composite = ptrtmnt[1]
    end

    if n > 2
        for i in 3:n
            composite = compose(composite, ptrtmnt[i])
        end
    end
    
    to_graphviz(composite, orientation=LeftToRight, outer_ports=false,
        node_attrs = Dict(
            :fontname => "Arial",
            :center => "true",
        )
    
    )
    
end

"""
**findFullMissingColumn** -*Function*.

    findFullMissingColumn(df::DataFrame; startCol::Int64 = 1) => Array{Float64,2}

Find column where all values are missing.
"""
function findFullMissingColumn(df::DataFrame; startCol::Int64 = 1) 
    
        
    m2test = Matrix(df[:,startCol:end])
    m2test = Array{Union{Missing, Float64},2}(m2test);
    
    nRows = size(m2test)[1]
    
    # find cols with missing 
    idxColMissing = findall(x -> x!=0, (sum(ismissing.(m2test), dims = 1)[:]))
    
    # Initialize
    vFullMissing = zeros(length(idxColMissing)) 
    
    #
    if length(idxColMissing) != 0
        for i in 1:length(idxColMissing)
            vFullMissing[i] = length(findall(m2test[:, idxColMissing[i]] .=== missing)) == nRows
        end
    end
    
    rslt = idxColMissing[vFullMissing .== 1]
    
    return rslt

end

"""
**getLipName** -*Function*.

    getLipName(x::String) => String

Returns the first name of lipid ID in order to be identified by a lipid shorthand nomenclature grammar.  
"""
function getLipName(x::String)
    frstLip = strip(split(x, "|")[1])
    frstLip = replace(frstLip, "_" => "/")
    
    return strip(split(frstLip, "+")[1])[3:end]   
end

"""
    deviation_contrast_matrix(k::Integer, ref::Integer)

Constructs a deviation contrast matrix for a given number of levels `k` 
and a specified reference level `ref`.

# Arguments
- `k::Integer`: The total number of levels (columns in the resulting matrix).
- `ref::Integer`: The reference level (row in the matrix where `-1` is inserted).

# Returns
- A matrix of size `k×(k-1)` where:
  - All columns form an identity matrix (of size `(k-1)x(k-1)`),
  - A row of `-1` is added at the reference level.

# Example
```julia
julia> deviation_contrast_matrix(4, 2)
4×3 Matrix{Float64}:
  1.0   0.0   0.0
 -1.0  -1.0  -1.0
  0.0   1.0   0.0
  0.0   0.0   1.0
"""
function deviation_contrast_matrix(k::Integer, ref::Integer)
    mat = 1*Matrix(I, k-1, k-1)    
    if ref == 1
        return vcat(-1*ones(Int, 1, k-1), mat)
    elseif ref == k
        return vcat(mat , -1*ones(Int, 1, k-1))
    else
        return vcat(
            mat[1:(ref-1), :],
            -1*ones(Int, 1, k-1),
            mat[(ref):end, :]
        )     
    end 
end