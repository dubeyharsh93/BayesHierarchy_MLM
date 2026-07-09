#=
Synopsosis: `utils.jl`.

Functions list:


- tstat2pval:
    Returns p-values from the t-statistics and degree of freedom.

- pval2qval:
    Returns q-values from the p-values.
    
=#

"""
**tstat2pval** -*Function*

    tstat2pval(mTstats::AbstractMatrix, df::Int64; istwotailed = true) => Matrix

Returns p-values from the t-statistics and degree of freedom.

***Arguments***

- `mTstats` vector.
- `df` degree of freedom.
- `istwotailed` indicate if two-tailed test or one tail test.
"""
function tstat2pval(mTstats::AbstractMatrix, df::Int64; istwotailed = true)
        
    mTstats = permutedims(mTstats)
    
    # initialize pValue matrix
    q, p = size(mTstats) # size of B' 
    pVals = zeros(q, p)
    
    if istwotailed
       a = 2
    else
        a = 1
    end

    for i in 1:p
       pVals[:,i] = ccdf(TDist(df), abs.(mTstats[:,i])).*a
    end
    
    return permutedims(pVals)
end


"""
**pval2qval** -*Function*

    pval2qval(mTstats::AbstractMatrix, df::Int64; istwotailed = true) => Matrix

Returns q-values from the p-values.

"""
function pval2qval(vPvals::AbstractVector)
        
    # Step 1: Estimate proportion of truly null hypothesis
    n = length(vPvals)
    Π₀ = estimate(vPvals, StoreyBootstrap())
    
    # Step 2: sort original p-values
    sort!(vPvals)

    # Step 3: calculate the q-value for the largest p-value
    vQvals = copy(vPvals)
    vQvals[end] = vQvals[end]*Π₀

    # Step 4: calculate the rest of the q-values such as
    for i in n-1:-1:1
       vQvals[i] = minimum([(Π₀*n*vPvals[i])/i, vQvals[i+1]]) 
    end
    
    return vQvals
end


function pval2qval2(vPvals::AbstractVector)
    
    # Step 1: Estimate proportion of truly null hypothesis
    n = length(vPvals)
    Π₀ = estimate(vPvals, StoreyBootstrap())
    
    vQvals = adjust(vPvals, BenjaminiHochbergAdaptive(Π₀))
        
    return vQvals
end




function pval2qval(mPvals::AbstractMatrix)
    
    mPvals = permutedims(mPvals)
    
    # initialize mQvals matrix
    q, p = size(mPvals) # size of B' 
    mQvals = zeros(q, p)
    
    
    
    for i in 1:p
       mQvals[:,i] = pval2qval2(mPvals[:, i])
    end
    
    return permutedims(mQvals)
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