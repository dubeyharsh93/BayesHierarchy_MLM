############################################################
# train_test_split.jl
#
# Purpose:
#   Split data into training and testing sets.
############################################################

using Random

function train_test_split_rows(Y, X; train_frac = 0.70, seed = 2026, shuffle = true)
    n = size(Y, 1)
    @assert size(X, 1) == n "X and Y must have same number of rows."

    rng = MersenneTwister(seed)
    idx = collect(1:n)

    if shuffle
        Random.shuffle!(rng, idx)
    end

    n_train = floor(Int, train_frac * n)

    idx_train = idx[1:n_train]
    idx_test  = idx[(n_train + 1):end]

    return (
        Y_train = Matrix(Y[idx_train, :]),
        X_train = Matrix(X[idx_train, :]),
        Y_test = Matrix(Y[idx_test, :]),
        X_test = Matrix(X[idx_test, :]),
        idx_train = idx_train,
        idx_test = idx_test
    )
end