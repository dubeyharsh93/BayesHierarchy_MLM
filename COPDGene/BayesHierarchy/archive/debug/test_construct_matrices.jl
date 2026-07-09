using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "construct_matrices.jl"))

obj = construct_copdgene_matrices()
println("X size: ", size(obj.X))
println("Y size: ", size(obj.Y))
println("Z size: ", size(obj.Z))
println("Number of coefficient names: ", length(obj.coef_names))
println("Coefficient names:")
println(obj.coef_names)
println("Non-site covariate indices:")
println(obj.idx_covar)