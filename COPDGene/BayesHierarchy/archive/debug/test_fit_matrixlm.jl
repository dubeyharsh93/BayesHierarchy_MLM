using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "construct_matrices.jl"))
include(joinpath(@__DIR__, "fit_matrixlm.jl"))

obj = construct_copdgene_matrices()

mlm_res = fit_copdgene_matrixlm(obj.X, obj.Y, obj.Z)

println("MatrixLM coefficient size: ", size(mlm_res.coef))
println("MatrixLM t-stat size: ", size(mlm_res.tstat))
println("MatrixLM SE size: ", size(mlm_res.se))

println("First few coefficient names:")
println(obj.coef_names)

println("First 5 estimates for Age:")
k_age = findfirst(==("Age"), obj.coef_names)
println(mlm_res.coef[k_age, 1:5])

println("First 5 SEs for Age:")
println(mlm_res.se[k_age, 1:5])