using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__,"gibbs_src", "construct_sams_matrices.jl"))

obj = construct_sams_matrices()

println("X size: ", size(obj.X))
println("Y size: ", size(obj.Y))
println("Z size: ", size(obj.Z))
println("Coefficient names: ", obj.coef_names)
println("H: ", obj.H)
println("First 20 subclass labels: ", obj.subclass_of_met[1:20])
println("Group counts: ", [sum(obj.subclass_of_met .== h) for h in 1:obj.H])