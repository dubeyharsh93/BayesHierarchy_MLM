using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "gibbs_src", "construct_pans_matrices.jl"))

obj = construct_pans_matrices()

println("X size: ", size(obj.X))
println("Y size: ", size(obj.Y))
println("Z size: ", size(obj.Z))
println("Coefficient names:")
for (i, nm) in enumerate(obj.coef_names)
    println(i, " : ", nm)
end

println("H: ", obj.H)
println("G: ", obj.G)
println("First 20 subclass_of_met: ", obj.subclass_of_met[1:20])
println("First 20 superclass_of_met: ", obj.superclass_of_met[1:20])
println("Sum subclass_of_met: ", sum(obj.subclass_of_met))
println("Min subclass_of_met: ", minimum(obj.subclass_of_met))
println("Max subclass_of_met: ", maximum(obj.subclass_of_met))