using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "gibbs_src", "construct_sim_matrices.jl"))

obj = simulate_bilinear_identity_1level(
    n = 60,
    m = 300,
    p = 10,
    H = 4,
    prop_sub = [0.55, 0.25, 0.15, 0.05],
    theta0_scale = 0.2f0,
    tau_v = 0.12f0,
    tau_w = 0.08f0,
    sigma_y = 1.0f0,
    T = Float32,
    seed = 1
)

println("X size: ", size(obj.X))
println("Y size: ", size(obj.Y))
println("Z size: ", size(obj.Z))
println("B_true size: ", size(obj.B_true))
println("H: ", obj.H)
println("Subclass counts: ", [sum(obj.subclass_of_met .== h) for h in 1:obj.H])
println("Coefficient names: ", obj.coef_names)