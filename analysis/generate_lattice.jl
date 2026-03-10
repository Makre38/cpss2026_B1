include(joinpath(@__DIR__, "..", "src", "CPSS2026Ising.jl"))

using .CPSS2026Ising

const PARAMS_FILE = joinpath(@__DIR__, "params.dat")

params = read_parameters(PARAMS_FILE)
lattice = generate_lattice(params.L)

println("L = ", params.L)
println("T = ", params.temperature)
println("J = ", params.coupling)
println("lattice =")
for row in eachrow(lattice)
    println(join(row, " "))
end
println("magnetization = ", magnetization(lattice))
println("energy = ", energy(lattice; J = params.coupling))
