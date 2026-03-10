include("src/CPSS2026B1.jl")

using .CPSS2026B1

const PARAMS_FILE = "params.dat"

params = read_parameters(PARAMS_FILE)
beta = 1 / params.temperature

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

accepted = metropolis_sweep!(lattice; β = beta, J = params.coupling)
println("accepted flips = ", accepted)
println("updated magnetization = ", magnetization(lattice))
println("updated energy = ", energy(lattice; J = params.coupling))
