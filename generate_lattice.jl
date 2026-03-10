include("src/CPSS2026B1.jl")

using .CPSS2026B1

const PARAMS_FILE = "params.dat"
const BETA = 0.4

lattice = generate_lattice_from_file(PARAMS_FILE)

println("L = ", size(lattice, 1))
println("lattice =")
for row in eachrow(lattice)
    println(join(row, " "))
end
println("magnetization = ", magnetization(lattice))
println("energy = ", energy(lattice))

accepted = metropolis_sweep!(lattice; β = BETA)
println("accepted flips = ", accepted)
println("updated magnetization = ", magnetization(lattice))
println("updated energy = ", energy(lattice))
