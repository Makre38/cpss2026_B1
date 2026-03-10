include("src/CPSS2026B1.jl")

using .CPSS2026B1

const PARAMS_FILE = "params.dat"

function main()
    params = read_parameters(PARAMS_FILE)
    beta = 1 / params.temperature
    lattice = generate_lattice(params.L)

    println("Simulation parameters:")
    println("  L = ", params.L)
    println("  T = ", params.temperature)
    println("  beta = ", beta)
    println("  J = ", params.coupling)
    println("  n_therm = ", params.n_therm)
    println("  n_measure = ", params.n_measure)
    println()
    println("Initial observables:")
    println("  magnetization = ", magnetization(lattice))
    println("  energy = ", energy(lattice; J = params.coupling))
end

main()
