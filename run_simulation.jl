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

    println()
    println("Thermalization:")
    for sweep in 1:params.n_therm
        accepted = metropolis_sweep!(lattice; β = beta, J = params.coupling)
        println("  sweep ", sweep, ": accepted = ", accepted)
    end

    println()
    println("Measurements:")
    println("sweep,accepted,magnetization,energy")
    for sweep in 1:params.n_measure
        accepted = metropolis_sweep!(lattice; β = beta, J = params.coupling)
        m = magnetization(lattice)
        e = energy(lattice; J = params.coupling)
        println(sweep, ",", accepted, ",", m, ",", e)
    end
end

main()
