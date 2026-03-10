include("src/CPSS2026B1.jl")

using .CPSS2026B1

const PARAMS_FILE = "params.dat"
const OUTPUT_DIR = "output"
const SNAPSHOT_DIR = joinpath(OUTPUT_DIR, "snapshots")
const OBSERVABLES_FILE = joinpath(OUTPUT_DIR, "observables.dat")

function write_lattice_snapshot(path::AbstractString, lattice)
    open(path, "w") do io
        for row in eachrow(lattice)
            println(io, join(row, " "))
        end
    end
end

function initialize_observables_file(path::AbstractString)
    open(path, "w") do io
        println(io, "# sweep accepted magnetization energy")
    end
end

function append_observable(path::AbstractString, sweep::Integer, accepted::Integer, magnetization_value, energy_value)
    open(path, "a") do io
        println(io, sweep, " ", accepted, " ", magnetization_value, " ", energy_value)
    end
end

function main()
    params = read_parameters(PARAMS_FILE)
    beta = 1 / params.temperature
    lattice = generate_lattice(params.L)
    mkpath(SNAPSHOT_DIR)
    initialize_observables_file(OBSERVABLES_FILE)

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
    write_lattice_snapshot(joinpath(SNAPSHOT_DIR, "sweep_0000.dat"), lattice)

    println()
    println("Thermalization:")
    for sweep in 1:params.n_therm
        accepted = metropolis_sweep!(lattice; β = beta, J = params.coupling)
        println("  sweep ", sweep, ": accepted = ", accepted)
    end

    println()
    println("Measurements:")
    for sweep in 1:params.n_measure
        accepted = metropolis_sweep!(lattice; β = beta, J = params.coupling)
        m = magnetization(lattice)
        e = energy(lattice; J = params.coupling)
        append_observable(OBSERVABLES_FILE, sweep, accepted, m, e)
        write_lattice_snapshot(joinpath(SNAPSHOT_DIR, "sweep_" * lpad(string(sweep), 4, '0') * ".dat"), lattice)
    end

    println("  observables written to ", OBSERVABLES_FILE)
    println("  snapshots written to ", SNAPSHOT_DIR)
end

main()
