module CPSS2026Ising

export AbstractBoundaryCondition,
    PeriodicBoundary,
    SimulationParameters,
    read_parameters,
    read_lattice_size,
    generate_lattice,
    generate_lattice_from_file,
    metropolis_step!,
    metropolis_sweep!,
    magnetization,
    energy

abstract type AbstractBoundaryCondition end

struct PeriodicBoundary <: AbstractBoundaryCondition end

struct SimulationParameters
    L::Int
    temperature::Float64
    coupling::Float64
    n_therm::Int
    n_measure::Int
    snapshot_interval::Int
end

function parse_parameter_file(path::AbstractString)
    parameters = Dict{String, String}()

    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "#") && continue
        occursin("=", line) || throw(ArgumentError("Invalid parameter line: $(raw_line)"))

        key, value = strip.(split(line, "=", limit = 2))
        parameters[key] = value
    end

    return parameters
end

function read_parameters(path::AbstractString)
    parameters = parse_parameter_file(path)
    haskey(parameters, "L") || error("Could not find lattice size L in $(path).")
    haskey(parameters, "T") || error("Could not find temperature T in $(path).")
    haskey(parameters, "J") || error("Could not find coupling J in $(path).")
    haskey(parameters, "n_therm") || error("Could not find n_therm in $(path).")
    haskey(parameters, "n_measure") || error("Could not find n_measure in $(path).")
    haskey(parameters, "snapshot_interval") || error("Could not find snapshot_interval in $(path).")

    return SimulationParameters(
        parse(Int, parameters["L"]),
        parse(Float64, parameters["T"]),
        parse(Float64, parameters["J"]),
        parse(Int, parameters["n_therm"]),
        parse(Int, parameters["n_measure"]),
        parse(Int, parameters["snapshot_interval"]),
    )
end

function read_lattice_size(path::AbstractString)
    return read_parameters(path).L
end

function generate_lattice(L::Integer; rng = nothing)
    L > 0 || throw(ArgumentError("L must be positive."))
    spins = Int8[-1, 1]
    return isnothing(rng) ? rand(spins, L, L) : rand(rng, spins, L, L)
end

function generate_lattice_from_file(path::AbstractString; rng = nothing)
    L = read_lattice_size(path)
    return generate_lattice(L; rng = rng)
end

function magnetization(lattice::AbstractMatrix{<:Integer})
    return sum(lattice)
end

function neighbor_index(index::Integer, size::Integer, ::PeriodicBoundary)
    return mod1(index, size)
end

function local_field(
    lattice::AbstractMatrix{<:Integer},
    row::Integer,
    col::Integer;
    boundary::AbstractBoundaryCondition = PeriodicBoundary(),
)
    nrows, ncols = size(lattice)
    nrows == ncols || throw(ArgumentError("lattice must be square."))
    1 <= row <= nrows || throw(BoundsError(lattice, (row, col)))
    1 <= col <= ncols || throw(BoundsError(lattice, (row, col)))

    row_minus = neighbor_index(row - 1, nrows, boundary)
    row_plus = neighbor_index(row + 1, nrows, boundary)
    col_minus = neighbor_index(col - 1, ncols, boundary)
    col_plus = neighbor_index(col + 1, ncols, boundary)

    return lattice[row_minus, col] + lattice[row_plus, col] + lattice[row, col_minus] + lattice[row, col_plus]
end

function energy_change_for_flip(
    lattice::AbstractMatrix{<:Integer},
    row::Integer,
    col::Integer;
    J::Real = 1,
    boundary::AbstractBoundaryCondition = PeriodicBoundary(),
)
    spin = lattice[row, col]
    return 2 * J * spin * local_field(lattice, row, col; boundary = boundary)
end

function metropolis_step!(
    lattice::AbstractMatrix{<:Integer};
    β::Real,
    J::Real = 1,
    boundary::AbstractBoundaryCondition = PeriodicBoundary(),
    rng = nothing,
)
    nrows, ncols = size(lattice)
    nrows == ncols || throw(ArgumentError("lattice must be square."))

    row = isnothing(rng) ? rand(1:nrows) : rand(rng, 1:nrows)
    col = isnothing(rng) ? rand(1:ncols) : rand(rng, 1:ncols)
    ΔE = energy_change_for_flip(lattice, row, col; J = J, boundary = boundary)

    threshold = isnothing(rng) ? rand() : rand(rng)
    accepted = ΔE <= 0 || threshold < exp(-β * ΔE)
    if accepted
        lattice[row, col] = -lattice[row, col]
    end

    return accepted, row, col, ΔE
end

function metropolis_sweep!(
    lattice::AbstractMatrix{<:Integer};
    β::Real,
    J::Real = 1,
    boundary::AbstractBoundaryCondition = PeriodicBoundary(),
    rng = nothing,
)
    accepted = 0
    trials = length(lattice)

    for _ in 1:trials
        move_accepted, _, _, _ = metropolis_step!(lattice; β = β, J = J, boundary = boundary, rng = rng)
        accepted += move_accepted
    end

    return accepted
end

function energy(lattice::AbstractMatrix{<:Integer}; J::Real = 1, boundary::AbstractBoundaryCondition = PeriodicBoundary())
    nrows, ncols = size(lattice)
    nrows == ncols || throw(ArgumentError("lattice must be square."))

    total = zero(promote_type(eltype(lattice), typeof(J)))

    for row in 1:nrows
        row_plus = neighbor_index(row + 1, nrows, boundary)
        for col in 1:ncols
            col_plus = neighbor_index(col + 1, ncols, boundary)
            spin = lattice[row, col]
            total -= J * spin * lattice[row_plus, col]
            total -= J * spin * lattice[row, col_plus]
        end
    end

    return total
end

end
