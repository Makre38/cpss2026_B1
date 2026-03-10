module CPSS2026B1

using Random

export AbstractBoundaryCondition,
    PeriodicBoundary,
    read_lattice_size,
    generate_lattice,
    generate_lattice_from_file,
    local_field,
    energy_change_for_flip,
    magnetization,
    energy

abstract type AbstractBoundaryCondition end

struct PeriodicBoundary <: AbstractBoundaryCondition end

function read_lattice_size(path::AbstractString)
    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "#") && continue

        if occursin("=", line)
            key, value = strip.(split(line, "=", limit = 2))
            key == "L" || continue
            return parse(Int, value)
        end

        return parse(Int, line)
    end

    error("Could not find lattice size L in $(path).")
end

function generate_lattice(L::Integer; rng::AbstractRNG = Random.default_rng())
    L > 0 || throw(ArgumentError("L must be positive."))
    spins = Int8[-1, 1]
    return rand(rng, spins, L, L)
end

function generate_lattice_from_file(path::AbstractString; rng::AbstractRNG = Random.default_rng())
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
