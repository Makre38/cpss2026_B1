module CPSS2026B1

using Random

export read_lattice_size, generate_lattice, generate_lattice_from_file

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

end
