using Test

include(joinpath(@__DIR__, "..", "src", "CPSS2026Ising.jl"))

using .CPSS2026Ising

@testset "parameter parsing" begin
    mktemp() do path, io
        write(io, "L = 4\nT = 2.5\nJ = 1.0\nn_therm = 100\nn_measure = 200\nsnapshot_interval = 10\n")
        close(io)

        params = read_parameters(path)
        @test params.L == 4
        @test params.temperature == 2.5
        @test params.coupling == 1.0
        @test params.n_therm == 100
        @test params.n_measure == 200
        @test params.snapshot_interval == 10
    end
end

@testset "lattice generation" begin
    lattice = generate_lattice(8)
    @test size(lattice) == (8, 8)
    @test all(spin -> spin in (-1, 1), lattice)
end

@testset "observables" begin
    lattice = Int8[1 1; 1 1]
    @test magnetization(lattice) == 4
    @test energy(lattice; J = 1.0) == -8.0

    checkerboard = Int8[1 -1; -1 1]
    @test magnetization(checkerboard) == 0
    @test energy(checkerboard; J = 1.0) == 8.0
end

@testset "flip energy consistency" begin
    lattice = Int8[1 1; 1 -1]
    row, col = 2, 2
    before = energy(lattice; J = 1.0)
    delta = CPSS2026Ising.energy_change_for_flip(lattice, row, col; J = 1.0)
    lattice[row, col] = -lattice[row, col]
    after = energy(lattice; J = 1.0)
    @test after - before == delta
end

@testset "metropolis invariants" begin
    lattice = Int8[1 1; 1 1]
    accepted = metropolis_sweep!(lattice; β = 0.4, J = 1.0)
    @test 0 <= accepted <= length(lattice)
    @test size(lattice) == (2, 2)
    @test all(spin -> spin in (-1, 1), lattice)
end
