include("src/CPSS2026B1.jl")

using .CPSS2026B1

const PARAMS_FILE = "params.dat"

lattice = generate_lattice_from_file(PARAMS_FILE)

println("L = ", size(lattice, 1))
println("lattice =")
for row in eachrow(lattice)
    println(join(row, " "))
end
