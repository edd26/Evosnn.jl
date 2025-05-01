using DrWatson
@quickactivate "Evosnn"

using Test
using Evosnn

# Run test suite
println("Starting tests")
ti = time()

@testset "Structures construction and class methods" begin

    "test/test-Parameters.jl" |> projectdir |> include
    "test/test-Individual.jl" |> projectdir |> include
    "test/test-ComputationalUnit.jl" |> projectdir |> include

    "test/test-Neuron.jl" |> projectdir |> include
end

@testset "Methods accross structures" begin
    "test/test-utils.jl" |> projectdir |> include
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
