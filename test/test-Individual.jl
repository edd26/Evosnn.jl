using DrWatson, Test
@quickactivate "Evosnn"


@testset "testing Individual structure" begin

    # "Individual.jl" |> srcdir |> include

    ind = Individual()
    @test ind.fitness == 0.0
    @test ind.rewardn == 0.0
    @test ind.penaltyn == 0.0
    @test ind.fdr == 0.0
    @test ind.precision == 0.0
    @test ind.reward == 0.0
    @test ind.penalty == 0.0
    @test ind.absWeightSum == 0.0
    @test ind.rank == 0
    @test ind.totalCorrPatterns == 0
    @test isempty(ind.indMatrix)
    @test isempty(ind.inputNeurons)
    @test isempty(ind.interNeurons)
    @test isempty(ind.outputNeurons)
    @test ind.noOfInputs == 0
    @test ind.noOfinterNeurons == 0
    @test ind.noOfOutputNeurons == 0
    @test ind.noOfNodesInNetwork == 0
    @test isempty(ind.gaussNoiseVector)
    @test isempty(ind.missIdentifiedPatterns)
end

@testset "testing Individual methods" begin

    # "Individual.jl" |> srcdir |> include

    @testset "sumOfConnectionWeights Tests" begin
        # Create a sample Individual for testing
        ind_matrix = [
            0.0 0.2 0.0 0.5;
            0.1 0.0 0.0 0.4;
            0.0 0.3 0.0 0.0
        ]

        ind = Individual(;
            indMatrix=ind_matrix,
        )

        # Test the result
        @test sumOfConnectionWeights(ind) == 1.5
    end
end
