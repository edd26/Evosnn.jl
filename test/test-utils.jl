using DrWatson, Test
@quickactivate "Evosnn"


@testset "testing Individual and Parameters methods" begin
    # "Prameters.jl" |> srcdir |> include
    # "Individual.jl" |> srcdir |> include

    # Test function for printIndividualMatrix
    @testset "printIndividualMatrix Tests" begin
        # Define a mock parameters structure for testing
        # "Parameters.jl" |> srcdir |> include

        # Mock parameters instance
        params = Parameters()

        # Create a sample Individual for testing
        individual = Individual()

        # Mock directory and file name
        saving_dir = "./"
        orgi_file_name = "test_file"
        file_core_name = "ind_Matrix_s"

        # Call the function
        printIndividualMatrix(individual, params, 0, 0, saving_dir, orgi_file_name;
            file_core_name=file_core_name, file_extension=".txt",)

        # Check if the file was created and has content
        file_path = joinpath(saving_dir, "$(file_core_name)10_$(orgi_file_name).txt")
        @info pwd()
        @test isfile(file_path)

        # Read back the file and check content format
        content = read(file_path, String)
        @test occursin("Adj. matrix of the network", content)
        @test occursin("run:0	gen:0	fitness =", content)
        @test occursin("Total correct patterns in the sequence:", content)
        @test occursin("Sum of all:", content)

        # Clean up
        rm(file_path)
    end


    # TODO investigate Why there are two functions for output?
    # Test function for printIndividualMatrix
    # @testset "printIndividualMatrix Tests" begin
    #     # Create a sample Individual for testing
    #     individual = Individual()
    #
    #     # Mock parameters instance
    #     params = Parameters()
    #     # Mock directory and file name
    #     saving_dir = "./"
    #     orgi_file_name = "test_file"
    #     file_core_name = "ind_Matrix_s"
    #
    #     # Call the function
    #     printIndividualMatrix(individual, params, 0, 0;
    #         saving_dir=saving_dir, file_core_name=orgi_file_name,)
    #
    #     file_path = joinpath(saving_dir, "$(file_core_name)10_$(orgi_file_name).txt")
    #     @test isfile(file_path)
    #
    #     # Read back the file and check content format
    #     content = read(file_path, String)
    #     @test occursin("Adj. matrix of the network", content)
    #     @test occursin("run:0	gen:0	fitness =", content)
    #     @test occursin("Total corr =", content)
    #     @test occursin("identified corr =", content)
    #     @test occursin("wrong =", content)
    #     @test occursin("FDR =", content)
    #     @test occursin("Precision =", content)
    #     @test occursin("fitness:", content)
    #     @test occursin("reward:", content)
    #     @test occursin("penalty:", content)
    #     @test occursin("rewardn (a.k.a. identified corr):", content)
    #     @test occursin("penalty (a.k.a. wrong):", content)
    #
    #     # Clean up
    #     rm(file_path)
    # end



    @testset "removeLowWeights Tests" begin
        # Create a sample Individual for testing
        ind_matrix = [
            0.1 0.6 0.2;
            0.4 0.8 0.05;
            0.9 0.02 0.3
        ]

        individual = Individual(;
            indMatrix=ind_matrix,
        )

        # Mock parameters instance
        params = Parameters(;
            minConnectionWeight=0.05,
            maxConnectionWeight=0.8,
        )
        # params = Parameters(maxConnectionWeight=0.05)

        # Execute the function
        new_weights = removeLowWeights(individual, params)

        # Expected result
        expected_weights = [
            0.1 0.6 0.2;
            0.4 0.0 0.0;
            0.0 0.0 0.3
        ]

        # Test the result
        @test new_weights == expected_weights
    end
end
