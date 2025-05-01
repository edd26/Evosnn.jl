module IndividualModule

export Individual,
    initializeMatrix,
    getIndividualMatrix,
    printIndividualMatrix,
    outputNetworkActivity, readIndividualMatrix,
    makeRandomIndividual, makeIndividualWithFixedInputOutputConnections,
    makeIndividualWithFixedInputOutputConnections,
    loadEvolvedTop,
    loadIndividualNetwork,
    makeIndividualWithFixedLoops,
    replicate,
    replicateWithGaussianNoise,
    replicateInputconnection_WithGaussianNoise,
    replicateInput_and_Switch_connections,
    replicateinterConnections,
    randomizeinterconnections,
    replicateExceptLoops,
    replicateByOnlyReducingWeights,
    deleteWeight,
    mutateWeight,
    mutateSign,
    deleteIndMatrix,
    networkStep,
    activateOutput,
    setInput,
    setGap,
    resetIndividual,
    sumOfConnectionWeights

using ..NeuronModule: Neuron
using ..MyTypes: NeuronType

# TODO verify correctness with C++
struct Individual
    fitness::Float64
    rewardn::Float64
    penaltyn::Float64
    fdr::Float64
    precision::Float64
    reward::Float64
    penalty::Float64
    absWeightSum::Float64
    rank::Int
    totalCorrPatterns::Int
    indMatrix::Matrix{Float64}
    inputNeurons::Vector{Neuron}
    interNeurons::Vector{Neuron}
    outputNeurons::Vector{Neuron}
    noOfInputs::Int
    noOfinterNeurons::Int
    noOfOutputNeurons::Int
    noOfNodesInNetwork::Int
    gaussNoiseVector::Vector{Float64}
    missIdentifiedPatterns::Vector{String}

    function Individual(;
        fitness::Float64=0.0,
        rewardn::Float64=0.0,
        penaltyn::Float64=0.0,
        fdr::Float64=0.0,
        precision::Float64=0.0,
        reward::Float64=0.0,
        penalty::Float64=0.0,
        absWeightSum::Float64=0.0,
        rank::Int=0,
        totalCorrPatterns::Int=0,
        indMatrix::Matrix{Float64}=zeros(Float64, 0, 0),
        inputNeurons::Vector{Neuron}=Vector{Neuron}(),
        interNeurons::Vector{Neuron}=Vector{Neuron}(),
        outputNeurons::Vector{Neuron}=Vector{Neuron}(),
        noOfInputs::Int=0,
        noOfinterNeurons::Int=0,
        noOfOutputNeurons::Int=0,
        noOfNodesInNetwork::Int=0,
        gaussNoiseVector::Vector{Float64}=Vector{Float64}(),
        missIdentifiedPatterns::Vector{String}=Vector{String}()
    )
        new(
            fitness, rewardn, penaltyn, fdr, precision,
            reward, penalty, absWeightSum, rank, totalCorrPatterns,
            indMatrix, inputNeurons, interNeurons, outputNeurons,
            noOfInputs, noOfinterNeurons, noOfOutputNeurons, noOfNodesInNetwork,
            gaussNoiseVector, missIdentifiedPatterns
        )
    end
end

function initializeMatrix(ind::Individual, size::Int)
end

function getIndividualMatrix(ind::Individual)
end

function printIndividualMatrix(ind::Individual, gNo::Int, irun::Int)
end


function outputNetworkActivity(ind::Individual, output_file::String)
end

function readIndividualMatrix(ind::Individual, fname::String)
end

function makeRandomIndividual(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
end

function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
end

function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int, file_to_load::String)
end

function loadEvolvedTop(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
end

function loadIndividualNetwork(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int, name::String)
end

function makeIndividualWithFixedLoops(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
end

function replicate(ind::Individual)
end

function replicateWithGaussianNoise(ind::Individual)
end

function replicateInputconnection_WithGaussianNoise(ind::Individual)
end

function replicateInput_and_Switch_connections(ind::Individual)
end

function replicateinterConnections(ind::Individual)
end

function randomizeinterconnections(ind::Individual)
end

function replicateExceptLoops(ind::Individual)
end

function replicateByOnlyReducingWeights(ind::Individual)
end

function deleteWeight(ind::Individual)
end

function mutateWeight(ind::Individual)
end

function mutateSign(ind::Individual)
end


function deleteIndMatrix(ind::Individual)
end

function networkStep(ind::Individual, stepNo::Int64)
end

function activateOutput(ind::Individual, stepNo::Int64)
end

function setInput(ind::Individual, inputSignal::Char, index::Int64)
end

function setGap(ind::Individual)
end

function resetIndividual(ind::Individual)
end


function sumOfConnectionWeights(individual::Individual)
    return individual.indMatrix .|> abs |> sum

    # absWeightSum::Float64 = 0.0
    # for i in 1:size(individual.indMatrix, 1)
    #     for j in (individual.noOfInputs+1):size(individual.indMatrix, 2)
    #         if individual.indMatrix[i][j] != 0.0
    #             absWeightSum += abs()
    #         end
    #     end
    # end
    # return absWeightSum
end

end # module
