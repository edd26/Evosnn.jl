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
using ..CompUnit: hasSpike, updateVoltage, updateAdaptation, updateExcitatoryCond, updateInhibitoryCond, resetVoltage!, resetAdaptation
using ..ParamsModule: Parameters
using ..MyTypes
using ..UtilityFunctions: getRandomValue

# TODO verify correctness with C++
mutable struct Individual
    fitness::Float64
    reward::Float64
    rewardn::Float64

    penalty::Float64
    penaltyn::Float64
    fdr::Float64
    precision::Float64

    rank::Int
    totalCorrPatterns::Int

    absWeightSum::Float64

    noOfInputs::Int
    noOfinterNeurons::Int
    noOfOutputNeurons::Int
    noOfNodesInNetwork::Int

    inputNeurons::Vector{Neuron}
    interNeurons::Vector{Neuron}
    outputNeurons::Vector{Neuron}

    gaussNoiseVector::Vector{Float64}

    indMatrix::Matrix{Float64}
    missIdentifiedPatterns::Vector{String}

    function Individual(;
        fitness::Float64=0.0,
        reward::Float64=0.0,
        rewardn::Float64=0.0,
        penalty::Float64=0.0,
        penaltyn::Float64=0.0,
        fdr::Float64=0.0,
        precision::Float64=0.0,
        rank::Int=0,
        totalCorrPatterns::Int=0,
        absWeightSum::Float64=0.0,
        noOfInputs::Int,
        noOfinterNeurons::Int,
        noOfOutputNeurons::Int,
        noOfNodesInNetwork::Int=0,
        inputNeurons::Vector{Neuron}=Vector{Neuron}(),
        interNeurons::Vector{Neuron}=Vector{Neuron}(),
        outputNeurons::Vector{Neuron}=Vector{Neuron}(),
        gaussNoiseVector::Vector{Float64}=Vector{Float64}(),
        indMatrix::Matrix{Float64}=zeros(Float64, 0, 0),
        missIdentifiedPatterns::Vector{String}=Vector{String}()
    )
        noOfNodesInNetwork = if noOfNodesInNetwork == 0
            noOfInputs + noOfinterNeurons + noOfOutputNeurons
        else
            noOfNodesInNetwork
        end
        new(
            fitness,
            reward,
            rewardn,
            penalty,
            penaltyn,
            fdr,
            precision,
            rank,
            totalCorrPatterns,
            absWeightSum,
            noOfInputs,
            noOfinterNeurons,
            noOfOutputNeurons,
            noOfNodesInNetwork,
            inputNeurons,
            interNeurons,
            outputNeurons,
            gaussNoiseVector,
            indMatrix,
            missIdentifiedPatterns,
        )
    end
end

function initializeMatrix(ind::Individual, size::Int)
    @error "Not implemented!"
end

function getIndividualMatrix(ind::Individual)
    @error "Not implemented!"
end

function printIndividualMatrix(ind::Individual, gNo::Int, irun::Int)
    @error "Not implemented!"
end


function outputNetworkActivity(ind::Individual, output_file::String)
    @error "Not implemented!"
end

function readIndividualMatrix(ind::Individual, fname::String)
end

function makeRandomIndividual(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int, file_to_load::String)
end

function loadEvolvedTop(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

function loadIndividualNetwork(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int, name::String)
    @error "Not implemented!"
end

function makeIndividualWithFixedLoops(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

function replicate(ind::Individual)
end

function replicateWithGaussianNoise(ind::Individual)
    @error "Not implemented!"
end

function replicateInputconnection_WithGaussianNoise(ind::Individual)
    @error "Not implemented!"
end

function replicateInput_and_Switch_connections(ind::Individual)
    @error "Not implemented!"
end

function replicateinterConnections(ind::Individual)
    @error "Not implemented!"
end

function randomizeinterconnections(ind::Individual)
    @error "Not implemented!"
end

function replicateExceptLoops(ind::Individual)
    @error "Not implemented!"
end

function replicateByOnlyReducingWeights(ind::Individual)
    @error "Not implemented!"
end

function deleteWeight(ind::Individual)
    @error "Not implemented!"
end

function mutateWeight(ind::Individual)
    @error "Not implemented!"
end

function mutateSign(ind::Individual)
    @error "Not implemented!"
end


function deleteIndMatrix(ind::Individual)
    @error "Not implemented!"
end

function networkStep(ind::Individual, stepNo::Int64)
end

function activateOutput(ind::Individual, stepNo::Int64)
end

function setInput(ind::Individual, inputSignal::Char, index::Int64)
end

function setGap(ind::Individual)
    @error "Not implemented!"
end

function resetIndividual(ind::Individual)
    @error "Not implemented!"
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
