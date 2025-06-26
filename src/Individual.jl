module IndividualModule

export Individual,
    initializeMatrix,
    getIndividualMatrix,
    printIndividualMatrix,
    outputNetworkActivity,
    readIndividualMatrix,
    makeRandomIndividual,
    makeIndividualWithFixedInputOutputConnections,
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
    setInput!,
    setGap,
    resetIndividual,
    sumOfConnectionWeights,
    copy_individual

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


function copy_individual(ind::Individual)::Individual
    return Individual(
        fitness=ind.fitness,
        reward=ind.reward,
        rewardn=ind.rewardn,
        penalty=ind.penalty,
        penaltyn=ind.penaltyn,
        fdr=ind.fdr,
        precision=ind.precision,
        rank=ind.rank,
        totalCorrPatterns=ind.totalCorrPatterns,
        absWeightSum=ind.absWeightSum,
        noOfInputs=ind.noOfInputs,
        noOfinterNeurons=ind.noOfinterNeurons,
        noOfOutputNeurons=ind.noOfOutputNeurons,
        noOfNodesInNetwork=ind.noOfNodesInNetwork,
        inputNeurons=deepcopy(ind.inputNeurons),
        interNeurons=deepcopy(ind.interNeurons),
        outputNeurons=deepcopy(ind.outputNeurons),
        gaussNoiseVector=copy(ind.gaussNoiseVector),
        indMatrix=copy(ind.indMatrix),
        missIdentifiedPatterns=copy(ind.missIdentifiedPatterns)
    )
end



function outputNetworkActivity(ind::Individual, outputfile::String)
    open(outputfile, "w") do ofs
        alphabet = "ABCDEFGHIJKLMNOPQRSTUWXYZ"
        signal_sequence = alphabet[1:length(ind.inputNeurons)]

        header = "#No\t Out"
        values = "-1\t -70"

        for c in signal_sequence
            header *= "\t $c"
            values *= "\t -100"
        end

        for i in 1:length(signal_sequence)
            header *= "\t N$(i-1)"
            values *= "\t -70"
        end

        println(ofs, header)
        println(ofs, values)

        println(length(ind.outputNeurons[1].voltageBuffer))

        for i in 1:length(ind.outputNeurons[1].voltageBuffer)
            line = string(i - 1, "\t", ind.outputNeurons[1].voltageBuffer[i], "\t")

            for inpN in 1:length(ind.inputNeurons)
                voltage = ind.inputNeurons[inpN].voltageBuffer[i]
                value = voltage == 0 ? -(inpN - 1) * 10 : -100
                line *= string(value, "\t")
            end

            for intN in 1:length(ind.interNeurons)
                line *= string(ind.interNeurons[intN].voltageBuffer[i], "\t")
            end

            println(ofs, line)
        end
    end
end


function readIndividualMatrix(fname::String, noOfNodesInNetwork::Int)
    indMatrix_components = Vector{Vector{Float64}}()
    if !isfile(fname)
        @info "Ups..."
        ("File \"$(fname)\" does not exists!" |> ErrorException |> throw)
    end

    open(fname, "r") do file
        for line in eachline(file)
            row = Float64[]
            for value in split(line, "\t")
                if value != "" && value != "\b"
                    push!(row, parse(Float64, value))
                end
            end
            if !isempty(row)
                push!(indMatrix_components, row)
            end
        end
    end

    all([length(v) for v in indMatrix_components] .== length(indMatrix_components)
    ) || "Not all of the components in the loaded matrix were of the same size" |> DimensionMismatch |> throw

    # Reshape loaded vectors into the matrix
    indMatrix = zeros(noOfNodesInNetwork, noOfNodesInNetwork)
    for (r_index, r) in enumerate(indMatrix_components)
        for (col_index, v) in enumerate(r)
            indMatrix[r_index, col_index] = v
        end
    end

    return indMatrix
end

function makeRandomIndividual(ind::Individual, ntype::Neuron, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int)
    @error "Not implemented!"
end

# TODO add ! to the function name to indicate changing of the argument
function makeIndividualWithFixedInputOutputConnections(ind::Individual, ntype::NeuronType, noInputs::Int, nointerNeurons::Int, noOutputs::Int, ind_matrix::Matrix{T}) where {T<:Number}# file_to_load::String)
    noOfInputs = noInputs
    # ind.noOfinterNeurons = nointerNeurons
    # ind.noOfOutputNeurons = noOutputs

    input_neuron_voltage = 0.0
    inter_neuron_voltage = -70.0
    output_neuron_voltage = -70.0

    ind.inputNeurons = [Neuron(ntype; voltage=input_neuron_voltage) for i in 1:noOfInputs]
    ind.interNeurons = [Neuron(ntype; voltage=inter_neuron_voltage) for i in 1:nointerNeurons]
    ind.outputNeurons = [Neuron(ntype; voltage=output_neuron_voltage) for i in 1:noOutputs]

    # ind.noOfNodesInNetwork = noInputs + nointerNeurons + noOutputs
    # initializeMatrix(noOfNodesInNetwork)
    ind.indMatrix = ind_matrix

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

function replicate(ind::Individual, params::Parameters)
    new_ind = copy_individual(ind)

    for i in 1:ind.noOfInputs
        for j in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0.0, 1.0)
                if params.mutationProb > randValue
                    new_ind.indMatrix[i, j] += getRandomValue(-1.0 * params.mutationStength, params.mutationStength)
                end
            end
        end
    end

    for i in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
        for j in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0.0, 1.0)
                if params.mutationProb > randValue
                    new_ind.indMatrix[i, j] += getRandomValue(-1.0 * params.mutationStength, params.mutationStength)
                end
            end
        end
    end

    for i in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
        for j in (ind.noOfInputs+ind.noOfinterNeurons+1):ind.noOfNodesInNetwork
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0.0, 1.0)
                if params.mutationProb > randValue
                    new_ind.indMatrix[i, j] += getRandomValue(-1.0 * params.mutationStength, params.mutationStength)
                end
            end
        end
    end
    return new_ind
end




function replicate_might_not_work!(ind::Individual, params::Parameters)
    randValue = 0.0

    # input connects inter-neurons
    for i in 1:ind.noOfInputs
        for j in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0, 1)

                # Mutate Weights
                if params.mutationProb > randValue
                    ind.indMatrix[i, j] += getRandomValue(-params.mutationStength, params.mutationStength)
                end
            end
        end
    end

    # connections among inter-neurons
    for i in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
        for j in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0, 1)

                # Mutate Weights
                if params.mutationProb > randValue
                    ind.indMatrix[i, j] += getRandomValue(-params.mutationStength, params.mutationStength)
                end
            end
        end
    end

    # interneurons to output
    for i in (ind.noOfInputs+1):(ind.noOfInputs+ind.noOfinterNeurons)
        for j in (ind.noOfInputs+ind.noOfinterNeurons+1):ind.noOfNodesInNetwork
            if ind.indMatrix[i, j] != 0.0
                randValue = getRandomValue(0, 1)

                # Mutate Weight
                if params.mutationProb > randValue
                    ind.indMatrix[i, j] += getRandomValue(-params.mutationStength, params.mutationStength)
                end
            end
        end
    end
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

function networkStep!(ind::Individual, stepNo::Int64, params::Parameters)
    for n in 1:length(ind.interNeurons)
        currNeuron = ind.interNeurons[n]
        if hasSpike(currNeuron.cUnit, currNeuron.voltage)
            currNeuron.cUnit.hasSpiked = true

            currNeuron.voltage = resetVoltage!(currNeuron.cUnit)

            if params.neuronalType == ADX
                currNeuron.adaptation =
                    resetAdaptation(currNeuron.cUnit, currNeuron.adaptation)
            end
            currNeuron.refrectoryPeriod = currNeuron.cUnit.absRefractoryTime
            if params.writeNetworkActivity
                push!(currNeuron.voltageBuffer, currNeuron.voltage)
                push!(currNeuron.spikeBitmap, true)
            end
            currNeuron.cUnit.spikeCount += 1
        else
            if currNeuron.refrectoryPeriod <= 0
                tempVoltage = currNeuron.voltage

                currNeuron.voltage = updateVoltage(currNeuron.cUnit, currNeuron.voltage, currNeuron.exConductance, currNeuron.inConductance, currNeuron.adaptation, params.timeStep)
                if params.neuronalType == ADX
                    currNeuron.adaptation = updateAdaptation(currNeuron.cUnit, currNeuron.adaptation, tempVoltage, params.timeStep)
                end
                if params.gaussianNoiseOnVoltage
                    if isempty(gaussNoiseVector)
                        gaussNoiseVector = getGaussianValueWithGivenMeanAndSD(params.gMean, params.gStdDev, params.noiseVectorSize)
                    else
                        currNeuron.voltage += gaussNoiseVector[rand(1:length(gaussNoiseVector))]
                    end
                end
                if hasSpike(currNeuron.cUnit, currNeuron.voltage)
                    currNeuron.voltage = currNeuron.cUnit.Vspike
                end
            else
                currNeuron.refrectoryPeriod -= 1
            end
            if params.writeNetworkActivity
                push!(currNeuron.voltageBuffer, currNeuron.voltage)
                push!(currNeuron.spikeBitmap, false)
            end
        end
        currNeuron.exConductance = updateExcitatoryCond(currNeuron.cUnit, currNeuron.exConductance, params.timeStep)
        currNeuron.inConductance = updateInhibitoryCond(currNeuron.cUnit, currNeuron.inConductance, params.timeStep)

        for interConn in 1:length(ind.interNeurons)
            index1 = ind.noOfInputs + interConn
            index2 = ind.noOfInputs + n
            connWeight = ind.indMatrix[index1, index2]

            if hasSpike(ind.interNeurons[interConn].cUnit, ind.interNeurons[interConn].voltage)
                if connWeight > 0
                    currNeuron.exConductance += params.ge_gain * connWeight
                elseif connWeight < 0
                    currNeuron.inConductance += params.gi_gain * (-connWeight)
                end
            end
        end

        for inputConn in 1:length(ind.inputNeurons)
            index1 = inputConn
            index2 = ind.noOfInputs + n
            connWeight = ind.indMatrix[index1, index2]
            if hasSpike(ind.inputNeurons[inputConn].cUnit, ind.inputNeurons[inputConn].voltage)
                if connWeight > 0
                    currNeuron.exConductance += params.ge_gain * connWeight
                elseif connWeight < 0
                    currNeuron.inConductance += params.gi_gain * (-connWeight)
                end
            end
        end
    end
    activateOutput(ind, stepNo, params)
end

function activateOutput(ind::Individual, stepNo::Int64, params::Parameters)
    for out_n in 1:length(ind.outputNeurons)
        outputNeuron = ind.outputNeurons[out_n]
        if hasSpike(outputNeuron.cUnit, outputNeuron.voltage)
            # outputNeuron.voltage = 
            # resetVoltage!(outputNeuron.cUnit)
            outputNeuron.voltage = resetVoltage!(outputNeuron.cUnit)

            if params.neuronalType == ADX
                outputNeuron.adaptation = resetAdaptation(outputNeuron.cUnit, outputNeuron.adaptation)
            end
            outputNeuron.refrectoryPeriod = outputNeuron.cUnit.absRefractoryTime
            push!(outputNeuron.spikeBitmap, true)
            if params.writeNetworkActivity
                push!(outputNeuron.voltageBuffer, outputNeuron.voltage)
            end
            outputNeuron.cUnit.spikeCount += 1
        else
            if outputNeuron.refrectoryPeriod <= 0
                tempVoltage = outputNeuron.voltage
                outputNeuron.voltage = updateVoltage(outputNeuron.cUnit, outputNeuron.voltage, outputNeuron.exConductance, outputNeuron.inConductance, outputNeuron.adaptation, params.timeStep)
                if params.neuronalType == ADX
                    outputNeuron.adaptation = updateAdaptation(outputNeuron.cUnit, outputNeuron.adaptation, tempVoltage, params.timeStep)
                end
                if params.gaussianNoiseOnVoltage
                    if isempty(gaussNoiseVector)
                        gaussNoiseVector = getGaussianValueWithGivenMeanAndSD(params.gMean, params.gStdDev, params.noiseVectorSize)
                    else
                        outputNeuron.voltage += gaussNoiseVector[rand(1:length(gaussNoiseVector))]
                    end
                end
                if hasSpike(outputNeuron.cUnit, outputNeuron.voltage)
                    outputNeuron.voltage = outputNeuron.cUnit.Vspike
                end
            else
                outputNeuron.refrectoryPeriod -= 1
            end
            if params.writeNetworkActivity
                push!(outputNeuron.voltageBuffer, outputNeuron.voltage)
            end
            push!(outputNeuron.spikeBitmap, false)
        end
        outputNeuron.exConductance = updateExcitatoryCond(outputNeuron.cUnit, outputNeuron.exConductance, params.timeStep)
        outputNeuron.inConductance = updateInhibitoryCond(outputNeuron.cUnit, outputNeuron.inConductance, params.timeStep)

        for outConn in 1:length(ind.interNeurons)
            index1 = ind.noOfInputs + outConn
            index2 = ind.noOfInputs + ind.noOfinterNeurons + out_n
            connWeight = ind.indMatrix[index1, index2]
            if hasSpike(ind.interNeurons[outConn].cUnit, ind.interNeurons[outConn].voltage)
                if connWeight > 0
                    outputNeuron.exConductance += params.ge_gain * connWeight
                elseif connWeight < 0
                    outputNeuron.inConductance += params.gi_gain * (-connWeight)
                end
            end
        end
    end
end

function setInput!(ind::Individual, inputSignal::Char, index::Int64, writeNetworkActivity::Bool)
    # Initialize all input neurons with default voltage
    for neuron in ind.inputNeurons
        neuron.voltage = -100
        if writeNetworkActivity
            push!(neuron.voltageBuffer, -100)
        end
    end

    # Map input signals to neuron indices
    signalToNeuronIndex = Dict('A' => 1, 'B' => 2, 'C' => 3, 'D' => 4,
        'E' => 5, 'F' => 6, 'G' => 7, 'H' => 8,
        'I' => 9, 'J' => 10)

    # Set voltage for the corresponding neuron if the signal is valid
    if haskey(signalToNeuronIndex, inputSignal)
        neuronIndex = signalToNeuronIndex[inputSignal]
        # @info "For $(inputSignal) setting neuron index $(neuronIndex)"
        ind.inputNeurons[neuronIndex].voltage = 0.0
        if writeNetworkActivity
            ind.inputNeurons[neuronIndex].voltageBuffer[index] = 0.0
        end
        # else
        #     "Given input singnal char was not recognised: $(inputSignal)" |> Exception |> throw
    end
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
