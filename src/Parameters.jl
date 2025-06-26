module ParamsModule

export Parameters,
    get_nodes_in_network

using DrWatson: datadir

using ..MyTypes: NeuronType, LIF, ADX, IZH
using ..MyTypes: ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest, doFixedSequenceTest

struct Parameters
    # doEvolution::Bool
    # doAllSequencesTest::Bool
    executionMode::ExecutionMode
    getAlmostCorrectNWs::Bool
    writeNetworkActivity::Bool

    outputdir::String

    minConnectionWeight::Float64
    maxConnectionWeight::Float64
    timeStep::Float64
    neuronalType::NeuronType

    noOfInputs::Int
    noOfinterNeurons::Int
    noOfSignals::Int
    noOfOutputs::Int
    ge_gain::Float64
    gi_gain::Float64

    synapticDelay::Float64
    minWeightThreshold::Float64
    maxWeightRhreshold::Float64

    popSize::Int
    eliteCount::Int
    randomizeCount::Int
    randomizeEveryXGen::Int
    weightDeletionProb::Float64
    mutationProb::Float64
    signChangeProb::Float64
    mutationMean::Float64
    mutationStength::Float64
    maxGen::Int
    maxRuns::Int

    weightReductionProb::Float64
    weightReductionStrength::Float64

    noiseVectorSize::Int
    variationOnSignal::Int
    variationOnSilence::Int

    gaussianNoiseOnVoltage::Bool
    gMean::Float64
    gStdDev::Float64

    letterSize::Int
    silenctInterval::Int
    noOfLetters::Int
    reevaluateSeq::Int
    sequenceSize::Int

    function Parameters(;
        # doEvolution::Bool=false,
        # doAllSequencesTest::Bool=false,
        executionMode=doEvolution,
        getAlmostCorrectNWs::Bool=false,
        writeNetworkActivity::Bool=false,
        outputdir::String=nothing,
        minConnectionWeight::Float64=-10.0,
        maxConnectionWeight::Float64=10.0,
        timeStep::Float64=1.0,
        neuronalType::NeuronType=ADX,
        noOfInputs::Int=10,
        noOfinterNeurons::Int=10,
        noOfSignals::Int=10,
        noOfOutputs::Int=1,
        ge_gain::Float64=0.007,
        gi_gain::Float64=0.007,
        synapticDelay::Float64=4.0,
        minWeightThreshold::Float64=0.0,
        maxWeightRhreshold::Float64=0.0,
        popSize::Int=100,
        eliteCount::Int=5,
        randomizeCount::Int=5,
        randomizeEveryXGen::Int=10,
        weightDeletionProb::Float64=0.0,
        mutationProb::Float64=0.1,
        signChangeProb::Float64=0.0,
        mutationMean::Float64=0.5,
        mutationStength::Float64=0.7,
        maxGen::Int=10000,
        maxRuns::Int=10,
        weightReductionProb::Float64=0.2,
        weightReductionStrength::Float64=0.1,
        noiseVectorSize::Int=1000,
        variationOnSignal::Int=0,
        variationOnSilence::Int=0,
        gaussianNoiseOnVoltage::Bool=true,
        gMean::Float64=0.0,
        gStdDev::Float64=0.5,
        letterSize::Int=6,
        silenctInterval::Int=24,
        noOfLetters::Int=1500 * 4,
        reevaluateSeq::Int=10000
    )

        resolved_outputdir =
            if isnothing(outputdir)
                if (executionMode == doEvolution) || (ExecutionMode(executionMode))
                    "evolution"
                elseif (executionMode == doAllSequencesTest) || (ExecutionMode(executionMode))
                    "test_all_sequences"
                elseif executionMode == doLongSequencesTest || (ExecutionMode(executionMode))
                    "long_sequence"
                elseif executionMode == doFixedSequenceTest || (ExecutionMode(executionMode))
                    "fixed_sequence_activity"
                end
            else
                outputdir
            end

        new(
            executionMode, getAlmostCorrectNWs, writeNetworkActivity,
            resolved_outputdir, minConnectionWeight, maxConnectionWeight, timeStep, neuronalType,
            noOfInputs, noOfinterNeurons, noOfSignals, noOfOutputs, ge_gain, gi_gain,
            synapticDelay, minWeightThreshold, maxWeightRhreshold, popSize, eliteCount,
            randomizeCount, randomizeEveryXGen, weightDeletionProb, mutationProb,
            signChangeProb, mutationMean, mutationStength, maxGen, maxRuns,
            weightReductionProb, weightReductionStrength, noiseVectorSize,
            variationOnSignal, variationOnSilence, gaussianNoiseOnVoltage, gMean, gStdDev,
            letterSize, silenctInterval, noOfLetters, reevaluateSeq,
            noOfLetters * (letterSize + silenctInterval)
        )
    end
end

get_nodes_in_network(params::Parameters) =
    params.noOfInputs + params.noOfinterNeurons + params.noOfOutputs

end # module
