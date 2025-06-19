module ParamsModule

using DrWatson: datadir

export Parameters

using ..MyTypes: NeuronType, LIF, ADX, IZH
using ..MyTypes: ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest

struct Parameters
    doEvolution::Bool
    getAlmostCorrectNWs::Bool
    writeNetworkActivity::Bool
    doAllSequencesTest::Bool

    outputdir::String

    minConnectionWeight::Float64
    maxConnectionWeight::Float64
    timeStep::Float64
    neuronal_type::NeuronType

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
    variationOnSignal::Float64
    variationOnSilence::Float64

    gaussianNoiseOnVoltage::Bool
    gMean::Float64
    gStdDev::Float64

    letterSize::Int
    silenctInterval::Int
    noOfLetters::Int
    reevaluateSeq::Int
    sequenceSize::Int

    function Parameters(;
        doEvolution::Bool=false,
        getAlmostCorrectNWs::Bool=false,
        writeNetworkActivity::Bool=false,
        doAllSequencesTest::Bool=false,
        outputdir::String=datadir("outputs5/"),
        minConnectionWeight::Float64=-10.0,
        maxConnectionWeight::Float64=10.0,
        timeStep::Float64=1.0,
        neuronalType::NeuronType=ADX,
        noOfInputs::Int=10,
        noOfinterNeurons::Int=10,
        noOfSignals::Int=10,
        noOfOutputs::Int=1,
        ge_gain::Float64=0.007,
        gi_gain::Float64=0.007, synapticDelay::Float64=4.0,
        minWeightThreshold::Float64=0.0,
        maxWeightRhreshold::Float64=0.0, popSize::Int=100,
        eliteCount::Int=5,
        randomizeCount::Int=5,
        randomizeEveryXGen::Int=10,
        weightDeletionProb::Float64=0.0,
        mutationProb::Float64=0.3,
        signChangeProb::Float64=0.0,
        mutationMean::Float64=0.5,
        mutationStength::Float64=1.8,
        maxGen::Int=10000,
        maxRuns::Int=10, weightReductionProb::Float64=0.2,
        weightReductionStrength::Float64=0.1, noiseVectorSize::Int=1000,
        variationOnSignal::Float64=0.0,
        variationOnSilence::Float64=0.0, gaussianNoiseOnVoltage::Bool=false,
        gMean::Float64=0.0,
        gStdDev::Float64=0.5, letterSize::Int=6,
        silenctInterval::Int=24,
        noOfLetters::Int=1500,
        reevaluateSeq::Int=10000
    )
        new(
            executionMode, getAlmostCorrectNWs, writeNetworkActivity,
            outputdir, minConnectionWeight, maxConnectionWeight, timeStep, neuronalType,
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

end # module
