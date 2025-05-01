using DrWatson, Test
@quickactivate "Evosnn"


@testset "Basic Parameters initialisation" begin
    p = Parameters()

    @test p.doEvolution == false
    @test p.getAlmostCorrectNWs == false
    @test p.writeNetworkActivity == false
    @test p.doAllSequencesTest == false

    @test p.outputdir == datadir("outputs5/")
    @test p.minConnectionWeight == -10.0
    @test p.maxConnectionWeight == 10.0
    @test p.timeStep == 1.0
    @test p.neuronal_type == ADX

    @test p.noOfInputs == 10
    @test p.noOfinterNeurons == 10
    @test p.noOfSignals == 10
    @test p.noOfOutputs == 1
    @test p.ge_gain == 0.007
    @test p.gi_gain == 0.007

    @test p.synapticDelay == 4.0
    @test p.minWeightThreshold == 0.0
    @test p.maxWeightRhreshold == 0.0

    @test p.popSize == 100
    @test p.eliteCount == 5
    @test p.randomizeCount == 5
    @test p.randomizeEveryXGen == 10
    @test p.weightDeletionProb == 0.0
    @test p.mutationProb == 0.3
    @test p.signChangeProb == 0.0
    @test p.mutationMean == 0.5
    @test p.mutationStength == 1.8
    @test p.maxGen == 10000
    @test p.maxRuns == 10

    @test p.weightReductionProb == 0.2
    @test p.weightReductionStrength == 0.1

    @test p.noiseVectorSize == 1000
    @test p.variationOnSignal == 0.0
    @test p.variationOnSilence == 0.0

    @test p.gaussianNoiseOnVoltage == false
    @test p.gMean == 0.0
    @test p.gStdDev == 0.5

    @test p.letterSize == 6
    @test p.silenctInterval == 24
    @test p.noOfLetters == 1500
    @test p.reevaluateSeq == 10000
    @test p.sequenceSize == p.noOfLetters * (p.letterSize + p.silenctInterval)
end
