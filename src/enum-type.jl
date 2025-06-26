module MyTypes

export NeuronType, LIF, ADX, IZH,
    ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest, doFixedSequenceTest

@enum NeuronType begin
    LIF
    ADX
    IZH
end # enum NeuronType


@enum ExecutionMode begin
    doEvolution
    doAllSequencesTest
    doLongSequencesTest
    doFixedSequenceTest
end # enum ExecutionMode
end # module