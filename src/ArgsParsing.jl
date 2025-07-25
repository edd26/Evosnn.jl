
module EvosnnArgs
using ArgParse
"enum-type.jl" |> include

using .MyTypes: NeuronType, ExecutionMode
# using MyTypes

function parse_parameters_args()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--executionMode"
        help = "Execution mode; possible options are 0 for doEvolution, 1 for doAllSequencesTest, 2 for doLongSequencesTest, 3 for doFixedSequenceTest"
        arg_type = Int
        default = 0

        "--selected_size"
        help = "If set, then sets noInputs, noOfinterNeurons, noOfSignals, disregarding to what those args are set"
        arg_type = Union{Int,Nothing}
        default = nothing

        "--getAlmostCorrectNWs"
        arg_type = Bool
        default = false

        "--writeNetworkActivity"
        arg_type = Bool
        default = false

        "--outputdir"
        arg_type = String
        default = "outputs5/"

        "--minConnectionWeight"
        arg_type = Float64
        default = -10.0

        "--maxConnectionWeight"
        arg_type = Float64
        default = 10.0

        "--timeStep"
        arg_type = Float64
        default = 1.0

        "--neuronalType"
        help = "Type of neuron; use 0 for LIF, 1 for ADX, 2 for IZH"
        arg_type = Int
        default = 1


        "--noOfInputs"
        arg_type = Int
        default = 10

        "--noOfinterNeurons"
        arg_type = Int
        default = 10

        "--noOfSignals"
        arg_type = Int
        default = 10

        "--noOfOutputs"
        arg_type = Int
        default = 1

        "--ge_gain"
        arg_type = Float64
        default = 0.007

        "--gi_gain"
        arg_type = Float64
        default = 0.007

        "--synapticDelay"
        arg_type = Float64
        default = 4.0

        "--minWeightThreshold"
        arg_type = Float64
        default = 0.0

        "--maxWeightRhreshold"
        arg_type = Float64
        default = 0.0

        "--popSize"
        arg_type = Int
        default = 100

        "--eliteCount"
        arg_type = Int
        default = 5

        "--randomizeCount"
        arg_type = Int
        default = 5

        "--randomizeEveryXGen"
        arg_type = Int
        default = 10

        "--weightDeletionProb"
        arg_type = Float64
        default = 0.0

        "--mutationProb"
        arg_type = Float64
        default = 0.1

        "--signChangeProb"
        arg_type = Float64
        default = 0.0

        "--mutationMean"
        arg_type = Float64
        default = 0.5

        "--mutationStength"
        arg_type = Float64
        default = 0.7

        "--maxGen"
        arg_type = Int
        default = 10000

        "--maxRuns"
        arg_type = Int
        default = 10

        "--weightReductionProb"
        arg_type = Float64
        default = 0.2

        "--weightReductionStrength"
        arg_type = Float64
        default = 0.1

        "--noiseVectorSize"
        arg_type = Int
        default = 1000

        "--variationOnSignal"
        arg_type = Int
        default = 0

        "--variationOnSilence"
        arg_type = Int
        default = 0

        "--gaussianNoiseOnVoltage"
        arg_type = Bool
        default = true

        "--gMean"
        arg_type = Float64
        default = 0.0

        "--gStdDev"
        arg_type = Float64
        default = 0.5

        "--letterSize"
        arg_type = Int
        default = 6

        "--silenctInterval"
        arg_type = Int
        default = 24

        "--noOfLetters"
        arg_type = Int
        default = 100000 # 6000

        "--reevaluateSeq"
        arg_type = Int
        default = 10000
    end

    return parse_args(s; as_symbols=true)
end

function refactorArgs(args)
    if !isnothing(args[:selected_size])
        args[:noOfInputs] = args[:noOfinterNeurons] = args[:noOfSignals] = args[:selected_size]
    end
    pop!(args, :selected_size)

    if !(args[:executionMode] isa ExecutionMode)
        args[:executionMode] = ExecutionMode(args[:executionMode])
    end
    if !(args[:neuronalType] isa NeuronType)
        args[:neuronalType] = NeuronType(args[:neuronalType])
    end

    return args
end

end # module