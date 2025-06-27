module Evosnn

export NeuronType, LIF, ADX, IZH
export ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest, doFixedSequenceTest

export ComputationalUnit,
    resetVoltage!,
    hasSpike,
    resetAdaptation,
    updateAdaptation,
    updateVoltage,
    updateExcitatoryCond,
    updateInhibitoryCond

export Neuron

export Parameters,
    get_nodes_in_network

export Individual,
    makeIndividualWithFixedInputOutputConnections,
    networkStep!,
    activateOutput,
    outputNetworkActivity,
    readIndividualMatrix,
    copy_individual,
    randomizeinterconnections

export printIndividualMatrix,
    sumOfConnectionWeights,
    removeLowWeights

export get_abcdefgXXX_XXXdefghij_Sequence,
    getCorrectPatternsMarkers,
    # getCorrectPatternsMarkersABCDEFGHIJ,
    # getCorrectPatternsMarkersABCDEFGHI,
    # getCorrectPatternsMarkersABCDEFGH,
    # getCorrectPatternsMarkersABCDEFG,
    # getCorrectPatternsMarkersABCDEF,
    # getCorrectPatternsMarkersABCDE,
    # getCorrectPatternsMarkersABCD,
    # getCorrectPatternsMarkersABC,
    # getCorrectPatternsMarkersAB,
    get_abcdefgXXX_XXXdefghij_Sequence,
    get_abcdefXXX_XXXdefghi_Sequence,
    get_abcdeXXX_XXXfgh_Sequence,
    get_abcdXXX_XXXdefg_Sequence,
    get_abcXXX_XXXdef_Sequence,
    getABXXX_XXXDE_Sequence,
    getABCDSequence,
    checkOrCreateDirectory,
    get_signal_of_len

export PatternFrequencyPair,
    empty!,
    Ga,
    run_Ga!,
    run_Ga_parallel!,
    fitness,
    reEvaluateAllPerm5_5sig,
    reEvaluateAllPerm6_5sig,
    reEvaluateAllPerm7_5sig,
    reEvaluateAllPerm6_6sig,
    reEvaluateAllPerm7_6sig,
    reEvaluateAllPerm8_6sig,
    reEvaluateAllPerm7_7sig,
    reEvaluateAllPerm8_7sig,
    reEvaluateAllPerm9_7sig,
    reEvaluateAllPerm,
    reEvaluateUserDefinedSequence,
    reEvaluateOnLargeSequence

# ===-===-
include("enum-type.jl")
using .MyTypes: NeuronType, LIF, ADX, IZH
using .MyTypes: ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest, doFixedSequenceTest

include("ComputationalUnit.jl")
using .CompUnit: ComputationalUnit,
    resetVoltage!,
    hasSpike,
    resetAdaptation,
    updateAdaptation,
    updateVoltage,
    updateExcitatoryCond,
    updateInhibitoryCond

include("Neuron.jl")
using .NeuronModule: Neuron

include("Parameters.jl")
using .ParamsModule: Parameters, get_nodes_in_network

include("UtilityFunctions.jl")
using .UtilityFunctions

include("Individual.jl")
using .IndividualModule

include("Ga.jl")
using .GaModule

include("utils.jl")



# A note to self- if there is a list of dependencies, they have to be imported into main package file, then they have to be references in the subpackage via double dot

end # module Evosnn
