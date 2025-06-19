module Evosnn

export NeuronType, LIF, ADX, IZH
export ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest

export ComputationalUnit,
    resetVoltage!,
    hasSpike,
    resetAdaptation,
    updateAdaptation,
    updateVoltage,
    updateExcitatoryCond,
    updateInhibitoryCond

export Neuron

export Parameters

export Individual,
    makeIndividualWithFixedInputOutputConnections,
    networkStep!,
    activateOutput

export printIndividualMatrix,
    sumOfConnectionWeights,
    removeLowWeights

export get_abcdefhXXX_XXXdefghij_Sequence,
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
    get_abcdefhXXX_XXXdefghij_Sequence,
    get_abcdefXXX_XXXdefghi_Sequence,
    get_abcdeXXX_XXXfgh_Sequence,
    get_abcdXXX_XXXdefg_Sequence,
    getABXXX_XXXDE_Sequence,
    getABCDSequence,
    checkOrCreateDirectory

export PatternFrequencyPair,
    empty!,
    Ga,
    run_Ga!,
    fitness,
    reEvaluateAllPerm5_5sig,
    reEvaluateAllPerm6_5sig,
    reEvaluateAllPerm7_5sig,
    reEvaluateAllPerm6_6sig,
    reEvaluateAllPerm7_6sig,
    reEvaluateAllPerm8_6sig,
    reEvaluateAllPerm7_7sig,
    reEvaluateAllPerm8_7sig,
    reEvaluateAllPerm9_7sig

# ===-===-
include("enum-type.jl")
using .MyTypes: NeuronType, LIF, ADX, IZH
using .MyTypes: ExecutionMode, doEvolution, doAllSequencesTest, doLongSequencesTest

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
using .ParamsModule: Parameters

include("UtilityFunctions.jl")
using .UtilityFunctions

include("Individual.jl")
using .IndividualModule

include("Ga.jl")
using .GaModule

include("utils.jl")



# A note to self- if there is a list of dependencies, they have to be imported into main package file, then they have to be references in the subpackage via double dot

end # module Evosnn
