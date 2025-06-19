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
    makeIndividualWithFixedInputOutputConnections

export printIndividualMatrix,
    sumOfConnectionWeights,
    removeLowWeights

# export get_abcdefhXXX_XXXdefghij_Sequence

export Ga,
    run_Ga!

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

include("Individual.jl")
using .IndividualModule

include("utils.jl")


include("UtilityFunctions.jl")
using .UtilityFunctions

include("Ga.jl")
using .GaModule




# A note to self- if there is a list of dependencies, they have to be imported into main package file, then they have to be references in the subpackage via double dot

end # module Evosnn
