module Evosnn

export Parameters

export NeuronType, LIF, ADX, IZH

export ComputationalUnit, 
       resetVoltage!,
        hasSpike,
        resetAdaptation,
        updateAdaptation,
        updateVoltage,
        updateExcitatoryCond,
        updateInhibitoryCond

export Neuron
export Individual


export printIndividualMatrix,
    sumOfConnectionWeights,
    removeLowWeights


include("Parameters.jl")

include("Neuron.jl")
using .NeuronModule: Neuron

# include("enum-type.jl")
using .MyTypes: NeuronType, LIF, ADX, IZH

include("ComputationalUnit.jl")
using .CompUnit: ComputationalUnit, 
       resetVoltage!,
        hasSpike,
        resetAdaptation,
        updateAdaptation,
        updateVoltage,
        updateExcitatoryCond,
        updateInhibitoryCond


include("Individual.jl")
include("utils.jl")


end # module Evosnn
