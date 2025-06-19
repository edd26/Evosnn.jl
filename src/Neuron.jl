module NeuronModule

export Neuron

using ..MyTypes: NeuronType, LIF, ADX, IZH
using ..CompUnit

mutable struct Neuron
    voltage::Float64
    adaptation::Float64
    exConductance::Float64
    inConductance::Float64
    synapticDelayBuffer::Int
    delayCount::Int
    refrectoryPeriod::Int
    cUnit::ComputationalUnit
    spikeBitmap::Vector{Bool}  # Assuming spikeBitmap is an integer vector
    voltageBuffer::Vector{Float64}

    function Neuron(neuronalType::NeuronType; voltage::Float64=-70.0)
        # @warn "Those two neuronal types are the same upon construction"
        # voltage = -70.0
        adaptation = 0.0
        exConductance = 0.0
        inConductance = 0.0
        synapticDelayBuffer = 0
        delayCount = 0
        refrectoryPeriod = 0
        spikeBitmap = Bool[]
        voltageBuffer = Float64[]
        cUnit = ComputationalUnit(neuronalType)

        return new(
            voltage,
            adaptation,
            exConductance,
            inConductance,
            synapticDelayBuffer,
            delayCount,
            refrectoryPeriod,
            cUnit,
            spikeBitmap,
            voltageBuffer,
        )
        # elseif neuronalType == ADX
        #     return new(
        #         voltage,
        #         adaptation,
        #         exConductance,
        #         inConductance,
        #         synapticDelayBuffer,
        #         delayCount,
        #         refrectoryPeriod,
        #         spikeBitmap,
        #         voltageBuffer,
        #     )
        # else
        #     throw(ArgumentError("Unknown neuronal type"))
        # end
    end
end
end # module
