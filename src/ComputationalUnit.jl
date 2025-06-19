module CompUnit

export ComputationalUnit,
    resetVoltage!,
    hasSpike,
    resetAdaptation,
    updateAdaptation,
    updateVoltage,
    updateExcitatoryCond,
    updateInhibitoryCond

# import DrWatson: srcdir, @quickactivate
# @quickactivate "Evosnn"
# "enum-type.jl" |> srcdir |> include
using ..MyTypes: NeuronType, LIF, ADX, IZH

mutable struct ComputationalUnit
    inv_C::Float64
    inv_taum::Float64
    DeltaT::Float64
    aDiv1000::Float64
    b::Float64
    i_offset::Float64
    inv_tauw::Float64

    Vreset::Float64
    Vresting::Float64
    Vthreshold::Float64
    Vspike::Float64

    inv_tau_syn_E::Float64
    inv_tau_syn_I::Float64

    e_rev_E::Float64
    e_rev_I::Float64

    absRefractoryTime::Float64

    spikeCount::Int
    hasSpiked::Bool

    nModel::NeuronType

    function ComputationalUnit(neuronalType::NeuronType)
        if neuronalType == LIF
            inv_C = 1.0 / 0.2
            inv_taum = 1.0 / 20.0

            DeltaT = -1 # unused
            aDiv1000 = 1.0 # unused
            b = 0.0 # unused
            i_offset = 0.0
            inv_tauw = 1.0 / 100.0 # unused

            Vreset = -70.0
            Vresting = -70.0
            local Vthreshold = -50.0
            Vspike = 0.0

            inv_tau_syn_E = 1.0 / 5.0
            inv_tau_syn_I = 1.0 / 5.0

            e_rev_E = 0.0
            e_rev_I = -70.0

            absRefractoryTime = 0.0

            spikeCount = 0
            hasSpiked = false

            nModel = neuronalType

            new(inv_C,
                inv_taum,
                DeltaT,
                aDiv1000,
                b,
                i_offset,
                inv_tauw,
                Vreset,
                Vresting,
                Vthreshold,
                Vspike,
                inv_tau_syn_E,
                inv_tau_syn_I,
                e_rev_E,
                e_rev_I,
                absRefractoryTime,
                spikeCount,
                hasSpiked,
                nModel
            )
        elseif neuronalType == ADX

            inv_C = 1.0 / 0.20
            inv_taum = 1.0 / 20.0
            DeltaT = 2.0

            aDiv1000 = 2.0 / 1000.0
            b = 0.0
            i_offset = 0.0
            inv_tauw = 1.0 / 30.0

            Vreset = -58.0
            Vresting = -70.0
            Vthreshold = -50.0
            Vspike = 0.0

            inv_tau_syn_E = 1.0 / 5.0
            inv_tau_syn_I = 1.0 / 5.0

            e_rev_E = 0.0
            e_rev_I = -70.0

            absRefractoryTime = 0.0

            spikeCount = 0
            hasSpiked = false

            nModel = neuronalType

            new(inv_C,
                inv_taum,
                DeltaT,
                aDiv1000,
                b,
                i_offset,
                inv_tauw,
                Vreset,
                Vresting,
                Vthreshold,
                Vspike,
                inv_tau_syn_E,
                inv_tau_syn_I,
                e_rev_E,
                e_rev_I,
                absRefractoryTime,
                spikeCount,
                hasSpiked,
                nModel
            )
        else

            throw(ArgumentError("Unknown neuronal type"))
        end
    end
end

# # This is moved to constructor
# function createNeuron!(cu::ComputationalUnit, params::Parameters)
# end

function resetVoltage!(cu::ComputationalUnit)
    if cu.nModel == ADX
        # cu.Vreset = -58.0
        return -58.0
    elseif cu.nModel == LIF
        # cu.Vreset = -70.0
        return -70.0
    else
        throw(ArgumentError("Unknown neuronal type"))
    end
end

function hasSpike(cu::ComputationalUnit, voltage::Float64)
    VCut =
        if (cu.nModel == LIF)
            cu.Vthreshold
        elseif (cu.nModel == ADX)
            cu.Vspike
        else
            cu.Vthreshold
        end

    return voltage >= VCut
end

function resetAdaptation(cu::ComputationalUnit, Adaptation::Float64)
    return Adaptation + cu.b
end

"""
Parameter `time_step` was a global variable used, so on every usage, it will have to be provided as the last arguemnt
"""
function updateAdaptation(cu::ComputationalUnit, adaptation::Float64, currentVoltage::Float64, timeStep::Float64)
    if cu.nModel == ADX
        return adaptation + ((cu.aDiv1000 * (currentVoltage - cu.Vresting) - adaptation) * cu.inv_tauw) * timeStep
    else
        return adaptation  # default behavior for models not defined
    end
end

"""
Parameter `time_step` was a global variable used, so on every usage, it will have to be provided as the last arguemnt
"""
function updateVoltage(cu::ComputationalUnit, currentVoltage::Float64, excitatory_conductance::Float64, inhibitory_conductance::Float64, Adaptation::Float64, timeStep::Float64)
    voltageChange = if cu.nModel == LIF
        (
            (cu.Vresting - currentVoltage) * cu.inv_taum + (excitatory_conductance * (cu.e_rev_E - currentVoltage) +
                                                            inhibitory_conductance * (cu.e_rev_I - currentVoltage) +
                                                            cu.i_offset) *
                                                           cu.inv_C
        ) * timeStep
    elseif cu.nModel == ADX
        (
            (cu.Vresting - currentVoltage + cu.DeltaT * exp((currentVoltage - cu.Vthreshold) / cu.DeltaT)
            ) * cu.inv_taum + (
                excitatory_conductance * (cu.e_rev_E - currentVoltage) + inhibitory_conductance *
                                                                         (cu.e_rev_I - currentVoltage) + cu.i_offset - Adaptation) * cu.inv_C
        ) * timeStep
    elseif cu.nModel == IZH
        (0.04 * currentVoltage^2 + 5 * currentVoltage + 140 - Adaptation +
         (excitatory_conductance * (cu.e_rev_E - currentVoltage) +
          inhibitory_conductance * (cu.e_rev_I - currentVoltage) + cu.i_offset) * cu.inv_C
        ) * timeStep
    end
    return currentVoltage + voltageChange
end

"""
Parameter `time_step` was a global variable used, so on every usage, it will have to be provided as the last arguemnt
"""
function updateExcitatoryCond(cu::ComputationalUnit, current_excitatory_conductance::Float64, timeStep::Float64)
    return current_excitatory_conductance - current_excitatory_conductance * cu.inv_tau_syn_E * timeStep
end

"""
Parameter `time_step` was a global variable used, so on every usage, it will have to be provided as the last arguemnt
"""
function updateInhibitoryCond(cu::ComputationalUnit, current_inhibitory_conductance::Float64, timeStep::Float64)
    return current_inhibitory_conductance - current_inhibitory_conductance * cu.inv_tau_syn_I * timeStep
end


end # CompUnit
