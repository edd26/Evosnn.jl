using DrWatson, Test
@quickactivate "Evosnn"

@testset "Neuron Tests" begin
    # Define a mock Parameters structure for testing

    # "Individual.jl" |> srcdir |> include

    # Test for LIF type
    neuron_LIF = Neuron(LIF)

    @test neuron_LIF.voltage == -70.0
    @test neuron_LIF.adaptation == 0.0
    @test neuron_LIF.exConductance == 0.0
    @test neuron_LIF.inConductance == 0.0
    @test neuron_LIF.synapticDelayBuffer == 0
    @test neuron_LIF.delayCount == 0
    @test isempty(neuron_LIF.spikeBitmap)
    @test isempty(neuron_LIF.voltageBuffer)
    @test neuron_LIF.refrectoryPeriod == 0

    # Test for ADX type
    neuron_ADX = Neuron(ADX)

    @test neuron_ADX.voltage == -70.0
    @test neuron_ADX.adaptation == 0.0
    @test neuron_ADX.exConductance == 0.0
    @test neuron_ADX.inConductance == 0.0
    @test neuron_ADX.synapticDelayBuffer == 0
    @test neuron_ADX.delayCount == 0
    @test isempty(neuron_ADX.spikeBitmap)
    @test isempty(neuron_ADX.voltageBuffer)
    @test neuron_ADX.refrectoryPeriod == 0
end
