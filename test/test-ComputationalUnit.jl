using DrWatson, Test
@quickactivate "Evosnn"

@testset "Computational Unit Tests" begin
    # "ComputationalUnit.jl" |> srcdir |> include

    # Define a mock NeuronType for testing purposes

    # const 
    # LIF = NeuronType(LIF)
    # ADX = NeuronType(ADX)

    # Tests for ComputationalUnit

    function test_ComputationalUnit_LIF()
        cu = ComputationalUnit(LIF)

        @test cu.Vreset == -70.0
        @test cu.Vthreshold == -50.0
        @test cu.Vspike == 0.0
        @test cu.absRefractoryTime == 0.0
        @test cu.nModel == LIF
    end

    function test_ComputationalUnit_ADX()
        cu = ComputationalUnit(ADX)

        @test cu.Vreset == -58.0
        @test cu.Vthreshold == -50.0
        @test cu.Vspike == 0.0
        @test cu.absRefractoryTime == 0.0
        @test cu.DeltaT == 2.0
        @test cu.nModel == ADX
    end

    function test_resetVoltage!()
        cu_LIF = ComputationalUnit(LIF)
        cu_ADX = ComputationalUnit(ADX)

        resetVoltage!(cu_LIF)
        resetVoltage!(cu_ADX)

        @test cu_LIF.Vreset == -70.0
        @test cu_ADX.Vreset == -58.0
    end

    function test_hasSpike()
        cu_LIF = ComputationalUnit(LIF)
        cu_ADX = ComputationalUnit(ADX)

        @test hasSpike(cu_LIF, -49.0) == true
        @test hasSpike(cu_LIF, -51.0) == false
        @test hasSpike(cu_ADX, 0.0) == true
        @test hasSpike(cu_ADX, -70.0) == false
    end

    function test_updateAdaptation()
        cu_ADX = ComputationalUnit(ADX)
        Adaptation = 0.0
        current_voltage = 20.0
        time_step = 1.0 # ms

        updated_value = updateAdaptation(cu_ADX, Adaptation, current_voltage, time_step)
        @test updated_value == 0.006
    end

    function test_updateVoltage()
        cu_LIF = ComputationalUnit(LIF)
        cu_ADX = ComputationalUnit(ADX)

        currentVoltage = -65.0
        excitatory_conductance = 10.0
        inhibitory_conductance = 5.0
        Adaptation = 1.0
        time_step = 1.0 #ms

        new_voltage_LIF = updateVoltage(cu_LIF, currentVoltage, excitatory_conductance, inhibitory_conductance, Adaptation, time_step)
        new_voltage_ADX = updateVoltage(cu_ADX, currentVoltage, excitatory_conductance, inhibitory_conductance, Adaptation, time_step)

        @test new_voltage_LIF != currentVoltage
        @test new_voltage_ADX != currentVoltage
    end

    function test_updateExcitatoryCond()
        cu = ComputationalUnit(LIF)
        current_excitatory_conductance = 10.0
        time_step = 1.0 # ms

        new_cond = updateExcitatoryCond(cu, current_excitatory_conductance, time_step)

        @test new_cond < current_excitatory_conductance
    end

    function test_updateInhibitoryCond()
        cu = ComputationalUnit(LIF)
        current_inhibitory_conductance = 10.0
        time_step = 1.0 # ms

        new_cond = updateInhibitoryCond(cu, current_inhibitory_conductance, time_step)

        @test new_cond < current_inhibitory_conductance
    end

    # Run all tests
    test_ComputationalUnit_LIF()
    test_ComputationalUnit_ADX()
    test_resetVoltage!()
    test_hasSpike()
    test_updateAdaptation()
    test_updateVoltage()
    test_updateExcitatoryCond()
    test_updateInhibitoryCond()

end
