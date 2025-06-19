module GaModule
using Random

export PatternFrequencyPair, empty!
export Ga, run_Ga!, fitness,
    reEvaluateAllPerm5_5sig,
    reEvaluateAllPerm6_5sig,
    reEvaluateAllPerm7_5sig,
    reEvaluateAllPerm6_6sig,
    reEvaluateAllPerm7_6sig,
    reEvaluateAllPerm8_6sig,
    reEvaluateAllPerm7_7sig,
    reEvaluateAllPerm8_7sig,
    reEvaluateAllPerm9_7sig

using ..IndividualModule: Individual, setInput!, networkStep!, replicate!
using ..UtilityFunctions
using ..ParamsModule
using ..UtilityFunctions: insertGapsAndSetLetterSize

using Base.Threads

mutable struct PatternFrequencyPair
    pattern::String
    freq::Int
    PatternFrequencyPair(; pattern::String="", freq::Int=0) = new(pattern, freq)
end

function clear!(patternFrequencyPair::PatternFrequencyPair)
    patternFrequencyPair.pattern = ""
    patternFrequencyPair.pattern = 0
end

function empty!(pattFrqStructList::Vector{PatternFrequencyPair})
    pattFrqStructList = Vector{PatternFrequencyPair}()
end

struct Ga
    history
    keepGen0Elites
    signalSiquence::Vector{Char}
    correctIndecies::Vector{Int}

    function Ga()
        history = []
        keepGen0Elites = []
        signalSiquence = []
        correctIndecies = []
        return new(history, keepGen0Elites, signalSiquence, correctIndecies)
    end
end


function run_Ga!(ga::Ga, pop::Vector{Individual}, genNo::Int, hardPatternSeq::Vector{Char}, params::Parameters)
    alphabet::String = "ABCDEFGHIJKLMNOPQRSTUWXYZ"
    signal::String = alphabet[1:params.noOfSignals]
    noOfSeq::Int = 6
    step::Int = 0

    signalSiquence = Vector{Char}()
    correctIndecies = Vector{Int}()
    all_insertions = Int[]

    if params.noOfSignals == 10
        signalSiquence, all_insertions = get_abcdefhXXX_XXXdefghij_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDEFGHIJ(signalSiquence, signal)
    elseif params.noOfSignals == 9
        signalSiquence, all_insertions = get_abcdefXXX_XXXdefghi_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDEFGHI(signalSiquence, signal)
    elseif params.noOfSignals == 8
        signalSiquence, all_insertions = get_abcdeXXX_XXXfgh_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDEFGH(signalSiquence, signal)
    elseif params.noOfSignals == 7
        signalSiquence, all_insertions = get_abcdXXX_XXXdefg_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
    elseif params.noOfSignals == 6
        signalSiquence, all_insertions = get_abcXXX_XXXdef_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
    elseif params.noOfSignals == 5
        signalSiquence, all_insertions = getABXXX_XXXDE_Sequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
    elseif params.noOfSignals == 4
        signalSiquence, all_insertions = getABCDSequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
    elseif params.noOfSignals == 3
        signalSiquence, all_insertions = getABCSequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        # correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
        # if params.getAlmostCorrectNWs
        #     correctIndecies = getCorrectPatternsMarkersAB(signalSiquence, "BC")
        # end
    elseif params.noOfSignals == 2
        signalSiquence = getABSequence(signal, params.noOfLetters, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
        # correctIndecies = getCorrectPatternsMarkersAB(signalSiquence, signal)
    end
    correctIndecies = getCorrectPatternsMarkers(signalSiquence, signal)
    # expanded_all_insertions = [i * 30 - 29 for i in all_insertions]
    # @assert issubset(expanded_all_insertions, correctIndecies)

    # Above this line seems to be working ok for size 3 and 4
    for i in 1:length(pop)
        total_elements_in_signal = length(signalSiquence)
        # @info "Total elements in the sequence $(total_elements_in_signal )"
        for j in 1:total_elements_in_signal# ÷100
            # j % 100 == 0 && @info "At iteration $(j)"
            step += 1
            setInput!(pop[i], signalSiquence[j], j, params.writeNetworkActivity)
            networkStep!(pop[i], step, params)
        end

        pop[i].fitness = fitness(pop[i], signalSiquence, correctIndecies, params)

        step = 0
        pop[i].outputNeurons[1].spikeBitmap = Bool[]
        pop[i].outputNeurons[1].voltageBuffer = Float64[]
        for p in 1:length(pop[i].inputNeurons)
            pop[i].inputNeurons[p].voltageBuffer = Float64[]
        end
        for q in 1:length(pop[i].interNeurons)
            pop[i].interNeurons[q].voltageBuffer = Float64[]
        end
    end

    sort!(pop, by=x -> x.fitness)

    for k in params.eliteCount+1:length(pop)
        firstRand = rand(1:length(pop))
        secondRand = rand(1:length(pop))
        if pop[firstRand].fitness < pop[secondRand].fitness
            pop[k] = pop[firstRand]
        else
            pop[k] = pop[secondRand]
        end
    end

    for rep in params.eliteCount+1:length(pop)
        replicate!(pop[rep], params)
    end

    signalSiquence = []
    correctIndecies = []
end

function runMinimalWeight(ga::Ga, pop::Vector{Individual}, genNo::Int, params::Parameters)
    signal::String = "ABC"
    letterSize::Int = params.letterSize
    gapSize::Int = params.silenctInterval
    sequenceLength::Int = params.noOfLetters
    noOfSeq::Int = 6
    step::Int = 0
    reward::Float64 = 0.0
    penalty::Float64 = 0.0
    rewardAll::Float64 = 0.0
    penaltyAll::Float64 = 0.0
    signalSiquence::Vector{Char} = getRandomSequence(signal, sequenceLength, gapSize, letterSize)
    # correctIndecies::Vector{Int} = getCorrectPatternsMarkersABC(signalSiquence, signal)
    correctIndecies::Vector{Int} = getCorrectPatternsMarkers(signalSiquence, signal)

    for i in 1:length(pop)
        reward = 0.0
        penalty = 0.0
        rewardAll = 0.0
        penaltyAll = 0.0

        for j in 1:length(signalSiquence)
            step += 1
            pop[i].setInput!(signalSiquence[j], j)
            pop[i].networkStep(step)
        end

        rCountint = 0.0
        rCountSig = 0.0
        pCountint = 0.0
        pCountSig = 0.0

        k = 0
        sizeHistory = length(pop[i].outputNeurons[1].spikeBitmap)

        for j in correctIndecies
            while k < sizeHistory && k != j
                while signalSiquence[k] != 'Z' && k != j
                    if pop[i].outputNeurons[1].spikeBitmap[k]
                        pCountSig += 1.0
                        k += 1
                    else
                        k += 1
                    end
                end
                while signalSiquence[k] == 'Z' && k != j
                    if pop[i].outputNeurons[1].spikeBitmap[k]
                        pCountint += 1.0
                        k += 1
                    else
                        k += 1
                    end
                end
                if pCountint + pCountSig > 0.0
                    pCountint = 0.0
                    pCountSig = 0.0
                    penalty += 1.0
                end
            end
            while k < sizeHistory && signalSiquence[k] != 'Z'
                if pop[i].outputNeurons[1].spikeBitmap[k]
                    rCountSig += 1.0
                    k += 1
                else
                    k += 1
                end
            end
            while k < sizeHistory && signalSiquence[k] == 'Z'
                if pop[i].outputNeurons[1].spikeBitmap[k]
                    rCountint += 1.0
                    k += 1
                else
                    k += 1
                end
            end
            if rCountSig + rCountint > 0.0
                rCountSig = 0.0
                rCountint = 0.0
                reward += 1.0
            end
        end

        pop[i].rewardn = reward / length(correctIndecies)
        pop[i].penaltyn = penalty / (sequenceLength - length(correctIndecies))
        pop[i].reward = reward
        pop[i].penalty = penalty
        pop[i].fitness = 1 - (pop[i].rewardn - 4 * pop[i].penaltyn)

        step = 0
        pop[i].outputNeurons[1].spikeBitmap = []
        pop[i].outputNeurons[1].voltageBuffer = []
        for p in 1:length(pop[i].inputNeurons)
            pop[i].inputNeurons[p].voltageBuffer = []
        end
        for q in 1:length(pop[i].interNeurons)
            pop[i].interNeurons[q].voltageBuffer = []
        end
    end

    sort!(pop, by=x -> x.fitness)

    wi = 1
    while pop[wi].fitness < 0.005
        pop[wi].sumOfConnectionWeights()
        wi += 1
    end
    pop = pop[1:wi]

    if wi != 0
        for k in wi+1:params.popSize
            push!(pop, pop[rand(1:wi)])
        end

        sort!(pop, by=x -> x.absWeightSum)

        if genNo == 0
            for ke in 1:params.eliteCount
                push!(ga.keepGen0Elites, pop[ke])
            end
        else
            for ei in 1:length(ga.keepGen0Elites)
                unshift!(pop, ga.keepGen0Elites[ei])
            end
        end

        sort!(pop, by=x -> x.absWeightSum)
        pop = pop[1:params.popSize]

        for k in params.eliteCount+1:length(pop)
            firstRand = rand(1:length(pop))
            secondRand = rand(1:length(pop))
            if pop[firstRand].absWeightSum < pop[secondRand].absWeightSum
                pop[k] = pop[firstRand]
            else
                pop[k] = pop[secondRand]
            end
        end

        for rep in params.eliteCount+1:length(pop)
            pop[rep].replicateByOnlyReducingWeights()
        end
    else
        println("wi == 0")
        return false
    end

    return true
end

function reEvaluateUserDefinedSequence(ind::Individual)
    alphabet = 'A':'Z'
    println("enter a user defined string for plotting: ")
    signalSequenceUserDefined = readline()
    signalSiquence = getUserDefinedSequence(signalSequenceUserDefined, params.silenctInterval, params.letterSize)

    char_sequence = ["$c" for c in alphabet[1:params.noOfSignals]]
    correctIndecies = getCorrectPatternsMarkersABCDEFGHIJ(signalSiquence, char_sequence)

    step = 0
    for j in 1:length(signalSiquence)
        setInput!(ind, signalSiquence[j], j, params.writeNetworkActivity)
        networkStep!(ind, step, params.writeNetworkActivity)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
    signalSiquence = []
    correctIndecies = []
    ind.outputNeurons[1].spikeBitmap = []

    return signalSequenceUserDefined
end

# function reEvaluateTop10_onALargeSequence(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateTop10_onALargeSequence(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEF"
    sequenceLength = 1_000_000
    signalSiquence = getRandomSequenceGivenCorrSignal(signal, sequenceLength, params.silenctInterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
    step = 0

    for s in 1:length(signalSiquence)
        ind.setInput!(signalSiquence[s], s)
        networkStep!(ind, step, params.writeNetworkActivity)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
    println("\tFitness:\t", ind.fitness, "\tTPR:\t", ind.rewardn, "\tFDR:\t", ind.fdr)
    ind.outputNeurons[1].spikeBitmap = []
end

function reEvaluateOnLargeSequence(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair}, params::Parameters)
    # function reEvaluateOnLargeSequence(ind::Individual, pattFrqStructList::Vector, params::Parameters)
    alphabet = 'A':'Z'
    signal = ["$c" for c in alphabet[1:params.noOfSignals]]

    sequenceLength = params.reevaluateSeq
    signalSiquence = getRandomSequenceGivenCorrSignalUpgraded(signal, sequenceLength, params.silenctInterval, params.letterSize, signal)

    if params.noOfSignals == 10
        correctIndecies = getCorrectPatternsMarkersABCDEFGHIJ(signalSiquence, signal)
    elseif params.noOfSignals == 9
        correctIndecies = getCorrectPatternsMarkersABCDEFGHI(signalSiquence, signal)
    elseif params.noOfSignals == 8
        correctIndecies = getCorrectPatternsMarkersABCDEFGH(signalSiquence, signal)
    elseif params.noOfSignals == 7
        correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
    elseif params.noOfSignals == 6
        correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
    elseif params.noOfSignals == 5
        correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
    elseif params.noOfSignals == 4
        correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
    elseif params.noOfSignals == 3
        signalSiquence = getRandomSequenceGivenCorrSignal(signal, sequenceLength, params.silenctInterval, params.letterSize)
        correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
    end

    step = 0
    for s in 1:length(signalSiquence)
        setInput!(ind, signalSiquence[s], s, params.writeNetworkActivity)
        networkStep!(ind, step, params.writeNetworkActivity)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
    println("\tFitness:\t", ind.fitness, "\tTPR:\t", ind.rewardn, "\tFDR:\t", ind.fdr)
    ind.outputNeurons[1].spikeBitmap = []
end

# function reEvaluateTop10ind(top10ind::Vector{Individual}, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateTop10ind(top10ind::Vector{Individual}, pattFrqStructList::Vector)
    signal = "ABCDEF"
    sequenceLength = 500_000
    signalSiquence = Vector{Char}()

    for i in 1:6
        for j in 1:6
            for k in 1:6
                for l in 1:6
                    for m in 1:6
                        for n in 1:6
                            pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
                            patternReconizedCount = 0

                            for nw in top10ind
                                step = 0
                                for s in 1:length(signalSiquence)
                                    nw.setInput!(signalSiquence[s], s)
                                    nw.networkStep(step)
                                    step += 1
                                end

                                bitmapSize = length(nw.outputNeurons[1].spikeBitmap)
                                for sbp in bitmapSize-30:bitmapSize
                                    if nw.outputNeurons[1].spikeBitmap[sbp]
                                        patternReconizedCount += 1
                                        break
                                    end
                                end

                                nw.outputNeurons[1].spikeBitmap = []
                            end

                            if patternReconizedCount > 0
                                # tempPair = PatternFrequencyPair(pattern, patternReconizedCount)
                                tempPair = (pattern, patternReconizedCount)
                                push!(pattFrqStructList, tempPair)
                            end
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
end

function run_pattern_through_network!(signal::String, signalSiquence, pattern::String, ind::Individual, tempPair::PatternFrequencyPair, pattFrqStructList::Vector{PatternFrequencyPair}, params::Parameters)
    patternReconizedCount = 0
    max_reps = 1

    for _ in 1:max_reps
        step = 0
        for s in 1:length(signalSiquence)
            setInput!(ind, signalSiquence[s], s, params.writeNetworkActivity)
            networkStep!(ind, step, params)
            step += 1
        end

        bitmapSize = length(ind.outputNeurons[1].spikeBitmap)
        for sbp in 30:bitmapSize
            if ind.outputNeurons[1].spikeBitmap[sbp]
                patternReconizedCount += 1
                break
            end
        end

        ind.outputNeurons[1].spikeBitmap = []
    end

    if patternReconizedCount > 0
        tempPair.pattern = pattern
        tempPair.freq = patternReconizedCount
        push!(pattFrqStructList, tempPair)

        if occursin(signal, pattern)
            ind.reward += patternReconizedCount
        else
            ind.penalty += patternReconizedCount
        end
    end
end

# function reEvaluateAllPerm9_7sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm9_7sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm9_7sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEFG"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                for p in 1:signal_len
                                    for q in 1:signal_len
                                        tempPair = PatternFrequencyPair()

                                        pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])$(signal[p])$(signal[q])"
                                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                        # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEFG(signalSiquence, signal)
                                        run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                        signalSiquence = []
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm8_7sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm8_7sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm8_7sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEFG"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                for p in 1:signal_len
                                    tempPair = PatternFrequencyPair()

                                    pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])$(signal[p])"
                                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                    # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEFG(signalSiquence, signal)
                                    run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                    signalSiquence = []
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm7_7sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm7_7sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm7_7sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEFG"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                tempPair = PatternFrequencyPair()

                                pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])"
                                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEFG(signalSiquence, signal)
                                run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                signalSiquence = []
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm8_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm8_6sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm8_6sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEF"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                for p in 1:signal_len
                                    tempPair = PatternFrequencyPair()

                                    pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])$(signal[p])"
                                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                    # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEF(signalSiquence, signal)
                                    run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                    signalSiquence = []
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end



# function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm7_6sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEF"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                tempPair = PatternFrequencyPair()

                                pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])"
                                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEF(signalSiquence, signal)
                                run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                signalSiquence = []
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm6_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm6_6sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm6_6sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDEF"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            tempPair = PatternFrequencyPair()

                            pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                            # correctIndecies = getCorrectPatternsMarkernsMarkersABCDEF(signalSiquence, signal)
                            run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm7_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm7_5sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm7_5sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            for o in 1:signal_len
                                tempPair = PatternFrequencyPair()

                                pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])"
                                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                                # correctIndecies = getCorrectPatternsMarkernsMarkersABCDE(signalSiquence, signal)
                                run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                                # signalSiquence = []
                            end
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function run_evaluation(input_pattern, ind, pattFrqStructList, nesting_lvl)

#     if nesting_lvl ==1
#     for n in 1:signal_len
#         pattern = input_pattern*"$(signal[n])"
#         signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
#         # correctIndecies = getCorrectPatternsMarkernsMarkersABCDE(signalSiquence, signal)
#         run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
#         signalSiquence = []
#     end
# else
#     nesting_lvl-=1
# end
# end

# function reEvaluateAllPerm6_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm6_5sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm6_5sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            tempPair = PatternFrequencyPair()

                            pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)
                            # correctIndecies = getCorrectPatternsMarkernsMarkersABCDE(signalSiquence, signal)
                            run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                            # signalSiquence = []
                        end
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm5_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm5_5sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm5_5sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""


    permutation_part_len = 5
        end
    end
    # for i in 1:signal_len
    #     for j in 1:signal_len
    #         for k in 1:signal_len
    #             for l in 1:signal_len
    #                 for m in 1:signal_len
    #                     tempPair = PatternFrequencyPair()

    #                     pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
    #                     signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)

    #                     # Is the line below needed, if its resutls is never used?
    #                     # correctIndecies = getCorrectPatternsMarkernsMarkersABCDE(signalSiquence, signal)
    #                     run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
    #                     # signalSiquence = []
    #                 end
    #             end
    #         end
    #     end
    # end
    return pattFrqStructList
end

# function reEvaluateAllPerm4_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm4_4sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
function reEvaluateAllPerm4_4sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()
    signal = "ABCD"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    tempPair = PatternFrequencyPair()

                    pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])"
                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                    # correctIndecies = getCorrectPatternsMarkernsMarkersABCD(signalSiquence, signal)
                    run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                    # signalSiquence = []
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm5_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair}, params::Parameters)
function reEvaluateAllPerm5_4sig(ind::Individual, params::Parameters)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    # function reEvaluateAllPerm5_4sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCD"
    signal_len = length(signal)
    tempPair = PatternFrequencyPair()
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                        # correctIndecies = getCorrectPatternsMarkernsMarkersABCD(signalSiquence, signal)
                        run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                        signalSiquence = []
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm6_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm6_4sig(ind::Individual, pattFrqStructList::Vector)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABCD"
    signal_len = length(signal)
    tempPair = PatternFrequencyPair()
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                        # correctIndecies = getCorrectPatternsMarkernsMarkersABCD(signalSiquence, signal)
                        run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                        signalSiquence = []
                    end
                end
            end
        end
    end
    return pattFrqStructList
end


# function reEvaluateAllPerm5_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm5_3sig(ind::Individual, pattFrqStructList::Vector)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABC"
    signal_len = length(signal)
    tempPair = PatternFrequencyPair()
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    for m in 1:3
                        pattern = "C$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                        # # correctIndecies = getCorrectPatternsMarkernsMarkersABC(signalSiquence, signal)
                        # correctIndecies = getCorrectPatternsMarkernsMarker(signalSiquence, signal)
                        run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                        signalSiquence = []
                    end
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm4_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm4_3sig(ind::Individual, pattFrqStructList::Vector)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABC"
    signal_len = length(signal)
    tempPair = PatternFrequencyPair()
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    pattern = "C$(signal[i])$(signal[j])$(signal[k])$(signal[l])"
                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                    # # correctIndecies = getCorrectPatternsMarkernsMarkersABC(signalSiquence, signal)
                    # correctIndecies = getCorrectPatternsMarkernsMarkers(signalSiquence, signal)
                    run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                    signalSiquence = []
                end
            end
        end
    end
    return pattFrqStructList
end

# function reEvaluateAllPerm3_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm3_3sig(ind::Individual, pattFrqStructList::Vector)
    pattFrqStructList = Vector{PatternFrequencyPair}()

    signal = "ABC"
    signal_len = length(signal)
    tempPair = PatternFrequencyPair()
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                pattern = "C$(signal[i])$(signal[j])$(signal[k])"
                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize)
                # # correctIndecies = getCorrectPatternsMarkernsMarkersABC(signalSiquence, signal)
                # correctIndecies = getCorrectPatternsMarkernsMarkers(signalSiquence, signal)
                run_pattern_through_network!(signal, signalSiquence, pattern, ind, tempPair, pattFrqStructList, params)
                signalSiquence = []
            end
        end
    end
    return pattFrqStructList
end

function reEvaluateOnProblematicPatternsFound(ind::Individual, probSignalSequence::String, params::Parameters)
    signal = "ABCDEF"
    step = 0
    signalSiquence = getUserDefinedSequence(probSignalSequence, params.silenctInterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)

    for j in 1:length(signalSiquence)
        step += 1
        setInput!(ind, signalSiquence[j], j, params.writeNetworkActivity)
        networkStep!(ind, step, params.writeNetworkActivity)
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
end

function reEvaluateABCD(ind::Individual)
    signal = "ABCD"
    step = 0
    signalSiquence = getRandomSequenceABCD(signal, 500_000, params.silenctInterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)

    for j in 1:length(signalSiquence)
        step += 1
        setInput!(ind, signalSiquence[j], j, params.writeNetworkActivity)
        networkStep!(ind, step, params.writeNetworkActivity)
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
end


function fitness(ind::Individual, signalSiquence::Vector{Char}, correctIndecies::Vector{Int}, params::Parameters)
    ind.rewardn = 0
    ind.reward = 0
    ind.penaltyn = 0
    ind.penalty = 0
    patCount = -1
    reward = 0.0
    penalty = 0.0
    rCountint = 0.0
    rCountSig = 0.0
    pCountint = 0.0
    pCountSig = 0.0
    k = 1
    inCorrPattIndex = Int[]
    sizeHistory = length(ind.outputNeurons[1].spikeBitmap)

    for j in 1:length(correctIndecies)
        while k < sizeHistory && k != correctIndecies[j]
            while signalSiquence[k] != 'Z' && k != correctIndecies[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountSig += 1.0
                end
                k += 1
            end

            while signalSiquence[k] == 'Z' && k != correctIndecies[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountint += 1.0
                end
                k += 1
            end

            patCount += 1
            if pCountint + pCountSig > 0.0
                push!(inCorrPattIndex, k)
                pCountint = 0.0
                pCountSig = 0.0
                penalty += 1.0
            end
        end

        while k < sizeHistory && signalSiquence[k] != 'Z'
            if ind.outputNeurons[1].spikeBitmap[k]
                rCountSig += 1.0
            end
            k += 1
        end

        while k < sizeHistory && signalSiquence[k] == 'Z'
            if ind.outputNeurons[1].spikeBitmap[k]
                rCountint += 1.0
            end
            k += 1
        end

        patCount += 1
        if rCountSig + rCountint > 0.0
            rCountSig = 0.0
            rCountint = 0.0
            reward += 1.0
        end
    end

    repetition_size = (params.letterSize + params.silenctInterval)
    incPatt = ""
    # @info inCorrPattIndex

    """
    NOTE: this was removed because it throws BoundsError- tries to access signal sequence at position -28
    There was an issue with signal index (difference by 1), but the error is still there after this fix.

    variable `missIdentifiedPatterns` is not used here

    Main question is- why 'ii' loop iterates 7 times in total? This works correctly (not throwin an error) in C++ version of the code.
    TODO fix revealing the missclassification of incorrect pattern
    """
    reveal_misidentified_patterns = false
    if reveal_misidentified_patterns
        if params.noOfSignals >= 4
            for icp in 1:length(inCorrPattIndex)
                signalIndex = inCorrPattIndex[icp]

                # TODO why there is the 180 here? in case when the correct sequence is inserted at the beginning, this causes crash
                if signalIndex >= 181 # originally this was 180, but given the difference in indexing, this was incresead +1
                    # @info "Current signal sequence \n$(signalSiquence)"
                    for ii in 1:7
                        signalIndex -= repetition_size
                        incPatt *= signalSiquence[signalIndex+1]
                    end
                end
                incPatt = reverse(incPatt)
                push!(ind.missIdentifiedPatterns, incPatt)
                incPatt = ""
            end
        end
    end # reveal_misidentified_patterns

    ind.totalCorrPatterns = length(correctIndecies)
    ind.rewardn = reward / length(correctIndecies)
    ind.penaltyn = penalty / (params.noOfLetters - length(correctIndecies))
    ind.fdr = penalty / (penalty + reward)
    ind.precision = reward / (penalty + reward)
    ind.reward = reward
    ind.penalty = penalty

    final_reward = 1 - (ind.rewardn - 25 * ind.penaltyn)
    return final_reward
end


end # module
