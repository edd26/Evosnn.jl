module GaModule
using Random

export Ga, run_Ga!

using ..IndividualModule: Individual
using ..UtilityFunctions
using ..ParamsModule
using ..UtilityFunctions

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

    if params.noOfSignals == 10
        signalSiquence = get_abcdefhXXX_XXXdefghij_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDEFGHIJ(signalSiquence, signal)
    elseif params.noOfSignals == 9
        signalSiquence = get_abcdefXXX_XXXdefghi_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDEFGHI(signalSiquence, signal)
    elseif params.noOfSignals == 8
        signalSiquence = get_abcdeXXX_XXXfgh_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDEFGH(signalSiquence, signal)
    elseif params.noOfSignals == 7
        signalSiquence = get_abcdXXX_XXXdefg_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
    elseif params.noOfSignals == 6
        signalSiquence = get_abcXXX_XXXdef_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
    elseif params.noOfSignals == 5
        signalSiquence = getABXXX_XXXDE_Sequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
    elseif params.noOfSignals == 4
        signalSiquence = getABCDSequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
    elseif params.noOfSignals == 3
        signalSiquence = getABCSequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        if !isempty(hardPatternSeq)
            append!(signalSiquence, hardPatternSeq)
        end
        correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
        if params.getAlmostCorrectNWs
            correctIndecies = getCorrectPatternsMarkersAB(signalSiquence, "BC")
        end
    elseif params.noOfSignals == 2
        signalSiquence = getABSequence(signal, params.noOfLetters, params.silenctinterval, params.letterSize)
        correctIndecies = getCorrectPatternsMarkersAB(signalSiquence, signal)
    end

    for i in 1:length(pop)
        for j in 1:length(signalSiquence)
            step += 1
            pop[i].setInput(signalSiquence[j], j)
            pop[i].networkStep(step)
        end

        pop[i].fitness = fitness(pop[i], signalSiquence, correctIndecies)

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
        pop[rep].replicate()
    end

    signalSiquence = []
    correctIndecies = []
end

function runMinimalWeight(ga::Ga, pop::Vector{Individual}, genNo::Int, params::Parameters)
    signal::String = "ABC"
    letterSize::Int = params.letterSize
    gapSize::Int = params.silenctinterval
    sequenceLength::Int = params.noOfLetters
    noOfSeq::Int = 6
    step::Int = 0
    reward::Float64 = 0.0
    penalty::Float64 = 0.0
    rewardAll::Float64 = 0.0
    penaltyAll::Float64 = 0.0
    signalSiquence::Vector{Char} = getRandomSequence(signal, sequenceLength, gapSize, letterSize)
    correctIndecies::Vector{Int} = getCorrectPatternsMarkersABC(signalSiquence, signal)

    for i in 1:length(pop)
        reward = 0.0
        penalty = 0.0
        rewardAll = 0.0
        penaltyAll = 0.0

        for j in 1:length(signalSiquence)
            step += 1
            pop[i].setInput(signalSiquence[j], j)
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
    signalSiquence = getUserDefinedSequence(signalSequenceUserDefined, params.silenctinterval, params.letterSize)

    char_sequence = ["$c" for c in alphabet[1:params.noOfSignals]]
    correctIndecies = getCorrectPatternsMarkersABCDEFGHIJ(signalSiquence, char_sequence)

    step = 0
    for j in 1:length(signalSiquence)
        ind.setInput(signalSiquence[j], j)
        ind.networkStep(step)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies)
    signalSiquence = []
    correctIndecies = []
    ind.outputNeurons[1].spikeBitmap = []

    return signalSequenceUserDefined
end

# function reEvaluateTop10_onALargeSequence(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateTop10_onALargeSequence(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEF"
    sequenceLength = 1_000_000
    signalSiquence = getRandomSequenceGivenCorrSignal(signal, sequenceLength, params.silenctinterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
    step = 0

    for s in 1:length(signalSiquence)
        ind.setInput(signalSiquence[s], s)
        ind.networkStep(step)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies)
    println("\tFitness:\t", ind.fitness, "\tTPR:\t", ind.rewardn, "\tFDR:\t", ind.fdr)
    ind.outputNeurons[1].spikeBitmap = []
end

# function reEvaluateOnLargeSequence(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateOnLargeSequence(ind::Individual, pattFrqStructList::Vector)
    alphabet = 'A':'Z'
    signal = ["$c" for c in alphabet[1:params.noOfSignals]]

    sequenceLength = params.reevaluateSeq
    signalSiquence = getRandomSequenceGivenCorrSignalUpgraded(signal, sequenceLength, params.silenctinterval, params.letterSize, signal)

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
        signalSiquence = getRandomSequenceGivenCorrSignal(signal, sequenceLength, params.silenctinterval, params.letterSize)
        correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
    end

    step = 0
    for s in 1:length(signalSiquence)
        ind.setInput(signalSiquence[s], s)
        ind.networkStep(step)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies)
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
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
                            patternReconizedCount = 0

                            for nw in top10ind
                                step = 0
                                for s in 1:length(signalSiquence)
                                    nw.setInput(signalSiquence[s], s)
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

# function run_pattern_through_network(signal::String, pattern::String, ind::Individual, tempPair::PatternFrequencyPair, pattFrqStructList::Vector{PatternFrequencyPair})
function run_pattern_through_network(signal::String, pattern::String, ind::Individual, tempPair, pattFrqStructList::Vector)
    patternReconizedCount = 0
    max_reps = 1

    for _ in 1:max_reps
        step = 0
        for s in 1:length(signalSiquence)
            ind.setInput(signalSiquence[s], s)
            ind.networkStep(step)
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
function reEvaluateAllPerm9_7sig(ind::Individual, pattFrqStructList::Vector)
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
                                pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])"
                                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                                correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
                                run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                                signalSiquence = []
                            end
                        end
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm8_7sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm8_7sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEFG"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
                            run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm7_7sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm7_7sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEFG"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            pattern = "G$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDEFG(signalSiquence, signal)
                            run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm8_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm8_6sig(ind::Individual, pattFrqStructList::Vector)
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
                                pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])$(signal[o])"
                                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                                correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
                                run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                                signalSiquence = []
                            end
                        end
                    end
                end
            end
        end
    end
end



# function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEF"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm6_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm6_6sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDEF"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        pattern = "F$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm7_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm7_5sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""

    for i in 1:signal_len
        for j in 1:signal_len
            for k in 1:signal_len
                for l in 1:signal_len
                    for m in 1:signal_len
                        for n in 1:signal_len
                            pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
                            run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm6_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm6_5sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""

    for i in 1:5
        for j in 1:5
            for k in 1:5
                for l in 1:5
                    for m in 1:5
                        for n in 1:5
                            pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])$(signal[n])"
                            signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                            correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
                            run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                            signalSiquence = []
                        end
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm5_5sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm5_5sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCDE"
    signal_len = length(signal)
    pattern = ""

    for i in 1:5
        for j in 1:5
            for k in 1:5
                for l in 1:5
                    for m in 1:5
                        pattern = "E$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABCDE(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm4_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm4_4sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCD"
    signal_len = length(signal)
    pattern = ""

    for i in 1:4
        for j in 1:4
            for k in 1:4
                for l in 1:4
                    pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])"
                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                    correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
                    run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                    signalSiquence = []
                end
            end
        end
    end
end

# function reEvaluateAllPerm5_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm5_4sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCD"
    signal_len = length(signal)
    pattern = ""

    for i in 1:4
        for j in 1:4
            for k in 1:4
                for l in 1:4
                    for m in 1:4
                        pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm6_4sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm6_4sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABCD"
    signal_len = length(signal)
    pattern = ""

    for i in 1:4
        for j in 1:4
            for k in 1:4
                for l in 1:4
                    for m in 1:4
                        pattern = "D$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end


# function reEvaluateAllPerm5_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm5_3sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABC"
    signal_len = length(signal)
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    for m in 1:3
                        pattern = "C$(signal[i])$(signal[j])$(signal[k])$(signal[l])$(signal[m])"
                        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                        correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
                        run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                        signalSiquence = []
                    end
                end
            end
        end
    end
end

# function reEvaluateAllPerm4_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm4_3sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABC"
    signal_len = length(signal)
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    pattern = "C$(signal[i])$(signal[j])$(signal[k])$(signal[l])"
                    signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                    correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
                    run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                    signalSiquence = []
                end
            end
        end
    end
end

# function reEvaluateAllPerm3_3sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
function reEvaluateAllPerm3_3sig(ind::Individual, pattFrqStructList::Vector)
    signal = "ABC"
    signal_len = length(signal)
    pattern = ""

    for i in 1:3
        for j in 1:3
            for k in 1:3
                pattern = "C$(signal[i])$(signal[j])$(signal[k])"
                signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctinterval, params.letterSize)
                correctIndecies = getCorrectPatternsMarkersABC(signalSiquence, signal)
                run_pattern_through_network(signal, pattern, ind, tempPair, pattFrqStructList)
                signalSiquence = []
            end
        end
    end
end

function reEvaluateOnProblematicPatternsFound(ind::Individual, probSignalSequence::String)
    signal = "ABCDEF"
    step = 0
    signalSiquence = getUserDefinedSequence(probSignalSequence, params.silenctinterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)

    for j in 1:length(signalSiquence)
        step += 1
        ind.setInput(signalSiquence[j], j)
        ind.networkStep(step)
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies)
end

function reEvaluateABCD(ind::Individual)
    signal = "ABCD"
    step = 0
    signalSiquence = getRandomSequenceABCD(signal, 500_000, params.silenctinterval, params.letterSize)
    correctIndecies = getCorrectPatternsMarkersABCD(signalSiquence, signal)

    for j in 1:length(signalSiquence)
        step += 1
        ind.setInput(signalSiquence[j], j)
        ind.networkStep(step)
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies)
end



end # module
