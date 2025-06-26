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
    # noOfSeq::Int = 6
    step::Int = 0

    sequence_generation_func =
        if params.noOfSignals == 10
            get_abcdefgXXX_XXXdefghij_Sequence
        elseif params.noOfSignals == 9
            get_abcdefXXX_XXXdefghi_Sequence
        elseif params.noOfSignals == 8
            get_abcdeXXX_XXXfgh_Sequence
        elseif params.noOfSignals == 7
            get_abcdXXX_XXXdefg_Sequence
        elseif params.noOfSignals == 6
            get_abcXXX_XXXdef_Sequence
        elseif params.noOfSignals == 5
            getABXXX_XXXDE_Sequence
        elseif params.noOfSignals == 4
            getABCDSequence
        elseif params.noOfSignals == 3
            getABCSequence
        elseif params.noOfSignals == 2
            getABSequence
        end
    insertionWindowSize = 20
    randSequence = sequence_generation_func(signal, params.noOfLetters, params.letterSize)
    new_sequence, all_insertions = insertSequenceIntoLetterChain(signal, randSequence, insertionWindowSize)
    signalSequence = "B" * new_sequence[2:end]
    expanded_sequence = insertGapsAndSetLetterSize(
        signalSequence,
        params.silenctInterval,
        params.letterSize,
        params.variationOnSignal,
        params.variationOnSilence
    )
    if !isempty(hardPatternSeq)
        append!(expanded_sequence, hardPatternSeq)
    end

    correctIndices = getCorrectPatternsMarkers(expanded_sequence, signal)
    for i in 1:length(pop)
        # @info "Started running population element $(i)"
        total_elements_in_signal = length(expanded_sequence)
        # @info "Total elements in the sequence $(total_elements_in_signal )"
        for j in 1:total_elements_in_signal
            # j % 100 == 0 && @info "At iteration $(j)"
            step += 1
            setInput!(pop[i], expanded_sequence[j], j, params.writeNetworkActivity)
            networkStep!(pop[i], step, params)
        end

        pop[i].fitness = fitness(pop[i], expanded_sequence, correctIndices, params)
        # ind = pop[i]
        # pop[i].fitness = fitness_simplified(pop[i], signalSequence, correctIndices, params)

        # The code below won't be needed in a parallelised process, because pop elements will be copied
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

    # sort!(pop, by=x -> x.fitness, rev=true)
    sort!(pop, by=x -> x.fitness)

    @info "Fitness of i=1:\t$(pop[1].fitness)"
    @info "Fitness of i=2:\t$(pop[2].fitness)"
    @info "Fitness of i=3:\t$(pop[3].fitness)"
    @info "..."
    @info "Fitness of i=end-2:\t$(pop[end-2].fitness)"
    @info "Fitness of i=end-1:\t$(pop[end-1].fitness)"
    @info "Fitness of i=end:\t$(pop[end].fitness)"

    starting_index = params.eliteCount + 1
    # Binary tournament selection
    for k in starting_index:length(pop)
        firstRand = rand(starting_index:length(pop))
        secondRand = rand(starting_index:length(pop))
        # @info firstRand secondRand
        if firstRand == secondRand
            continue
        elseif pop[firstRand].fitness < pop[secondRand].fitness
            pop[k] = pop[firstRand]
        else
            pop[k] = pop[secondRand]
        end
    end

    for rep in params.eliteCount+1:length(pop)
        pop[rep] = replicate(pop[rep], params)
    end
end


function run_Ga_parallel!(ga::Ga, pop::Vector{Individual}, genNo::Int, hardPatternSeq::Vector{Char}, params::Parameters)
    alphabet::String = "ABCDEFGHIJKLMNOPQRSTUWXYZ"
    signal::String = alphabet[1:params.noOfSignals]
    # noOfSeq::Int = 6
    # step::Int = 0

    sequence_generation_func =
        if params.noOfSignals == 10
            get_abcdefgXXX_XXXdefghij_Sequence
        elseif params.noOfSignals == 9
            get_abcdefXXX_XXXdefghi_Sequence
        elseif params.noOfSignals == 8
            get_abcdeXXX_XXXfgh_Sequence
        elseif params.noOfSignals == 7
            get_abcdXXX_XXXdefg_Sequence
        elseif params.noOfSignals == 6
            get_abcXXX_XXXdef_Sequence
        elseif params.noOfSignals == 5
            getABXXX_XXXDE_Sequence
        elseif params.noOfSignals == 4
            getABCDSequence
        elseif params.noOfSignals == 3
            getABCSequence
        elseif params.noOfSignals == 2
            getABSequence
        end

    insertionWindowSize = 20
    randSequence = sequence_generation_func(signal, params.noOfLetters, params.letterSize)
    new_sequence, all_insertions = insertSequenceIntoLetterChain(signal, randSequence, insertionWindowSize)

    # The following is for debugging >>>
    pattern = "CABABCDEFC"
    new_sequence = pattern * new_sequence[(length(signal)+1):end]
    # The following is for debugging <<<
    signalSequence = "B" * new_sequence[2:end]

    expanded_sequence = insertGapsAndSetLetterSize(
        signalSequence,
        params.silenctInterval,
        params.letterSize,
        params.variationOnSignal,
        params.variationOnSilence
    )
    if !isempty(hardPatternSeq)
        append!(expanded_sequence, hardPatternSeq)
    end

    correctIndices = getCorrectPatternsMarkers(expanded_sequence, signal)

    # TODO: this is where we could parallelise code execution- elements of the population should be independent
    nthreads = Threads.nthreads()
    index_fitness_pairs = [Tuple[] for _ in 1:nthreads]

    # TODO: create a list per each thread, where a paris (id-> fintess) will be stored

    # Threads.@threads for i in 1:length(pop)
    # local_ind = deepcopy(pop[i])
    enumerated_population = [(i, local_ind) for (i, local_ind) in enumerate(pop)]
    Threads.@threads for p in enumerated_population
        i, local_ind = p
        tid = Threads.threadid()
        local_expanded_sequence = copy(expanded_sequence)
        local_correct_indices = copy(correctIndices)

        # @info "Started running population element $(i)"
        total_elements_in_signal = length(local_expanded_sequence)
        # @info "Total elements in the sequence $(total_elements_in_signal )"

        # step = 0
        for j in 1:total_elements_in_signal
            # j % 100 == 0 && @info "At iteration $(j)"
            # step += 1
            step = j
            setInput!(local_ind, local_expanded_sequence[j], j, params.writeNetworkActivity)
            networkStep!(local_ind, step, params)
        end

        list_in_thread = index_fitness_pairs[tid]
        population_fitenss = fitness(local_ind, expanded_sequence, local_correct_indices, params)
        push!(list_in_thread, (i, population_fitenss))
        local_ind.fitness = population_fitenss

        # ind = pop[i]
        # pop[i].fitness = fitness_simplified(pop[i], signalSequence, correctIndices, params)

        # The code below won't be needed in a parallelised process, because pop elements will be copied
        local_ind.outputNeurons[1].spikeBitmap = Bool[]
        local_ind.outputNeurons[1].voltageBuffer = Float64[]
        for p in 1:length(local_ind.inputNeurons)
            local_ind.inputNeurons[p].voltageBuffer = Float64[]
        end
        for q in 1:length(local_ind.interNeurons)
            local_ind.interNeurons[q].voltageBuffer = Float64[]
        end
    end
    merged_index_fitness_pairs = reduce(vcat, index_fitness_pairs)
    # println("\n", merged_index_fitness_pairs, "\n")

    for (index, f) = merged_index_fitness_pairs
        pop[index].fitness = f
    end

    # population_sorting_by_fitness = sortperm(merged_index_fitness_pairs, by=x -> x[2])
    # pop = pop[population_sorting_by_fitness]

    # sort!(pop, by=x -> x.fitness, rev=true)
    sort!(pop, by=x -> x.fitness)

    @info "Fitness of i=1:\t$(pop[1].fitness)"
    @info "Fitness of i=2:\t$(pop[2].fitness)"
    @info "Fitness of i=3:\t$(pop[3].fitness)"
    @info "..."
    @info "Fitness of i=end-2:\t$(pop[end-2].fitness)"
    @info "Fitness of i=end-1:\t$(pop[end-1].fitness)"
    @info "Fitness of i=end:\t$(pop[end].fitness)"

    starting_index = params.eliteCount + 1
    for k in starting_index:length(pop)
        firstRand = rand(starting_index:length(pop))
        secondRand = rand(starting_index:length(pop))
        # @info firstRand secondRand
        if firstRand == secondRand
            continue
        elseif pop[firstRand].fitness < pop[secondRand].fitness
            pop[k] = pop[firstRand]
        else
            pop[k] = pop[secondRand]
        end
    end

    for rep in params.eliteCount+1:length(pop)
        pop[rep] = replicate(pop[rep], params)
    end
    return pop
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

function reEvaluateUserDefinedSequence(ga::Ga, ind::Individual, params::Parameters, signalSequenceUserDefined::String)
    # randSequence = getUserDefinedSequence(signalSequenceUserDefined)
    # signalSiquence = insertGapsAndSetLetterSize(signalSequenceUserDefined, gap, letterSize, variationOnSignal, variationOnSilence)
    signalSiquence = insertGapsAndSetLetterSize(signalSequenceUserDefined, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)


    char_sequence = get_signal_of_len(params.noOfSignals)
    correctIndecies = getCorrectPatternsMarkers(signalSiquence, char_sequence)

    step = 0
    for j in 1:length(signalSiquence)
        setInput!(ind, signalSiquence[j], j, params.writeNetworkActivity)
        networkStep!(ind, step, params)
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
    correctIndecies = getCorrectPatternsMarkers(signalSiquence, signal)
    step = 0

    for s in 1:length(signalSiquence)
        ind.setInput!(signalSiquence[s], s)
        networkStep!(ind, step, params)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
    println("\tFitness:\t", ind.fitness, "\tTPR:\t", ind.rewardn, "\tFDR:\t", ind.fdr)
    ind.outputNeurons[1].spikeBitmap = []
end

function reEvaluateOnLargeSequence(ga::Ga, ind::Individual, params::Parameters)
    # function reEvaluateOnLargeSequence(ind::Individual, pattFrqStructList::Vector, params::Parameters)
    pattFrqStructList = PatternFrequencyPair[]
    alphabet = 'A':'Z'
    signal = ["$c" for c in alphabet[1:params.noOfSignals]]
    # TODO possible replacement of above with the following
    # char_sequence = get_signal_of_len(params.noOfSignals)

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
        networkStep!(ind, step, params)
        step += 1
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
    println("\tFitness:\t", ind.fitness, "\tTPR:\t", ind.rewardn, "\tFDR:\t", ind.fdr)
    ind.outputNeurons[1].spikeBitmap = []

    return pattFrqStructList
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

function run_pattern_through_network!(signalSiquence, ind::Individual, params::Parameters; max_reps::Int=1)
    patternReconizedCount = 0

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
    return patternReconizedCount
end

function evaluateRecognisedPatterns!(patternReconizedCount::Int, signal::String, pattern::String, ind::Individual, tempPair::PatternFrequencyPair)
    pattFrqStructList = Vector{PatternFrequencyPair}()

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
    return pattFrqStructList
end

@deprecate reEvaluateAllPerm9_7sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 9, 7)
@deprecate reEvaluateAllPerm8_7sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 8, 7)
@deprecate reEvaluateAllPerm7_7sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 7, 7)
@deprecate reEvaluateAllPerm8_6sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 8, 6)
# function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector{PatternFrequencyPair})
# function reEvaluateAllPerm7_6sig(ind::Individual, pattFrqStructList::Vector, params::Parameters)
@deprecate reEvaluateAllPerm7_6sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 7, 6)
@deprecate reEvaluateAllPerm6_6sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 6, 6)
@deprecate reEvaluateAllPerm7_5sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 7, 5)
@deprecate reEvaluateAllPerm6_5sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 6, 5)
@deprecate reEvaluateAllPerm5_5sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 5, 5)
@deprecate reEvaluateAllPerm4_4sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 4, 4)
@deprecate reEvaluateAllPerm5_4sig(ind::Individual, params::Parameters) reEvaluateAllPerm(ind, params, 5, 4)
@deprecate reEvaluateAllPerm6_4sig(ind::Individual, pattFrqStructList::Vector) reEvaluateAllPerm(ind, params, 6, 4)
@deprecate reEvaluateAllPerm5_3sig(ind::Individual, pattFrqStructList::Vector) reEvaluateAllPerm(ind, params, 5, 3)
@deprecate reEvaluateAllPerm4_3sig(ind::Individual, pattFrqStructList::Vector) reEvaluateAllPerm(ind, params, 4, 3)
@deprecate reEvaluateAllPerm3_3sig(ind::Individual, pattFrqStructList::Vector) reEvaluateAllPerm(ind, params, 3, 3)

function reEvaluateAllPerm(ind::Individual, params::Parameters, permutation_part_len::Int, signal_len::Int)
    signal = get_signal_of_len(signal_len)

    nthreads = Threads.nthreads()
    local_patterns_lists = [PatternFrequencyPair[] for _ in 1:nthreads]
    local_reward_lists = [0 for _ in 1:nthreads]
    local_penalty_lists = [0 for _ in 1:nthreads]

    all_indices = collect(Iterators.product([1:signal_len for _ in 1:permutation_part_len]...))
    # Threads.@threads for idxs in Iterators.product([1:signal_len for _ in 1:permutation_part_len]...)
    Threads.@threads for i in 1:length(all_indices)
        idxs = all_indices[i]

        tid = Threads.threadid()
        local_ind = deepcopy(ind)

        pattern = signal[signal_len]
        for z in idxs
            pattern *= signal[z]
        end
        # @info pattern
        if min(idxs...) == max(idxs...)
            println(pattern)
        end
        signalSiquence = insertGapsAndSetLetterSize(pattern, params.silenctInterval, params.letterSize, params.variationOnSignal, params.variationOnSilence)

        tempPair = PatternFrequencyPair()
        patternReconizedCount = run_pattern_through_network!(signalSiquence, local_ind, params)
        local_pattFrqStructList = evaluateRecognisedPatterns!(patternReconizedCount, signal, pattern, local_ind, tempPair)

        list_in_thread = local_patterns_lists[tid]
        for p in local_pattFrqStructList
            push!(list_in_thread, p)
        end

        local_reward_lists[tid] += local_ind.reward
        local_penalty_lists[tid] += local_ind.penalty
    end

    # pattFrqStructList = Vector{PatternFrequencyPair}()
    pattFrqStructList = reduce(vcat, local_patterns_lists)

    ind.reward = sum(local_reward_lists)
    ind.penalty = sum(local_penalty_lists)

    return pattFrqStructList
end


function reEvaluateOnProblematicPatternsFound(ind::Individual, probSignalSequence::String, params::Parameters)
    signal = "ABCDEF"
    step = 0
    # randSequence = getUserDefinedSequence(probSignalSequence)
    signalSiquene = insertGapsAndSetLetterSize(probSignalSequence, gap, letterSize, variationOnSignal, variationOnSilence)
    correctIndecies = getCorrectPatternsMarkersABCDEF(signalSiquence, signal)

    for j in 1:length(signalSiquence)
        step += 1
        setInput!(ind, signalSiquence[j], j, params.writeNetworkActivity)
        networkStep!(ind, step, params)
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
        networkStep!(ind, step, params)
    end

    ind.fitness = fitness(ind, signalSiquence, correctIndecies, params)
end


function fitness_simplified(
    ind::Individual,
    signalSequence::Vector{Char},
    correctIndices::Vector{Int},
    params::Parameters
)
    ind.rewardn = 0.0
    ind.reward = 0.0
    ind.penaltyn = 0.0
    ind.penalty = 0.0

    patCount = -1
    reward = 0.0
    penalty = 0.0
    rCountint = 0.0
    rCountSig = 0.0
    pCountint = 0.0
    pCountSig = 0.0

    k = 1
    inCorrPattIndex = Int[]

    sizeHistory = length(
        ind.outputNeurons[1].spikeBitmap
    )

    j = 1







    for j in 1:length(correctIndices)
        while k <= sizeHistory && k != correctIndices[j]
            while k <= sizeHistory && signalSequence[k] != 'Z' && k != correctIndices[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountSig += 1.0
                end
                k += 1
            end
            while k <= sizeHistory && signalSequence[k] == 'Z' && k != correctIndices[j]
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

        while k <= sizeHistory && signalSequence[k] != 'Z'
            if ind.outputNeurons[1].spikeBitmap[k]
                rCountSig += 1.0
            end
            k += 1
        end
        while k <= sizeHistory && signalSequence[k] == 'Z'
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

    if params.noOfSignals >= 4
        for signalIndex in inCorrPattIndex
            if signalIndex >= 180
                incPatt = ""
                for ii in 1:7
                    newIndex = signalIndex - 30
                    if newIndex >= 1
                        incPatt *= signalSequence[newIndex]
                        signalIndex = newIndex
                    else
                        break # Avoid accessing out-of-bounds
                    end
                end
                # ind.missIdentifiedPatterns = get(ind, :missIdentifiedPatterns, String[])
                push!(ind.missIdentifiedPatterns, reverse(incPatt))
            end
        end
        inCorrPattIndex = Int[]
    end
    ind.totalCorrPatterns = length(correctIndices)
    ind.reward = reward
    ind.rewardn = reward / length(correctIndices)
    ind.penalty = penalty
    ind.penaltyn = penalty / (params.noOfLetters - length(correctIndices))
    ind.fdr = penalty / (penalty + reward)
    ind.precision = reward / (penalty + reward)

    final_reward = 1.0 - (ind.rewardn - 25.0 * ind.penaltyn)
    return final_reward
end

function fitness(ind::Individual, signalSequence::Vector{Char}, correctIndices::Vector{Int}, params::Parameters)
    ind.rewardn = 0.0
    ind.reward = 0.0
    ind.penaltyn = 0.0
    ind.penalty = 0.0

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

    for j in 1:length(correctIndices)
        while k <= sizeHistory && k != correctIndices[j]
            while k <= sizeHistory && signalSequence[k] != 'Z' && k != correctIndices[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountSig += 1.0
                end
                k += 1
            end
            while k <= sizeHistory && signalSequence[k] == 'Z' && k != correctIndices[j]
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

        while k <= sizeHistory && signalSequence[k] != 'Z'
            if ind.outputNeurons[1].spikeBitmap[k]
                rCountSig += 1.0
            end
            k += 1
        end
        while k <= sizeHistory && signalSequence[k] == 'Z'
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

    if params.noOfSignals >= 4
        for signalIndex in inCorrPattIndex
            if signalIndex >= 180
                incPatt = ""
                for ii in 1:7
                    newIndex = signalIndex - 30
                    if newIndex >= 1
                        incPatt *= signalSequence[newIndex]
                        signalIndex = newIndex
                    else
                        break # Avoid accessing out-of-bounds
                    end
                end
                # ind.missIdentifiedPatterns = get(ind, :missIdentifiedPatterns, String[])
                push!(ind.missIdentifiedPatterns, reverse(incPatt))
            end
        end
        inCorrPattIndex = Int[]
    end
    ind.totalCorrPatterns = length(correctIndices)
    ind.reward = reward
    ind.rewardn = reward / length(correctIndices)
    ind.penalty = penalty
    ind.penaltyn = penalty / (params.noOfLetters - length(correctIndices))
    ind.fdr = penalty / (penalty + reward)
    ind.precision = reward / (penalty + reward)


    penalty_cooef = 10.0
    final_reward = 1.0 - (ind.rewardn - penalty_cooef * ind.penaltyn)
    return final_reward
end


function fitness_might_not_work(ind::Individual, signalSiquence::Vector{Char}, correctIndecies::Vector{Int}, params::Parameters)
    ind.rewardn = 0
    ind.reward = 0
    ind.penaltyn = 0
    ind.penalty = 0
    patCount = -1
    reward = 0.0
    penalty = 0.0
    rCountint = 0
    rCountSig = 0
    pCountint = 0
    pCountSig = 0
    k = 1
    inCorrPattIndex = Int[]
    sizeHistory = length(ind.outputNeurons[1].spikeBitmap)

    for j in 1:length(correctIndecies)
        while k < sizeHistory && k != correctIndecies[j]
            while signalSiquence[k] != 'Z' && k != correctIndecies[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountSig += 1
                end
                k += 1
            end

            while signalSiquence[k] == 'Z' && k != correctIndecies[j]
                if ind.outputNeurons[1].spikeBitmap[k]
                    pCountint += 1
                end
                k += 1
            end

            patCount += 1
            if pCountint + pCountSig > 0
                push!(inCorrPattIndex, k)
                pCountint = 0
                pCountSig = 0
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
                rCountint += 1
            end
            k += 1
        end

        patCount += 1
        if rCountSig + rCountint > 0
            rCountSig = 0
            rCountint = 0
            reward += 1
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
