module UtilityFunctions

export get_abcdefhXXX_XXXdefghij_Sequence

using Random
using Distributions

function getGaussianRandomNumber(de::MersenneTwister, mean::Float64, stddev::Float64)
    return rand(Normal(mean, stddev))
end

function getGaussianValueWithGivenMeanAndSD(mean::Float64, sd::Float64, seqSize::Int64)
    rng = MersenneTwister()
    Random.seed!(rng)
    nd = Normal(mean, sd)
    gaussVector = [rand(nd, rng) for _ in 1:seqSize]
    return gaussVector
end

function getRandomValue(LO::Float64, HI::Float64)
    return rand(Uniform(LO, HI))
end

function randomInputSequence(signals::String, sequenceLength::Int)
    return [rand(1:length(signals)) for _ in 1:sequenceLength]
end

function getGroundtruth(seq::Vector{Char}, corSeq::String)
    correctMarkers = Int[]
    for i in 1:length(seq)
        j = 1
        k = i
        while j <= length(corSeq) && k <= length(seq) && seq[k] == corSeq[j]
            j += 1
            k += 1
            if j > length(corSeq)
                push!(correctMarkers, i + j - 2)
                j = 1
                k = i
                break
            end
        end
    end
    return correctMarkers
end

function getRandomSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    randSequence = ""
    chSequence = "AAAABBBBCCCC"
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function insertGapsAndSetLetterSize(randStr::String, gap::Int, letterSize::Int)
    tempVector = Char[]
    for i in 1:length(randStr)
        varSignal = getUniformVariation(params.variationOnSignal)
        append!(tempVector, fill(randStr[i], letterSize + varSignal))
        varSilence = getUniformVariation(params.variationOnSilence)
        append!(tempVector, fill('Z', gap + varSilence))
    end
    return tempVector
end

















function getUniformVariation(num::Int)
    return rand(0:num)
end

function checkOrCreateDirectory(path::String)
    try
        if !isdir(path)
            mkpath(path)
            println("Created directory: ", path)
            return true
        elseif isdir(path)
            println("Output directory already exists: ", path)
            return false
        end
        throw(ErrorException("Path exists but is not a directory"))
    catch e
        println("Output directory not created: ", path)
        println("Error: ", e)
        return false
    end
end


function getPaternBeforeSpikeIndex(signalSiquence::Vector{Char}, index::Int64, patlength::Int)
    patt = ""
    for ii in 1:patlength
        index -= 30
        if index < 0
            break
        else
            patt *= signalSiquence[index]
        end
    end
    return reverse(patt)
end

function getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    sizePS = length(patSequence)
    for i in 1:sizePS
        j = 1
        while i <= sizePS && patSequence[i] == word[j] && j <= length(word)
            j += 1
            if j > length(word)
                push!(correctPatternsMarkers, i)
                while i <= sizePS && patSequence[i] != 'Z'
                    i += 1
                end
                while i <= sizePS && patSequence[i] == 'Z'
                    i += 1
                end
            end
        end
        while i <= sizePS && !in(patSequence[i], word)
            i += 1
        end
    end
    return correctPatternsMarkers
end

function getCorrectPatternsMarkersABCDEF(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    sizePS = length(patSequence)
    for i in 1:sizePS
        j = 1
        while i <= sizePS && j <= length(word) && patSequence[i] == word[j]
            j += 1
            if j > length(word)
                push!(correctPatternsMarkers, i)
                while i <= sizePS && patSequence[i] != 'Z'
                    i += 1
                end
                while i <= sizePS && patSequence[i] == 'Z'
                    i += 1
                end
            end
        end
    end
    return correctPatternsMarkers
end

function getCorrectPatternsMarkersABCD(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    sizePS = length(patSequence)

    for i in 1:sizePS
        j = 1
        while i <= sizePS && patSequence[i] == word[j] && j <= length(word)
            j += 1
            if j > length(word)
                push!(correctPatternsMarkers, i)
                while i <= sizePS && patSequence[i] != 'Z'
                    i += 1
                end
                while i <= sizePS && patSequence[i] == 'Z'
                    i += 1
                end
            end
        end
    end
    return correctPatternsMarkers
end

function generatePermutationsWithReplacement(str::String, current::String, length::Int, permutations::Vector{String}, prefix::String, suffix::String)
    if length(current) == length
        push!(permutations, prefix * current * suffix)
        return
    end
    for c in str
        generatePermutationsWithReplacement(str, current * c, length, permutations, prefix, suffix)
    end
end

function getRandomSequenceGivenCorrSignal(chSequence::String, size::Int, gap::Int, letterSize::Int, sequence::String)
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    insertSequenceIntoLetterChain(sequence, randSequence, 100)
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function getPermutations(sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String}, sequence::String)
    total_sequences = length(sequence_heads)
    all_permutations = Vector{Vector{String}}(undef, total_sequences)

    for i in 1:total_sequences
        some_vector = String[]
        generatePermutationsWithReplacement(permutation_bases[i], "", 3, some_vector, sequence_heads[i], sequence_tails[i])
        all_permutations[i] = some_vector
    end
    return all_permutations
end

function generateSequenceWithPermutation(size::Int, gap::Int, letterSize::Int, sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String}, sequence::String)
    total_sequences = length(sequence_heads)
    all_permutations = getPermutations(sequence_heads, permutation_bases, sequence_tails, sequence)

    totalSize = sum(length.(all_permutations))
    chSequence = String[]

    index = 1
    for perms in all_permutations
        for s in perms
            chSequence[index] = s
            index += 1
        end
    end

    randSequence = ""
    for _ in 1:size
        chosenSeq = chSequence[rand(1:totalSize)]
        randSequence *= chosenSeq
    end

    insertSequenceIntoLetterChain(sequence, randSequence, 100)
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function getRandomSequenceGivenCorrSignalUpgraded(chSequence::String, size::Int, gap::Int, letterSize::Int, sequence::String)
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    insertSequenceIntoLetterChain(sequence, randSequence, 100)
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function get_abcdXXX_XXXdefg_Sequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    # srand(time(0))
    sequence = "ABCDEFG"
    sequence_heads = ["ABCD", "", "AB", "ABC", "", ""]
    permutation_bases = ["EFG", "ABC", "CDE", "DEF", "ABC"]
    sequence_tails = ["", "DEFG", "FG", "G", "EDEF", "CDEF"]
    return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence)
end

function getABCDEFSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    # srand(time(0))
    randSequence = ""
    chSequence = ["ABCAAF", "ABCABF", "ABCACF", "ABCADF", "ABCAEF", "ABCBAF", "ABCBBF", "ABCBCF", "ABCBDF", "ABCBEF", "ABCCAF", "ABCCBF", "ABCCCF", "ABCCDF", "ABCCEF", "ABCDAF", "ABCDBF", "ABCDCF", "ABCDDF", "ABCDEF", "ABCEAF", "ABCEBF", "ABCECF", "ABCEDF", "ABCEEF"]
    arraySize = length(chSequence)

    for _ in 1:size
        chosenChar = chSequence[rand(1:arraySize)]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function getUserDefinedSequence(inputStr::String, gap::Int, letterSize::Int)
    randSequence = inputStr
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end


function getCorrectPatternsMarkersABC(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    sizePS = length(patSequence)

    for i in 1:sizePS
        j = 1
        while i <= sizePS && patSequence[i] == word[j] && j <= length(word)
            j += 1
            if j > length(word)
                push!(correctPatternsMarkers, i)
                while i <= sizePS && patSequence[i] != 'Z'
                    i += 1
                end
                while i <= sizePS && patSequence[i] == 'Z'
                    i += 1
                end
            end
        end
    end
    return correctPatternsMarkers
end

function getCorrectPatternsMarkersAB(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    sizePS = length(patSequence)

    for i in 1:sizePS
        j = 1
        while i <= sizePS && patSequence[i] == word[j] && j <= length(word)
            j += 1
            if j > length(word)
                push!(correctPatternsMarkers, i)
                while i <= sizePS && patSequence[i] != 'Z'
                    i += 1
                end
                while i <= sizePS && patSequence[i] == 'Z'
                    i += 1
                end
            end
        end
    end
    return correctPatternsMarkers
end

function getGroundtruthMarkers(signalSequence::Vector{Char}, corrTrip::String)
    correctMarkers = Int[]
    sizePS = length(signalSequence)
    word = corrTrip

    for i in 1:sizePS
        j = 1
        while j <= length(word) && i <= sizePS
            if signalSequence[i] == word[j]
                j += 1
                if j > length(word)
                    push!(correctMarkers, i - length(word) + 1)
                end
            end
            i += 1
        end
    end
    return correctMarkers
end

function getGaussianValue(size::Int, mean::Float64, stddev::Float64)
    return getGaussianValueWithGivenMeanAndSD(mean, stddev, size)
end

function getPermutedSequence(fileName::String, str::String, k::Int, gap::Int, lettersize::Int)
    open(fileName, "r") do file
        allPerm = read(file, String)
        return insertGapsAndSetLetterSize(allPerm, gap, lettersize)
    end
end


function allPossibleCombinations(seq::String, v::Vector{Int}, pos::Vector{Int}, n::Int)
    if n == length(v)
        for i in 1:n
            seq *= v[pos[i]]
        end
        println(seq)
        return
    end

    for i in 1:length(v)
        pos[n+1] = i
        allPossibleCombinations(seq, v, pos, n + 1)
    end
end

function getAllPossibleCombinationSequence(gap::Int, letterSize::Int)
    allPatterns = ""
    signals = "ABCDEF"
    pos = zeros(Int, length(signals))
    allPossibleCombinations(allPatterns, signals, pos, 0)
    return insertGapsAndSetLetterSize(allPatterns, gap, letterSize)
end

function get_abcdefhXXX_XXXdefghij_Sequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    # srand(time(0))
    sequence = "ABCDEFGHIJ"
    sequence_heads = ["ABCDEFG", "", "ABC", "ABCD"]
    permutation_bases = ["HIJ", "ABC", "DEF", "EFG"]
    sequence_tails = ["", "DEFGHIJ", "GHIJ", "HIJ"]
    return generateSequenceWithPermutation_with_doubling(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence)
end


function get_abcdefXXX_XXXdefghi_Sequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    # srand(time(0))
    sequence = "ABCDEFGHI"
    sequence_heads = ["ABCDEF", "", "ABC"]
    permutation_bases = ["GHI", "ABC", "DEF"]
    sequence_tails = ["", "DEFGHI", "GHI"]
    return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence)
end

function generateSequenceWithPermutation_with_doubling(size::Int, gap::Int, letterSize::Int, sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String}, sequence::String)
    total_sequences = length(sequence_heads)
    all_permutations = getPermutations(sequence_heads, permutation_bases, sequence_tails, sequence)

    totalSize = sum(length.(all_permutations))
    chSequence = String[]

    index = 1
    for perms in all_permutations
        for s in perms
            chSequence[index] = s
            index += 1
        end
    end

    randSequence = ""
    for _ in 1:size
        chosenSeq = chSequence[rand(1:totalSize)]
        randSequence *= chosenSeq
    end

    insertSequenceIntoLetterChain(sequence, randSequence, 100)
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end


function get_abcdeXXX_XXXfgh_Sequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int)
    # srand(time(0))
    sequence = "ABCDEFG"
    sequence_heads = ["ABCDE", "", "ABC"]
    permutation_bases = ["FGH", "ABC", "DEF"]
    sequence_tails = ["", "DEFGH", "GH"]
    return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence)
end

function getPermutationWithReplacement(sequence::String, size::Int, gap::Int, letterSize::Int)
    totalPermutations = Vector{String}()
    generatePermutationsWithReplacement(sequence, "", size, totalPermutations, "", "")
    randSequence = totalPermutations[rand(1:length(totalPermutations))]
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

function getRandomCharacterSequence(size::Int, gap::Int, letterSize::Int)
    chSequence = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize)
end

# function insertSequenceIntoLetterChain(sequence::String, randSequence::String, subSequenceLen::Int)
#     seqLength = length(sequence)
#     positionInSubSequence = 0
#
#     for i in 1:(length(randSequence)-subSequenceLen)
#         positionInSubSequence = rand(1:(subSequenceLen-seqLength-1))
#         for j in 1:seqLength
#             subSequenceLocation = i + positionInSubSequence + j - 1
#             randSequence[subSequenceLocation] = sequence[j]
#         end
#     end
# end

function insertSequenceIntoLetterChain(phrase::String, phrases_with_dobuling_vec::Vector{String}, randSequence::String, subSequenceLen::Int)
    half_subSequenceLen = div(subSequenceLen, 2)
    phrase_len = length(phrase)
    phrase_with_dobuling = phrases_with_dobuling_vec[1]
    phrase_with_doubling_len = length(phrase_with_dobuling)
    total_phrases = length(phrases_with_dobuling_vec)

    # for i in 1:sub(length(randSequence), subSequenceLen)
    for i in 1:(length(randSequence)-subSequenceLen)
        positionInSubSequence = rand(1:(half_subSequenceLen-phrase_len-1))
        position_in_sub_sequence_for_doubled_letter = rand(1:(half_subSequenceLen-phrase_with_doubling_len-2)) + half_subSequenceLen

        for j in 1:phrase_len
            subSequenceLocation = i + positionInSubSequence + j - 1
            randSequence[subSequenceLocation] = phrase[j]
        end

        for k in 1:phrase_with_doubling_len
            subSequenceLocation2 = i + position_in_sub_sequence_for_doubled_letter + k - 1
            randSequence[subSequenceLocation2] = phrase_with_dobuling[k]
        end
    end
end

# TODO decide which one is the correct one
# function insertSequenceIntoLetterChain(sequence::String, placeholders::Vector{String}, randSequence::String, subSequenceLen::Int)
#     seqLength = length(placeholders)
#
#     for i in 1:(length(randSequence)-subSequenceLen)
#         position = rand(1:(length(randSequence)-seqLength))
#         for j in 1:seqLength
#             randSequence[position+j-1] = placeholders[j]
#         end
#     end
# end


function generateSequenceWithReplacement(sequence::String, size::Int, gap::Int, letterSize::Int)
    generatedSequence = ""
    for _ in 1:size
        generatedSequence *= sequence[rand(1:length(sequence))]
    end
    return insertGapsAndSetLetterSize(generatedSequence, gap, letterSize)
end

function randomizeInputSequence(signals::String, sequenceLength::Int)
    seq = []
    for _ in 1:sequenceLength
        push!(seq, signals[rand(1:length(signals))])
    end
    return seq
end

function createFileIfNotExists(fileName::String)
    if !isfile(fileName)
        open(fileName, "w") do file
            write(file, "")
        end
        println("File created: ", fileName)
    else
        println("File already exists: ", fileName)
    end
end

function getTimestampedFileName(baseName::String)
    timestamp = string(now())
    timestamp = replace(timestamp, [" " => "_", ":" => "-", "." => "_"])
    return string(baseName, "_", timestamp, ".txt")
end

function logMessageToFile(fileName::String, message::String)
    open(fileName, "a") do file
        write(file, message * "\n")
    end
end

end #module
