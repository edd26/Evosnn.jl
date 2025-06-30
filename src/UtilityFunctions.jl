module UtilityFunctions

export get_abcdefgXXX_XXXdefghij_Sequence,
    # getCorrectPatternsMarkersABCDEFGHIJ,
    # getCorrectPatternsMarkersABCDEFGHI,
    # getCorrectPatternsMarkersABCDEFGH,
    # getCorrectPatternsMarkersABCDEFG,
    # getCorrectPatternsMarkersABCDEF,
    # getCorrectPatternsMarkersABCDE,
    # getCorrectPatternsMarkersABCD,
    # getCorrectPatternsMarkersABC,
    # getCorrectPatternsMarkersAB,
    # get_abcdefgXXX_XXXdefghij_Sequence,
    get_abcdefXXX_XXXdefghi_Sequence,
    get_abcdeXXX_XXXfgh_Sequence,
    get_abcdXXX_XXXdefg_Sequence,
    get_abcXXX_XXXdef_Sequence,
    getABXXX_XXXDE_Sequence,
    getCorrectPatternsMarkers,
    getABCDSequence,
    getABCSequence,
    getRandomValue,
    checkOrCreateDirectory,
    insertSequenceIntoLetterChain,
    insertGapsAndSetLetterSize,
    get_signal_of_len,
    getGaussianValueWithGivenMeanAndSD


using Random
using Distributions
using Combinatorics: combinations

Random.seed!(123)

function getGaussianRandomNumber(de::MersenneTwister, mean::Float64, stddev::Float64)
    return rand(Normal(mean, stddev))
end

function getGaussianValueWithGivenMeanAndSD(mean::Float64, sd::Float64, seqSize::Int64)
    rng = MersenneTwister()
    Random.seed!(rng)
    # nd = Normal(mean, sd)
    # gaussVector = [rand(nd, rng) for _ in 1:seqSize]
    gaussVector = rand(seqSize) * sd .+ mean
    return gaussVector
end

function getRandomValue(LO::T, HI::T) where {T<:Number}
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

function getRandomSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    randSequence = ""
    chSequence = "AAAABBBBCCCC"
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function insertGapsAndSetLetterSize(randStr::String, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    # Initialize result array with reasonable size
    tempVector = Char[]

    # Loop through each character in the string
    # TODO This can be optimized by creating an empty Char vector of predefined size
    # TODO Second stage would be to replace the inner for loops with assigning a range of values from 
    # vector to the rand sequence
    for i in 1:length(randStr)
        # Get signal variation
        varSignal = getUniformVariation(variationOnSignal)

        # Add repeated characters with signal variation
        for _ in 1:(letterSize+varSignal)
            push!(tempVector, randStr[i])
        end

        # Get silence variation
        varSilence = getUniformVariation(variationOnSilence)

        # Add gap characters with silence variation
        for _ in 1:(gap+varSilence)
            push!(tempVector, 'Z')
        end
    end

    return tempVector
end

function getUniformVariation(num::Int)
    return rand(0:num)
end

function checkOrCreateDirectory(path::String)
    try
        if !ispath(path)
            mkpath(path)
            println("Created path: ", path)
            return true
        elseif ispath(path)
            println("Output path already exists: ", path)
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

# getCorrectPatternsMarkersABC(patSequence::Vector{Char}, corrTrip::String)
# function getCorrectPatternsMarkersABC(patSequence::Vector{Char}, corrTrip::String)
#     correctPatternsMarkers = Int[]
#     word = collect(corrTrip)
#     sizePS = length(patSequence)
#     i = 0
#     word_size = 30
#     where_to_look_at = patSequence[1:word_size:end]

#     while i <= sizePS
#         i += 1
#         j = 1

#         while i <= sizePS && patSequence[i] == word[j]#  && j <= length(word)

#             if j == length(word)
#                 push!(correctPatternsMarkers, i)
#             end

#             # Skip past current pattern
#             i += word_size
#             # while i <= sizePS && patSequence[i] != 'Z'
#             #     i += 1
#             # end
#             # while i <= sizePS && patSequence[i] == 'Z'
#             #     i += 1
#             # end
#             j += 1
#             if j > length(word)
#                 break
#             end
#         end

#         # Skip partial matches
#         # while i <= sizePS && (patSequence[i] == word[2] || patSequence[i] == word[3])
#         #     i += 1
#         # end
#         # while i <= sizePS && patSequence[i] == 'Z'
#         #     i += 1
#         # end

#         if i >= sizePS
#             break
#         end
#         i -= 1  # Adjust index for next iteration
#     end

#     return correctPatternsMarkers
# end

# function getCorrectPatternsMarkersABCDEF(patSequence::Vector{Char}, corrTrip::String)
#     correctPatternsMarkers = Int[]
#     word = corrTrip
#     sizePS = length(patSequence)
#     for i in 1:sizePS
#         j = 1
#         while i <= sizePS && j <= length(word) && patSequence[i] == word[j]
#             j += 1
#             if j > length(word)
#                 push!(correctPatternsMarkers, i)
#                 while i <= sizePS && patSequence[i] != 'Z'
#                     i += 1
#                 end
#                 while i <= sizePS && patSequence[i] == 'Z'
#                     i += 1
#                 end
#             end
#         end
#     end
#     return correctPatternsMarkers
# end



# function getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String, endIndex::Int)
#     correctPatternsMarkers = Int[]
#     word = corrTrip
#     sizePS = length(patSequence)
#     i::Int = 1
#     did_increment = false

#     while i <= sizePS
#         j = 1
#         while i <= sizePS && patSequence[i] == word[j] && j <= length(word)
#             j += 1
#             if j > length(word)
#                 push!(correctPatternsMarkers, i)
#             end
#             while i <= sizePS && patSequence[i] != 'Z'
#                 i += 1
#             end
#             while i <= sizePS && patSequence[i] == 'Z'
#                 i += 1
#             end
#         end

#         while i <= sizePS && occursin(patSequence[i], word[2:endIndex])
#             i += 1
#             did_increment = true
#         end

#         while i <= sizePS && patSequence[i] == 'Z'
#             i += 1
#             did_increment = true
#         end

#         if did_increment
#             i -= 1
#             did_increment = false
#         end
#     end

#     return correctPatternsMarkers
# end

# getCorrectPatternsMarkersABCDEFGHIJ(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 10)

# getCorrectPatternsMarkersABCDEFGHI(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 9)

# getCorrectPatternsMarkersABCDEFGH(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 8)

# getCorrectPatternsMarkersABCDEFG(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 7)

# getCorrectPatternsMarkersABCDEF(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 6)

# getCorrectPatternsMarkersABCDE(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 5)

# getCorrectPatternsMarkersABCD(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 4)

# getCorrectPatternsMarkersABC(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 3)

# getCorrectPatternsMarkersAB(patSequence::Vector{Char}, corrTrip::String) =
#     getCorrectPatternsMarkers(patSequence, corrTrip, 2)
@deprecate getCorrectPatternsMarkersABCDEFGHI(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABCDEFGH(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABCDEFG(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABCDEF(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABCDE(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABCD(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
@deprecate getCorrectPatternsMarkersABC(patSequence::Vector{Char}, corrTrip::String) getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)

function getCorrectPatternsMarkers(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = collect(corrTrip)
    sizePS = length(patSequence)
    i = 1

    while i <= sizePS
        j = 1

        # Try to match the full word
        while i <= sizePS && j <= length(word) && patSequence[i] == word[j]
            j += 1
            i += 1
            if j > length(word)
                push!(correctPatternsMarkers, i - 1)  # Match found, record last index
            end

            # Skip non-'Z' characters
            while i <= sizePS && patSequence[i] != 'Z'
                i += 1
            end
            # Skip all 'Z' characters
            while i <= sizePS && patSequence[i] == 'Z'
                i += 1
            end
        end

        # Skip partial matches (letters from word[2] to end)
        did_increment = false
        while i <= sizePS && patSequence[i] in word[2:end]
            i += 1
            # did_increment = true
        end
        while i <= sizePS && patSequence[i] == 'Z'
            i += 1
            # did_increment = true
        end
        if did_increment
            i -= 1  # Adjust index for next iteration
        end
    end

    return correctPatternsMarkers
end



function getCorrectPatternsMarkers_my_adaptation(patSequence::Vector{Char}, corrTrip::String)
    correctPatternsMarkers = Int[]
    word = corrTrip
    letters_in_word = (combinations(word, length(word))|>collect)[1]
    sizePS = length(patSequence)
    i = 0
    letterSize = 6 # TODO set this as an argument for the function
    silenctInterval = 24# TODO set this as an argument for the function
    word_size = letterSize + silenctInterval
    where_to_look_at = patSequence[1:word_size:end]

    while i <= sizePS
        i += 1
        j = 1
        while i <= sizePS && patSequence[i] == word[j]#  && j <= length(word)
            if j == length(word)
                push!(correctPatternsMarkers, i)
            end
            # while i <= sizePS && patSequence[i] != 'Z'
            #     i += 1
            # end
            # while i <= sizePS && patSequence[i] == 'Z'
            #     i += 1
            # end
            i += word_size
            j += 1
            if j > length(word)
                break
            end
        end
        # Skip partial matches
        if i >= sizePS
            break
        end

        comparison = patSequence[i] .== letters_in_word[2:end]
        if i <= sizePS && any(comparison)
            # while i <= sizePS && (patSequence[i] == word[2] || patSequence[i] == word[3] || patSequence[i] == word[4])
            comparison = patSequence[i] .== letters_in_word
            i += letterSize
        end

        if i <= sizePS && patSequence[i] == 'Z'
            i += silenctInterval
        end

        if i >= sizePS
            break
        end
        i -= 1  # Adjust index for next iteration
    end
    return correctPatternsMarkers
end


function generatePermutationsWithReplacement(str::String, current::String, permutation_length::Int, permutations::Vector{String}, prefix::String, suffix::String)
    if length(current) == length
        push!(permutations, prefix * current * suffix)
        return
    end
    for c in str
        generatePermutationsWithReplacement(str, current * c, permutation_length, permutations, prefix, suffix)
    end
end

function getRandomSequenceGivenCorrSignal(chSequence::String, size::Int, gap::Int, letterSize::Int, sequence::String, variationOnSignal, variationOnSilence)
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    # insertSequenceIntoLetterChain!(sequence, randSequence, 100)
    new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, total_insertions)
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function getPermutations(sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String})
    total_sequences = length(sequence_heads)
    all_permutations = Vector{Vector{String}}(undef, total_sequences)

    for i in 1:total_sequences
        # some_vector = String[]
        # possible_permutations =
        # generatePermutationsWithReplacement(permutation_bases[i], "", 3, some_vector, sequence_heads[i], sequence_tails[i])
        chars = (combinations(permutation_bases[i], 3)|>collect)[1]
        all_permutations[i] = # some_vector
            sequence_heads[i] .* ([Iterators.product(fill(chars, 3)...)...] .|> y -> join(y, "")) .* sequence_tails[i]
    end
    return all_permutations
end

function generateSequenceWithPermutation(num_of_letters::Int, letterSize::Int, sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String})
    total_sequences = length(sequence_heads)
    all_permutations = getPermutations(sequence_heads, permutation_bases, sequence_tails)

    totalSize = sum(length.(all_permutations))
    # Generation of chSequence could be replaced with a vcat of all_permutations
    chSequence = fill("", totalSize)

    index = 1
    for perms in all_permutations
        for s in perms
            chSequence[index] = s
            index += 1
        end
    end

    # TODO It might be better to have randSequence as a vector of char, rather a string- a string is a constructino as whole, while this is not needed
    single_seq_len = length(chSequence[1])
    how_many_sequences_to_sample = max(ceil(Int, num_of_letters / single_seq_len), 1) + 1
    @info "single_seq_len = $(single_seq_len)"
    @info "num_of_letters = $(num_of_letters)"
    @info "how_many_sequences_to_sample = $(how_many_sequences_to_sample)"

    randSequence = ""
    for _ in 1:how_many_sequences_to_sample
        randSequence *= rand(chSequence)
    end

    @assert length(randSequence) >= num_of_letters
    return randSequence[1:num_of_letters]
end

function getRandomSequenceGivenCorrSignalUpgraded(chSequence::String, size::Int, gap::Int, letterSize::Int, sequence::String, variationOnSignal, variationOnSilence)
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end

    total_insertions = 10
    new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, total_insertions)

    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function get_abcXXX_XXXdef_Sequence(correctTriplet::String, size::Int, letterSize::Int,)
    # srand(time(0))
    # sequence = "ABCDEF"
    sequence_heads = ["ABC", "", "A", "AB",]
    permutation_bases = ["DEF", "ABC", "BCD", "CDE",]
    sequence_tails = ["", "DEF", "EF", "F"]
    # return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)
    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )

    return randSequence
end

function get_abcdXXX_XXXdefg_Sequence(correctTriplet::String, size::Int, letterSize::Int,)# variationOnSignal, variationOnSilence)
    # srand(time(0))
    sequence = "ABCDEFG"
    sequence_heads = ["ABCD", "", "AB", "ABC", ""]
    permutation_bases = ["EFG", "ABC", "CDE", "DEF", "ABC",]
    sequence_tails = ["", "DEFG", "FG", "G", "EDEFG",]
    # return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)
    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )
    return randSequence

    # insertionWindowSize = 10
    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)

    # expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)
    # return expanded_sequence, all_insertions
end

function getABXXX_XXXDE_Sequence(correctTriplet::String, size::Int, letterSize::Int,)# variationOnSignal, variationOnSilence)
    sequence = "ABCDE"
    sequence_heads = ["AB", ""]
    permutation_bases = ["CDE", "ABC"]
    sequence_tails = ["", "DE"]
    # return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)

    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )

    return randSequence

    # insertionWindowSize = 10
    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)

    # expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)
    # return expanded_sequence, all_insertions
end

function getABCDEFSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    # srand(time(0))
    randSequence = ""
    chSequence = ["ABCAAF", "ABCABF", "ABCACF", "ABCADF", "ABCAEF", "ABCBAF", "ABCBBF", "ABCBCF", "ABCBDF", "ABCBEF", "ABCCAF", "ABCCBF", "ABCCCF", "ABCCDF", "ABCCEF", "ABCDAF", "ABCDBF", "ABCDCF", "ABCDDF", "ABCDEF", "ABCEAF", "ABCEBF", "ABCECF", "ABCEDF", "ABCEEF"]
    arraySize = length(chSequence)

    for _ in 1:size
        chosenChar = chSequence[rand(1:arraySize)]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function getABCDSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence; do_initial_version::Bool=false)
    if do_initial_version
        chSequence = ["ABCD", "ABCD", "ABCD", "ABCD", "ABCD", "BACD", "ACCD", "CACD",
            "ADCD", "DACD", "BCCD", "CBCD", "CDCD", "DCCD", "BDCD", "DBCD", "AACD",
            "BBCD", "CCCD", "DDCD", #// 20
            "ABAB", "ABBA", "ABAC", "ABCA", "ABAD", "ABDA", "ABBC", "ABCB",
            "ABCD", "ABDC", "ABBD", "ABDB", "ABAA", "ABBB", "ABCC", "ABDD",
            "BBCD", "BBCD", "CBCD", "DBCD", #// 20
            "ADABCD", "ADABCD", "ADABCD", "ABABCD", "ABABCD", "ABABCD",# // 6
            "ABBCD", "ABCCD", "ABCDD", "ABBCD", "ABCCD", "ABCDD", #// 6
            "ABDCD", "ADBCD", "ACBCD", "ABDCD", "ADBCD", "ACBCD", #// 6
            "ABDDCD", "ACDBCD", "ABBCD", "ACDBCD", "ACBBCD", #/*5*/
            "ACDBCD", "ABDD", "ADCBCD", "ABDDD", "ABDDDD", #/*5*/
            "ABBCCD", "ABCCD", "ABCCCD", "ABCCD", "ABCCD", "ABCCCD", # 6
            "ACCBCD", "ACCBCD", "ABBCD", "ABCCBCCD", "ACBCD", "ABCBCD" # 6
        ] # 80
        arraySize = length(chSequence)
        randSequence = ""
        for _ in 1:size
            chosenChar = chSequence[rand(1:arraySize)]
            randSequence *= chosenChar
        end
        return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
    else
        sequence::String = "ABCD"
        sequence_heads = ["", "A",]
        permutation_bases = ["ABC", "BCD"]
        sequence_tails = ["D", ""]

        randSequence = generateSequenceWithPermutation(
            size,
            letterSize,
            sequence_heads,
            permutation_bases,
            sequence_tails,
        )

        insertionWindowSize = 10
        new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)
        # positionInSubSequence = 3
        # new_sequence =
        #     new_sequence[1:(positionInSubSequence-1)] * sequence * new_sequence[(positionInSubSequence+length(sequence)):end]

        # # For debuggin- remove once done >>>
        # pattern = "AABABCDACC"
        # new_sequence = pattern * new_sequence[11:end]
        # # for i in 1:10
        # #     new_sequence[i] = pattern[i]
        # # end
        # # For debuggin- remove once done <<<

        expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)
        return expanded_sequence, all_insertions
    end
end


function getABCSequence(correctTriplet::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence; do_initial_version::Bool=false)
    # srand(time(0))
    randSequence = ""
    if do_initial_version
        chSequence = ["AAA", "AAB", "AAC", "ABA", "ABB", "ABC", "ACA", "ACB", "ACC",
            "BAA", "BAB", "BAC", "BBA", "BBB", "BBC", "BCA", "BCB", "BCC",
            "CAA", "CAB", "CAC", "CBA", "CBB", "CBC", "CCA", "CCB", "CCC",
            "ABC", "ABC", "ABC", "ABC", "ABC", "ABC", "ABBC", "ABBBC", "ABBBBC",
            "AAC", "AAAC", "AABAC", "ABAAAAC", "BBABBC", "BCABBC", "CCABBC", "CBABBC", "ABABBC",
            "ACABBC", "BAABBC", "AAABBC", "CAABBC", "ACCBBC"]
        arraySize = length(chSequence)

        for _ in 1:size
            chosenChar = chSequence[rand(1:arraySize)]
            randSequence *= chosenChar
        end
        return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
    else
        sequence::String = "ABC"
        sequence_heads = ["",]
        permutation_bases = ["ABC",]
        sequence_tails = ["",]

        randSequence = generateSequenceWithPermutation(
            size,
            letterSize,
            sequence_heads,
            permutation_bases,
            sequence_tails,
        )

        insertionWindowSize = 10
        new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)
        positionInSubSequence = 3
        new_sequence =
            new_sequence[1:(positionInSubSequence-1)] * sequence * new_sequence[(positionInSubSequence+length(sequence)):end]
        expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)

        return expanded_sequence, all_insertions
    end

end


function getUserDefinedSequence(inputStr::String)
    # TODO Why is this function needed?
    randSequence = inputStr
    return randSequence
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

function getPermutedSequence(fileName::String, str::String, k::Int, gap::Int, lettersize::Int, variationOnSignal, variationOnSilence)
    open(fileName, "r") do file
        allPerm = read(file, String)
        return insertGapsAndSetLetterSize(allPerm, gap, lettersize, variationOnSignal, variationOnSilence)
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

function getAllPossibleCombinationSequence(gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    allPatterns = ""
    signals = "ABCDEF"
    pos = zeros(Int, length(signals))
    allPossibleCombinations(allPatterns, signals, pos, 0)
    return insertGapsAndSetLetterSize(allPatterns, gap, letterSize, variationOnSignal, variationOnSilence)
end

function get_abcdefgXXX_XXXdefghij_Sequence(correctTriplet::String, size::Int, letterSize::Int)#, variationOnSignal, variationOnSilence)
    # srand(time(0))
    sequence = "ABCDEFGHIJ"
    sequence_heads = ["ABCDEFG", "", "ABC", "ABCD"]
    permutation_bases = ["HIJ", "ABC", "DEF", "EFG"]
    sequence_tails = ["", "DEFGHIJ", "GHIJ", "HIJ"]
    # return generateSequenceWithPermutation_with_doubling(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)
    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )
    return randSequence

    # insertionWindowSize = 50
    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)
    # positionInSubSequence = 3
    # new_sequence =
    #     new_sequence[1:(positionInSubSequence-1)] * sequence * new_sequence[(positionInSubSequence+length(sequence)):end]
    # expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)

    # return expanded_sequence, all_insertions
end

function get_abcdefXXX_XXXdefghi_Sequence(correctTriplet::String, size::Int, letterSize::Int)#, variationOnSignal, variationOnSilence)
    # srand(time(0))
    sequence = "ABCDEFGHI"
    sequence_heads = ["ABCDEF", "", "ABC"]
    permutation_bases = ["GHI", "ABC", "DEF"]
    sequence_tails = ["", "DEFGHI", "GHI"]
    # return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)
    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )

    # insertionWindowSize = 50
    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)
    # positionInSubSequence = 3
    # new_sequence =
    #     new_sequence[1:(positionInSubSequence-1)] * sequence * new_sequence[(positionInSubSequence+length(sequence)):end]
    # expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)

    # return expanded_sequence, all_insertions
end


function get_abcdeXXX_XXXfgh_Sequence(correctTriplet::String, size::Int, letterSize::Int)#, variationOnSignal, variationOnSilence;)
    # srand(time(0))
    sequence = "ABCDEFG"
    sequence_heads = ["ABCDE", "", "ABC"]
    permutation_bases = ["FGH", "ABC", "DEF"]
    sequence_tails = ["", "DEFGH", "GH"]

    # return generateSequenceWithPermutation(size, gap, letterSize, sequence_heads, permutation_bases, sequence_tails, sequence, variationOnSignal, variationOnSilence)
    randSequence = generateSequenceWithPermutation(
        size,
        letterSize,
        sequence_heads,
        permutation_bases,
        sequence_tails,
    )
    return randSequence

    # insertionWindowSize = 50
    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, insertionWindowSize)
    # positionInSubSequence = 3
    # new_sequence =
    #     new_sequence[1:(positionInSubSequence-1)] * sequence * new_sequence[(positionInSubSequence+length(sequence)):end]
    # expanded_sequence = insertGapsAndSetLetterSize(new_sequence, gap, letterSize, variationOnSignal, variationOnSilence)

    # return expanded_sequence, all_insertions
end


function generateSequenceWithPermutation_with_doubling(size::Int, gap::Int, letterSize::Int, sequence_heads::Vector{String}, permutation_bases::Vector{String}, sequence_tails::Vector{String}, sequence::String, variationOnSignal, variationOnSilence)
    "Not implemented" |> ErrorException |> throw
    # total_sequences = length(sequence_heads)
    # all_permutations = getPermutations(sequence_heads, permutation_bases, sequence_tails, variationOnSignal)

    # totalSize = sum(length.(all_permutations))
    # chSequence = String[]

    # index = 1
    # for perms in all_permutations
    #     for s in perms
    #         chSequence[index] = s
    #         index += 1
    #     end
    # end

    # randSequence = ""
    # for _ in 1:size
    #     chosenSeq = chSequence[rand(1:totalSize)]
    #     randSequence *= chosenSeq
    # end

    # new_sequence, all_insertions = insertSequenceIntoLetterChain(sequence, randSequence, 100)
    # return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function getPermutationWithReplacement(sequence::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    totalPermutations = Vector{String}()
    generatePermutationsWithReplacement(sequence, "", size, totalPermutations, "", "")
    randSequence = totalPermutations[rand(1:length(totalPermutations))]
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function getRandomCharacterSequence(size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    chSequence = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    randSequence = ""
    for _ in 1:size
        chosenChar = chSequence[rand(1:length(chSequence))]
        randSequence *= chosenChar
    end
    return insertGapsAndSetLetterSize(randSequence, gap, letterSize, variationOnSignal, variationOnSilence)
end

function insertSequenceIntoLetterChain(sequence::String, randSequence::String, insertionWindowSize::Int)

    sequence_elements = collect(sequence)
    net_input_elements = collect(randSequence)
    seqLength = length(sequence)
    # positionInSubSequence = 0
    all_insertions = Int[]

    length_before_insertion = length(net_input_elements)
    # TODO it might be more efficient to convert this to a vector, then after its done convert into a string
    for i in 1:insertionWindowSize:(length_before_insertion-insertionWindowSize)
        positionInSubSequence = rand(1:(insertionWindowSize-seqLength-1))
        insertion_start = i + positionInSubSequence
        insertion_end = i + positionInSubSequence + seqLength - 1
        # a = length(insertion_start:insertion_end)
        where_to_put = insertion_start:insertion_end
        net_input_elements[where_to_put] .= sequence_elements
        # new_sequence[1:insertion_start] * sequence * new_sequence[insertion_end:end]
        push!(all_insertions, insertion_end - 1)
    end
    @assert length_before_insertion == length(net_input_elements)

    new_sequence = String(net_input_elements)
    return new_sequence, all_insertions
end

function insertSequenceIntoLetterChain(phrase::String, phrases_with_dobuling_vec::Vector{String}, randSequence::String, insertionWindowSize::Int)
    half_insertionWindowSize = div(insertionWindowSize, 2)
    phrase_len = length(phrase)
    phrase_with_dobuling = phrases_with_dobuling_vec[1]
    phrase_with_doubling_len = length(phrase_with_dobuling)
    # total_phrases = length(phrases_with_dobuling_vec)

    # for i in 1:sub(length(randSequence), insertionWindowSize)
    for i in 1:(length(randSequence)-insertionWindowSize)
        positionInSubSequence = rand(1:(half_insertionWindowSize-phrase_len-1))
        position_in_sub_sequence_for_doubled_letter = rand(1:(half_insertionWindowSize-phrase_with_doubling_len-2)) + half_insertionWindowSize

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
# function insertSequenceIntoLetterChain(sequence::String, placeholders::Vector{String}, randSequence::String, insertionWindowSize::Int)
#     seqLength = length(placeholders)
#
#     for i in 1:(length(randSequence)-insertionWindowSize)
#         position = rand(1:(length(randSequence)-seqLength))
#         for j in 1:seqLength
#             randSequence[position+j-1] = placeholders[j]
#         end
#     end
# end


function generateSequenceWithReplacement(sequence::String, size::Int, gap::Int, letterSize::Int, variationOnSignal, variationOnSilence)
    generatedSequence = ""
    for _ in 1:size
        generatedSequence *= sequence[rand(1:length(sequence))]
    end
    return insertGapsAndSetLetterSize(generatedSequence, gap, letterSize, variationOnSignal, variationOnSilence)
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


function get_signal_of_len(signal_len::Int)
    return 'A':'Z' |> y -> y[1:signal_len] |> collect |> z -> join(z, "")
end

end #module