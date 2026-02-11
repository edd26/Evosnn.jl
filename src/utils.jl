import DrWatson: srcdir, @quickactivate
@quickactivate "Evosnn"

import DelimitedFiles
# include(srcdir("Individual.jl"))
# include(srcdir("Parameters.jl"))

# TODO verify correctness with C++
function printIndividualMatrix(individual::Individual, params::Parameters, gNo::Int, iRun::Int, saving_dir::String, orgi_file_name::String;
    file_core_name::String="ind_Matrix_s", file_extension::String=".txt", column_separator::Char='\t')
    file_parameters = string(params.noOfSignals)
    ispath(saving_dir) || mkpath(saving_dir)
    final_path = joinpath(saving_dir, file_core_name * file_parameters * "_" * orgi_file_name)

    if !occursin( file_extension, final_path)
      final_path *= file_extension
    end

    @info "Saving dir: $(saving_dir)"
    @info "Final path,: $(final_path)"

    open(final_path, "a") do ofs
        println(ofs, "Adj. matrix of the network")

        if column_separator == '\t'
            for i in 1:individual.noOfNodesInNetwork
                for j in 1:individual.noOfNodesInNetwork
                    print(ofs, individual.indMatrix[i, j], column_separator)
                end
                println(ofs)
            end
        elseif column_separator == ' '
            max_size = 0
            for i in 1:individual.noOfNodesInNetwork
                for j in 1:individual.noOfNodesInNetwork
                    str_element = string(individual.indMatrix[i, j])
                    max_size = max(max_size, length(str_element))
                end
            end

            separator = column_separator^max_size

            for i in 1:individual.noOfNodesInNetwork
                for j in 1:individual.noOfNodesInNetwork
                    if individual.indMatrix[i, j] == 0
                        local_separator = separator
                    else
                        str_element = string(individual.indMatrix[i, j])
                        new_size = max_size - length(str_element)
                        if individual.indMatrix[i, j] < 0
                            new_size += 1
                        end
                        local_separator = join(fill(" ", (new_size + 2)))
                    end
                    print(ofs, individual.indMatrix[i, j], local_separator)
                end
                println(ofs)
            end
        else
            "Did not recognise the separator of columns." |> ErrorException |> throw
        end

        if gNo == 0 && iRun == 0
            println(ofs, "\nPerformance of the nw on a random sequence of length ", params.reevaluateSeq, " signals")
        end

        # println(ofs, "run:", iRun, "\tgen:", gNo, "\tfitness = ", individual.fitness, "\tTotal corr = ", individual.totalCorrPatterns, "\tidentified corr = ", individual.reward, "\tResponded to incorrect patterns = ", individual.penalty, "\tFDR = ", individual.fdr, "\tPrecision = ", individual.precision)
        println(ofs, "run:", iRun, "\tgen:", gNo)
        println(ofs, "\tfitness = ", individual.fitness,)
        println(ofs, "\tTotal corr = ", individual.totalCorrPatterns)
        println(ofs, "\tidentified corr = ", individual.reward)
        println(ofs, "\tResponded to incorrect patterns = ", individual.penalty)
        println(ofs, "")
        println(ofs, "Total correct patterns in the sequence:", individual.totalCorrPatterns)

        println(ofs, "\nConfusion matrix:")
        TP = individual.reward
        # TODO ask: why there is here reevaluateSeq and not noOfLetters?
        # TN = (params.reevaluateSeq - individual.totalCorrPatterns) - individual.penalty
        TN = (params.noOfLetters - individual.totalCorrPatterns) - individual.penalty
        FP = individual.penalty
        FN = individual.totalCorrPatterns - individual.reward

        println(ofs, "True positives:", TP)
        println(ofs, "True negatives:", TN)
        println(ofs, "False positives:", FP)
        println(ofs, "False negatives:", FN)
        println(ofs, "Sum of all:", TP + TN + FP + FN)
        println(ofs, "_____________________________________________________________________________\n")
    end
end

function export_matrix_to_tsv(individual::Individual, params::Parameters, saving_dir::String, orgi_file_name::String;
    file_core_name::String="matrix_export", file_extension::String=".txt", name_suffix::String="")

    file_parameters = string(params.noOfSignals)
    ispath(saving_dir) || mkpath(saving_dir)

    local_file_name  = replace(orgi_file_name, ".txt"=>"")

    # final_path = joinpath(saving_dir,local_file_name * "_"* file_core_name * file_parameters *  file_extension)
    final_path = joinpath(saving_dir,local_file_name * file_extension)
    @info "Export file as 'tsv': $(final_path)"

    export_mat = individual.indMatrix

    open(final_path, "w") do io
        DelimitedFiles.writedlm(io, export_mat, '\t')
    end
end


# TODO verify correctness with C++
function printIndividualMatrix(individual::Individual, params::Parameters, gNo::Int, iRun::Int;
    saving_dir::String="outputs3/", file_core_name="ind_Matrix_s", file_extension=".txt")

    file_parameters = string(params.noOfSignals)

    final_path = joinpath(saving_dir, file_core_name * file_parameters)

    if !occursin(final_path, file_extension)
      final_path *= file_extension
    end

    open(final_path, "a") do ofs
        println(ofs, "Adj. matrix of the network")

        max_size = 0
        for i in 1:individual.noOfNodesInNetwork
            for j in 1:individual.noOfNodesInNetwork
                str_element = string(individual.indMatrix[i, j])
                max_size = max(max_size, length(str_element))
            end
        end

        separator = "\t"

        for i in 1:individual.noOfNodesInNetwork
            for j in 1:individual.noOfNodesInNetwork
                print(ofs, individual.indMatrix[i, j], separator)
            end
            println(ofs)
        end

        if gNo == 0 && iRun == 0
            println(ofs, "\nPerformance of the nw on a random sequence of length 10 thousand signals")
        end

        println(ofs, "run:", iRun, "\tgen:", gNo, "\tfitness = ", individual.fitness, "\tTotal corr = ", individual.totalCorrPatterns, "\tidentified corr = ", individual.reward, "\twrong = ", individual.penalty, "\tFDR = ", individual.fdr, "\tPrecision = ", individual.precision)
        println(ofs, "\nIndividual parameters:")
        println(ofs, "fitness:", individual.fitness)
        println(ofs, "reward:", individual.reward)
        println(ofs, "penalty:", individual.penalty)
        println(ofs, "rewardn (a.k.a. identified corr):", individual.rewardn)
        println(ofs, "penalty (a.k.a. wrong):", individual.penaltyn)
        println(ofs, "_____________________________________________________________________________\n")
    end
end


function removeLowWeights(ind::Individual, params::Parameters)
    where_to_keep = params.minConnectionWeight .< ind.indMatrix .&& ind.indMatrix .< params.maxConnectionWeight
    new_weights = copy(ind.indMatrix)
    new_weights[.!where_to_keep] .= 0
    return new_weights
end
