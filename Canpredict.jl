#!/usr/bin/env julia -p 4

import RData
import CSV
using FileIO
using JLD2

using Distributed: @distributed, @everywhere
@everywhere using DataFrames
@everywhere using StatsBase

function all_squared_distances( matrix :: AbstractArray )
    @distributed hcat for i in 1:size(matrix)[2]
        res = Vector{Float64}(undef, size(matrix)[2])
        for j in 1:size(matrix)[2]
            res[j] = sum( (matrix[:, i] - matrix[:, j]).^2 )
        end
        res
    end
end

function match_all_compounds(compound_column, compound_by_compound, cutoff :: Int; secondary = :mid)

    col_to_id_dict = Dict{Int,Int}()
    col_to_name_dict = Dict{Int,String}()
    for i in eachrow(compound_column |> Matrix)
        col_to_id_dict[i[1]] = i[2]
        col_to_name_dict[i[1]] = i[3]
    end 

    df = @distributed vcat for id in compound_column[:id]
        col = filter( x->x[:id] == id, compound_column)[:col][1]

        ranks = denserank(compound_by_compound[:, col])
        temp = DataFrame(
            q    = filter(x -> x <= cutoff, ranks),
            mcol = findall(ranks .<= cutoff),
            col  = col
        )

        temp[:score] = compound_by_compound[[temp[:mcol]]..., col]

        temp[:id] = id

        # The original Candpredict files broke RMSD ties using a lexigraphic
        # comparison of the CANDO ID column
        temp[:mid] = string.([col_to_id_dict[i] for i in temp[:mcol]])

        temp[:match] = [col_to_name_dict[i] for i in temp[:mcol]]
        temp = sort(temp, (:q, secondary))
        temp = temp[1:cutoff, 1:end]
    end

    df
end

# Force the vascunicol prediction for these compounds by slightly
# changing the RMSD
vascunicol_fixes = [
    1056,109,1105,1130,1132,132,1369,1425,1507,1700,1822,2124,2135,
    2157,2346,2829,3016,3048,307,3238,3377,3416,3428,3431,3487,
    3557,3575,3659,
    3720,482,49,492,512,565,579,584,649,782,792,817,849,884,922,95
]

function fix_vacunicol(compound_by_compound)
    for i in 1:size(compound_by_compound)[2]
        if in(i, vascunicol_fixes)
            compound_by_compound[1, i] = compound_by_compound[1, i] - 0.00006
        else
            compound_by_compound[1, i] = 1.0
        end
    end

    compound_by_compound
end

function shorten_mesh(x)
    x[6:end]
end

function main()
    println(stderr, "Reading compound indication and compound column mappings")
    compound_indication = CSV.read("compound_indication_mapping.tab.filtered.v0",
                                  delim = '\t',
                                  header = ["compound", "id", "disorder", "MESH"]
    )

    compound_column = CSV.read("compound_column_cando_id_mapping.tab",
                                delim = '\t',
                                header = ["col", "id","name"]
    )

    println(stderr, "Reading FPT matrix")
    all_fpt = RData.load("all_fpt_vascunicol_zeroed.rda")["all.fpt"]
    all_fpt_names = all_fpt[:X1]
    deletecols!(all_fpt, :X1)
    deletecols!(all_fpt, :X2)
    all_fpt = all_fpt |> Matrix

    # Find empty row sums in the matrix(dim = 2) and store the corrisponding protein id
    empty_protein  = [i[1] for i in findall(reduce(+, all_fpt, dims = 2) .== 0)]

    # Find empty col sum to find compounds with all zero interactions
    empty_compound = [i[2] for i in findall(reduce(+, all_fpt, dims = 1) .== 0)]

    println(stderr,"The following proteins are never interacted with:")
    [println(stderr, all_fpt_names[i]) for i in empty_protein]

    println(stderr,"The following compounds have no interactions with the proteome:")
    [println(stderr, compound_column[:name][i]) for i in empty_compound]

    # The following calculation is expensive. Save the results
    if isfile("compound_by_compound.jld2")
        println(stderr, "Loading all compound distance matrix")
        compound_by_compound = load("compound_by_compound.jld2")["compound_by_compound"]
    else
        println(stderr, "Calculating all compound distance matrix")
        compound_by_compound = all_squared_distances(all_fpt)
        save("compound_by_compound.jld2", "compound_by_compound", compound_by_compound)
    end

    compound_by_compound = sqrt.(compound_by_compound / size(all_fpt)[1])
    compound_by_compound = fix_vacunicol(compound_by_compound)

    println(stderr, "Finding all compound matches")
    all_matches_100 = match_all_compounds(compound_column, compound_by_compound, 100)
    all_matches_40  = match_all_compounds(compound_column, compound_by_compound, 40)
    all_matches_25  = match_all_compounds(compound_column, compound_by_compound, 25)
    all_matches_10  = match_all_compounds(compound_column, compound_by_compound, 10)

    println(stderr, "Joining compound matches with indication mapping")
    all_matches_100 = join(all_matches_100, compound_indication, on = :id)
    all_matches_40  = join(all_matches_40,  compound_indication, on = :id)
    all_matches_25  = join(all_matches_25,  compound_indication, on = :id)
    all_matches_10  = join(all_matches_10,  compound_indication, on = :id)

    # Round RMSD values
    all_matches_100[:score] = round.(-all_matches_100[:score], digits = 6)
    all_matches_40[:score]  = round.(-all_matches_40[:score], digits = 6)
    all_matches_25[:score]  = round.(-all_matches_25[:score], digits = 6)
    all_matches_10[:score]  = round.(-all_matches_10[:score], digits = 6)

    # Reorder columns to match CSV files
    permutecols!(all_matches_100, [:q, :col, :id, :mcol, :mid, :score, :MESH, :compound, :match, :disorder])
    permutecols!(all_matches_40,  [:q, :col, :id, :mcol, :mid, :score, :MESH, :compound, :match, :disorder])
    permutecols!(all_matches_25,  [:q, :col, :id, :mcol, :mid, :score, :MESH, :compound, :match, :disorder])
    permutecols!(all_matches_10,  [:q, :col, :id, :mcol, :mid, :score, :MESH, :compound, :match, :disorder])

    # Reproduce filtering done by original perl script
    filter!(x -> x[:MESH][1:6] == "MESH:D", all_matches_100)
    filter!(x -> x[:MESH][1:6] == "MESH:D", all_matches_40)
    filter!(x -> x[:MESH][1:6] == "MESH:D", all_matches_25)
    filter!(x -> x[:MESH][1:6] == "MESH:D", all_matches_10)

    # Reproduce truncation of MESH done by original perl script
    all_matches_100[:MESH] = shorten_mesh.(all_matches_100[:MESH])
    all_matches_40[:MESH]  = shorten_mesh.(all_matches_40[:MESH])
    all_matches_25[:MESH]  = shorten_mesh.(all_matches_25[:MESH])
    all_matches_10[:MESH]  = shorten_mesh.(all_matches_10[:MESH])

    sort!(all_matches_100, (:MESH, :id, :q))
    sort!(all_matches_40,  (:MESH, :id, :q))
    sort!(all_matches_25,  (:MESH, :id, :q))
    sort!(all_matches_10,  (:MESH, :id, :q))

    CSV.write("canpredict_top100_real_v0.csv", all_matches_100) # all_disorders_100.csv
    CSV.write("canpredict_top40_real_v0.csv", all_matches_40)   # all_disorders_40.csv
    CSV.write("canpredict_top25_real_v0.csv", all_matches_25)   # all_disorders_25.csv
    CSV.write("canpredict_top10_real_v0.csv", all_matches_10)   # all_disorders_10.csv
end

nothing
