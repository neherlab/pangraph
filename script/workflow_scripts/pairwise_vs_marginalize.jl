##
using PanGraph
using PanGraph.Graphs
using JSON


struct Interval # closed sequence interval
    l::Int # left edge
    r::Int # right edge
    L::Int # genome length
    bs::Vector{String} # block name(s)
end

# length of the interval, keeping into account periodic boundary conditions
Base.length(I::Interval)::Int = I.l ≤ I.r ? I.r - I.l + 1 : (I.r + I.L) - I.l + 1

# given a set of left segment edges, and the total genome length,
#  finds the set of corresponding right edges
Redge_from_Ledge(begs::Vector{Int}, L::Int)::Vector{Int} =
    circshift(mod.(begs .- 2, L) .+ 1, -1)

# given a graph and an isolate, returns the interval representation of the isolate
function to_intervals(G::Graph, isolate::String)::Vector{Interval}
    L = length(sequence(G, isolate))
    seq = G.sequence[isolate]
    begs = seq.position[1:end-1]
    ends = Redge_from_Ledge(begs, L)
    @assert length(begs) == length(seq.node)
    Is = [Interval(b, e, L, [n.block.uuid]) for (n, b, e) in zip(seq.node, begs, ends)]
    @assert length(Is) == length(seq)
    @assert all(length.(Is) .> 0)
    return Is
end

# merges two interval partitions of the same isolate, retaining names of blocks in both partitions
function merge(I1::Vector{Interval}, I2::Vector{Interval})::Vector{Interval}
    L = I1[1].L
    @assert all([I.L for I in I1] .== L)
    @assert all([I.L for I in I2] .== L)

    Idict = Dict()
    Pos = Dict()
    for s ∈ [:r, :l], (n, I) ∈ enumerate([I1, I2])
        # dictionary [position -> interval] for each pair (r/l, 1/2)
        Idict[(s, n)] = Dict(getfield(i, s) => i for i in I)
        Pos[(s, n)] = sort(collect(keys(Idict[(s, n)])))
    end

    # find the block corresponding to one particular position
    function get_blocks(n::Int, pos::Int)::Vector{String}
        key = (:l, n)
        Ps = Pos[key]
        w = searchsorted(Ps, pos).stop
        p = (w == 0) ? Ps[end] : Ps[w]
        return Idict[key][p].bs
    end

    new_lPs = sort(unique(vcat(Pos[(:l, 1)], Pos[(:l, 2)])))
    new_rPs = Redge_from_Ledge(new_lPs, L)

    Is = Interval[]
    for (lp, rp) in zip(new_lPs, new_rPs)
        bs1, bs2 = get_blocks(1, lp), get_blocks(2, lp)
        @assert get_blocks(1, rp) == bs1 "$(get_blocks(1, rp)) !== $bs1"
        @assert get_blocks(2, rp) == bs2 "$(get_blocks(2, rp)) !== $bs2"
        i = Interval(lp, rp, L, vcat(bs1, bs2))
        push!(Is, i)
    end

    return Is
end

# create a dictionary (block uiid => (#block in ι for ι in isolates))
# julia> block_map(Gpj)
# Dict{String, Vector{Int64}} with 348 entries:
#   "QBKMHPTRYB" => [0, 1]
#   "SGXMRKQJZL" => [1, 1]
#   "IIUFJPIWIY" => [1, 1]
#   "AJGMEUFYLW" => [0, 1]
#   ...
function block_map(G::Graph)

    # define counter type
    BlockCounter = Dict{String,Int}

    # meta-counter
    Counters = Dict{String,BlockCounter}()

    # count block occurrences per isolate
    for (iso, seq) in G.sequence
        BC = BlockCounter()
        for n in seq.node
            id = n.block.uuid
            if id ∉ keys(BC)
                BC[id] = 0
            end
            BC[id] += 1
        end
        Counters[iso] = BC
    end

    # add missing ids to the counter (0-count blocks)
    for id in keys(G.block), iso in keys(G.sequence)
        if id ∉ keys(Counters[iso])
            Counters[iso][id] = 0
        end
    end
    return Counters
end


# given a list of intervals and a boolean function (Interval -> true/false) returns
# stats on the size of true and false sizes of intervals
function bool_function_length(Is::Vector{Interval}, fbool::Function)
    v = [fbool(i.bs...) for i in Is]
    l = [length(i) for i in Is]
    return Dict("L" => sum(l[v]), "N" => sum(v))
end

if abspath(PROGRAM_FILE) == @__FILE__

    # capture arguments
    path_pairwise, path_project, outfile = ARGS

    Gpw = open(unmarshal, path_pairwise)
    Gpj = open(unmarshal, path_project)

    BMw = block_map(Gpw)
    BMj = block_map(Gpj)


    # dictionary containing results, to later be exported
    results = Dict{String,Dict{String,Dict{String,Int64}}}()
    function pushres!(iso, kind, val)
        (iso ∉ keys(results)) && (results[iso] = Dict{String,Dict{String,Int64}}())
        results[iso][kind] = val
    end

    # dictionary (isolate => merged intervals)
    Im = Dict{String,Vector{Interval}}()

    # cycle through the two possible merging
    iso1, iso2 = collect(keys(Gpw.sequence))
    for (iso, niso) in [[iso1, iso2], [iso2, iso1]]

        # create interval merging
        Iw = to_intervals(Gpw, iso)
        Ij = to_intervals(Gpj, iso)
        Im[iso] = merge(Iw, Ij)

        # function that recover the counts of each block and applies 
        # a boolean function of it. Saves the results on the result dictionary
        function query!(count_f::Function, label::String)
            bool_f = (bw, bj) -> let
                cwi = BMw[iso][bw] # count of block on isolate (pairwise)
                cji = BMj[iso][bj] # count of block on isolate (projection)
                cwn = BMw[niso][bw] # count of block on other isolate (pairwise)
                cjn = BMj[niso][bj] # count of block on other isolate (projection)
                @assert cwi > 0
                @assert cji > 0
                return count_f(cwi, cji, cwn, cjn)
            end
            pushres!(iso, label, bool_function_length(Im[iso], bool_f))
        end

        # total n. of blocks and length
        query!((wi, ji, wn, jn) -> true, "total")

        # --- shared or private ---

        # fraction of private on pairwise
        query!((wi, ji, wn, jn) -> (wn == 0), "private")

        # fraction of private on projection
        query!((wi, ji, wn, jn) -> (jn == 0), "private")

        # fraction of agreed
        query!((wi, ji, wn, jn) -> ~((wn > 0) ⊻ (jn > 0)), "agree")

        # disagree
        query!((wi, ji, wn, jn) -> (wn > 0) ⊻ (jn > 0), "disagree")

        # private agreed
        query!((wi, ji, wn, jn) -> (wn == 0) & (jn == 0), "private agree")

        # shared agreed
        query!((wi, ji, wn, jn) -> (wn > 0) & (jn > 0), "shared agree")

        # shared only pairwise
        query!((wi, ji, wn, jn) -> (wn > 0) & (jn == 0), "shared only pairwise")

        # shared only projection
        query!((wi, ji, wn, jn) -> (wn == 0) & (jn > 0), "shared only projection")

        # --- duplicated or not ---

        # duplcated on pairwise
        query!((wi, ji, wn, jn) -> (wi > 1), "duplicated pairwise")

        # duplcated on projection
        query!((wi, ji, wn, jn) -> (ji > 1), "duplicated projection")

        # agree duplicated
        query!((wi, ji, wn, jn) -> ~((wi > 1) ⊻ (ji > 1)), "agree duplicated")

        # disagree duplicated
        query!((wi, ji, wn, jn) -> (wi > 1) ⊻ (ji > 1), "disagree duplicated")

        # both duplicated
        query!((wi, ji, wn, jn) -> (wi > 1) & (ji > 1), "both duplicated")

        # no duplicated
        query!((wi, ji, wn, jn) -> ~((wi > 1) & (ji > 1)), "both duplicated")

    end

    open(outfile, "w") do io
        JSON.print(io, results)
    end

end