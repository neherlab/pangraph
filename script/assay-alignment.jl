module Accuracy

using JLD2

using PanGraph
using PanGraph.Graphs
using Hungarian
using Statistics

# ------------------------------------------------------------------------
# cost computation

δ(q, r; L=100) = mod(q-r, L)
distance(qry, ref; L=100) = min.(δ.(qry, ref'; L=L), δ.(ref', qry; L=L))

function cutoff(position, len, low; match=nothing)
    measurable = copy(position)

    @label loop
    D = mod.(measurable.- measurable', len)
    i = findfirst(0 .< D .< low)
    i !== nothing || @goto endloop

        if match === nothing
            deleteat!(measurable, i.I[2])
        else
            m1 = minimum(distance(measurable[i.I[1]], match; L=len))
            m2 = minimum(distance(measurable[i.I[2]], match; L=len))
            if m1 < m2
                deleteat!(measurable, i.I[2])
            else
                deleteat!(measurable, i.I[1])
            end
        end

    @goto  loop
    @label endloop

    return measurable
end

function nearestbreaks(graph)
    return [
        let
            L = length(sequence(known))
            guess = graph.guess.sequence[name]

            x = cutoff(guess.position,L,100)
            y = cutoff(known.position,L,100; match=x)
            d = distance(x,y; L=L)

            _, cost = hungarian(d)
            cost / min(length(x),length(y))
        end for (name,known) in graph.known.sequence
    ]
end

struct Partition
    dx :: Int64
    b1 :: Graphs.Block
    b2 :: Graphs.Block
end

function breakpoints(known, guess)
    L = length(sequence(known))
    #  x1->     x2->     x3->
    #  | node 1 | node 2 |
    bps = [(x,j,p.node[i].block,p.node[i]) for (j,p) in enumerate([known,guess]) for (i,x) in enumerate(p.position[1:end-1])]
    sort!(bps, by=first)

    # get the start block
    blocks = Array{Graphs.Block,1}(undef, 2)
    i, j = (bps[1][2] == 1) ? (1, 2) : (2, 1)

    blocks[i] = bps[1][3]
    for i in length(bps)-1:-1:1
        if bps[i][2] == 2
            blocks[j] = bps[i][3]
            break
        end
    end

    # edge case, if i = 1 and i = 2 correspond to equal break points
    part = Partition[]; i = 2
    if bps[i][1] == bps[i-1][1]
        blocks[bps[i][2]] = bps[i][3]
        i += 1
    end

    while i <= length(bps)
        dx = bps[i][1] - bps[i-1][1]
        push!(part, Partition(dx, blocks[1], blocks[2]))
        if (i <= length(bps)-1) && (bps[i][1] == bps[i+1][1])
            blocks[bps[i+0][2]] = bps[i+0][3]
            blocks[bps[i+1][2]] = bps[i+1][3]
            i += 2
        else
            blocks[bps[i][2]] = bps[i][3]
            i += 1
        end
    end
    # edge case: have to close the circular genome manually
    dx = δ(bps[1][1], bps[end][1]; L=L)
    push!(part, Partition(dx, blocks[1], blocks[2]))

    return part
end

function tiling(breaks)
    map = Dict{Tuple{Graphs.Block,Graphs.Block}, Int64}()
    for bp in breaks
        for p in bp
            key = (p.b1, p.b2)
            if key in keys(map)
                map[key] += p.dx
            else
                map[key] = p.dx
            end
        end
    end
    return map
end

function mutualentropy(graph)
    self  = (breakpoints(known, known) for (name, known) in graph.known.sequence if name == "isolate_2")
    cross = (breakpoints(known, graph.guess.sequence[name]) for (name, known) in graph.known.sequence if name == "isolate_2")

    function entropy(tile)
        len = sum(values(tile))
        return -sum(v/len * log(v/len) for v in values(tile))
    end

    return entropy(tiling(cross)) - entropy(tiling(self))
end

function compare(path)
	graph = (
		known = open(unmarshal,path.known),
		guess = open(unmarshal,path.guess),
	)

    μ = [Graphs.Blocks.diversity(b) for b in values(graph.guess.block)]
    l = [Graphs.Blocks.length(b) for b in values(graph.guess.block)]
    return (
        nearestbreaks(graph),
        mutualentropy(graph),
        sum(μ.*l) ./ sum(l),
        graph
    )
end

# ------------------------------------------------------------------------
# main point of entry

function usage()
	println("usage: julia script/assay-alignment.jl <known.json> <guess.json> <output.jld2>")
	exit(2)
end

if abspath(PROGRAM_FILE) == @__FILE__
	length(ARGS) == 5 || usage()

	costs, tiles, dists, input = compare((known=ARGS[1], guess=ARGS[2]))

    hgt, snp, output = ARGS[3:5]

    jldopen(output, "a+") do database
        group = JLD2.Group(database, "$(hgt)/$(snp)")
        # group["input"] = input # NOTE: increases the stored data massively
        group["costs"] = costs
        group["tiles"] = tiles
        group["dists"] = dists
    end
end

end
