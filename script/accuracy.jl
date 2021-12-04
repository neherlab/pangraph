module Accuracy

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

function compare(path)
	graph = (
		known = open(unmarshal,path.known),
		guess = open(unmarshal,path.guess),
	)

    costs = [
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

    @show mean(costs), std(costs)
    return costs
end

# ------------------------------------------------------------------------
# main point of entry

function usage()
	println("usage: julia script/accuracy.jl <known.json> <guess.json>")
	exit(2)
end

if abspath(PROGRAM_FILE) == @__FILE__
	length(ARGS) == 2 || usage()

	compare((known=ARGS[1], guess=ARGS[2]))
end

end
