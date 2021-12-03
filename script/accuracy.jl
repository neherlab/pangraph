module Accuracy

using PanGraph
using Hungarian

# ------------------------------------------------------------------------
# cost computation

δ(q, r; L=100) = mod(q-r, L)
distance(qry, ref; L=100) = δ.(qry, ref'; L=L)

function compare(path)
	graph = (
		known = open(Graphs.unmarshal,path.known),
		guess = open(Graphs.unmarshal,path.guess),
	)

	for (name,known) in graph.known.sequence
		guess = graph.guess.sequence[name]
		dists = distance(guess.position[1:end-1],known.position[1:end-1]; L=known.position[end])
        @show guess.position, known.position
        assignment, cost = hungarian(dists)
        @show cost
	end
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
