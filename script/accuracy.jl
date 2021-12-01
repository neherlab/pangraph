module Accuracy

using PanGraph

# ------------------------------------------------------------------------
# munkres assignment algorithm

struct Index
	row :: Int
	col :: Int
end
Index() = Index(0,0)
valid(i::Index) = i.row > 0 && i.col > 0

struct Cost{T <: Real}
	matrix :: Matrix{T}
	δrow   :: Vector{T}
	δcol   :: Vector{T}
end

function Cost(matrix::Matrix{T}) where T <: Real
	return Cost{T}(
			matrix,
			zeros(T, size(matrix,1)),
			zeros(T, size(matrix,2)),
	)
end

# method extensions
Base.eltype(cost::Cost) = Base.eltype(cost.matrix)
Base.size(cost::Cost, args...) = Base.size(cost.matrix, args...)

function Base.getindex(cost::Cost, i, j)
	@inbounds return cost.matrix[i,j] - cost.δrow[i] - cost.δcol[j]
end

function iszero(cost::Cost, i, j)
	@inbounds return cost.matrix[i,j] == cost.δrow[i] + cost.δcol[j]
end

const Mark = (
	Star  = 1,
	Prime = 2
)

# remove minimum from each row
function step1!(cost)
	cost.δrow = vec(minimum(cost.matrix, dims=2))
	return [findall(j->iszero(cost,i,j), 1:size(cost,2)) for i in 1:size(cost,1)]
end

# find any singleton row zeros
function step2!(cost, mask, cover)
	for i in 1:size(cost,1)
		for j in 1:size(cost,2)
			if iszero(cost,i,j) && !cover.row[i] && !cover.col[j]
				mask[i,j] = Mark.Star
				cover.row[i] = true
				cover.col[j] = true
			end
		end
	end

	cover.row .= false
	cover.col .= false
end

# cover each column with a starred zero. if everything is covered, we are done
function step3!(mask, cover)
	n, m = size(mask)

	for i in 1:n
		for j in 1:m
			if mask[i,j] == Mark.Star
				cover.col[j] = true
			end
		end
	end

	count = sum(cover.col)
	return ((count >= m) || (count >= n)) ? 7 : 4
end

function step4!(cost, mask, cover, zeros)
end

function step5!(mask, cover, start)
end

function step6!(mask, cover, start)
end

function assign(data)
	n, m = size(data)
	flip = false
	if n > m
		data = Array(data')
		flip = true
		n, m = size(data)
	end

	cost  = Cost(data)
	mask  = zeros(Int8,n,m)
	cover = (
		row = zeros(Bool,n),
		col = zeros(Bool,m),
	)

	zeros = step1!(cost)


	flip && return [findfirst(mask[:,i] .== Mark.Star) for i in 1:size(mask,2)]
	return [findfirst(mask[i,:] .== Mark.Star) for i in 1:size(mask,1)]
end

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
		D = distance(guess.position[1:end-1],known.position[1:end-1]; L=known.position[end])
	end
end

# ------------------------------------------------------------------------
# main point of entry

function usage()
	println("usage: julia src/accuracy.jl <known.json> <guess.json>")
	exit(2)
end

if abspath(PROGRAM_FILE) == @__FILE__
	length(ARGS) == 2 || usage()

	compare((known=ARGS[1], guess=ARGS[2]))
end

end
