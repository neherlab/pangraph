module Paths

import Base:
    length, show

using ..Nodes
using ..Blocks

import ..Graphs: 
    pair, sequence, reverse_complement,
    Counter, add!

export Path
export count_isolates

mutable struct Path
    name     :: String
    node     :: Array{Node{Block}}
    offset   :: Union{Int,Nothing}
    circular :: Bool
end

# --------------------------------
# constructors

Path(name::String,node::Node{Block};circular::Bool=false) = Path(name,[node],circular ? 0 : nothing,circular)

# --------------------------------
# operators

pair(p::Path) = p.name => p
show(io::IO, p::Path) = Base.show(io, (name=p.name, blocks=p.node))
length(p::Path) = length(p.node)

function sequence(p::Path)
    seq = join(String(sequence(n.block, n)) for n ∈ p.node)
    if p.offset !== nothing
        seq = String(circshift(Array{UInt8,1}(seq), p.offset))
    end
    return seq
end

# used for block merging
function Base.replace!(p::Path, old::Block, new::Array{Block}, orientation::Bool)
    indices = Int[]
    inserts = Array{Node{Block}}[]

    oldseq = sequence(p)

    for (i, n₁) in enumerate(p.node)
        n₁.block != old && continue

        @show n₁.block.mutate[n₁]
        @show n₁.block.insert[n₁]
        @show n₁.block.delete[n₁]

        push!(indices, i)

        nodes = ((n₁.strand==orientation) 
                     ? [Node{Block}(nb;strand=true) for nb in new] 
                     : [Node{Block}(nb;strand=false) for nb in reverse(new)])

        for n₂ in nodes
            swap!(n₂.block, n₁, n₂)
        end

        push!(inserts, nodes)
    end

    # reverse so that lower indices don't shift upper while we iterate
    reverse!(indices)
    reverse!(inserts)

    for (index, nodes) in zip(indices, inserts)
        splice!(p.node, index, nodes)
    end

    newseq = sequence(p)
    if oldseq != newseq
        badloci = Int[]
        for i ∈ 1:min(length(newseq),length(oldseq))
            if newseq[i] != oldseq[i]
                push!(badloci, i)
            end
        end
        left, right = max(badloci[1]-10, 1), min(badloci[1]+10, length(newseq))

        @show orientation
        @show indices
        println([(n.block.uuid,n.strand) for n in p.node])

        println("--> length:           ref($(length(oldseq))) <=> seq($(length(newseq)))")
        println("--> number of nodes:  $(length(p.node))")
        println("--> cumulative len:   $(cumsum([length(n.block,n) for n in p.node]))")
        println("--> offset:           $(p.offset)")
        println("--> window:           $(left):$(badloci[1]):$(right)")
        println("--> old:              $(oldseq[left:right])") 
        println("--> new:              $(newseq[left:right])") 

        node = p.node[indices[1]]

        @show node.block.mutate[node]
        @show node.block.insert[node]
        @show node.block.delete[node]

        @show new

        error("bad splicing")
    end
end

# XXX: should we create an interval data structure that unifies both cases?
# XXX: wrap as an iterator instead of storing the whole array in memory?
intervals(starts, stops) = [if start ≤ stop 
                                start:stop 
                            else 
                                stop:start 
                            end for (start,stop) in zip(starts,stops)]

function intervals(starts, stops, gap, len)
    if (stops[1]-starts[1]) == gap
        return intervals(starts, stops)
    elseif (len-starts[end]+stops[1]+1) == gap
        return [(starts[end]:len; 1:stops[1]);[start:stop for (start,stop) in zip(starts[1:end-1], stops[2:end])]]
    else
        error("unrecognized gap pattern in block replacement of path") 
    end
end

# used for detransitive
const Link = NamedTuple{(:block, :strand), Tuple{Block, Bool}}
function Base.replace!(p::Path, old::Array{Link}, new::Block)
    unzip(a) = map(x->getfield.(a,x), fieldnames(eltype(a)))

    next = p.circular ? (x) -> (mod(x-0,length(p.node)) + 1) : (x) -> (x == length(p.node) ? nothing : x+1)
    prev = p.circular ? (x) -> (mod(x-2,length(p.node)) + 1) : (x) -> (x == 1 ? nothing : x-1)

    matches = findall((n)->n.block == old[1].block, p.node)

    interval, strand = unzip(
        map(matches) do start
            parity  = p.node[start].strand == old[1].strand
            advance = parity ? next : prev

            stop, x = start, advance(start)
            for (blk, s) ∈ old[2:end]
                if x === nothing 
                    error("ran off genome")
                end
                if blk != p.node[x].block 
                    @show old
                    @show x, start, old[1].strand, parity
                    @show (p.node[x].block, p.node[x].strand)
                    @show (p.node[start].block, p.node[start].strand)
                    error("bad interval match")
                end
                if p.node[x].strand != (parity ? s : ~s)
                    @show old
                    @show x, start, old[1].strand, parity
                    @show (p.node[x].block, p.node[x].strand)
                    @show (p.node[start].block, p.node[start].strand)
                    error("bad strandedness")
                end
                stop, x = x, advance(x)
            end

            if parity
                stop ≥ start && return (interval=start:stop, strand=true)           # simple case: | -(--)- |
                !p.circular  && error("invalid circular interval on linear path")   # broken case: | -)--(- |

                return (interval=(start:length(p.node), 1:stop), strand=true)
            else
                stop ≤ start && return (interval=stop:start, strand=false)          # simple case: | -(--)- |
                !p.circular  && error("invalid circular interval on linear path")   # broken case: | -)--(- |

                error("REVERSE")

                return (interval=(1:start, stop:length(p.node)), strand=false)
            end
        end
    )

    oldnodes(i::AbstractArray)                      = p.node[i]
    oldnodes(i::Tuple{AbstractArray,AbstractArray}) = [p.node[i[1]]; p.node[i[2]]]

    splice!(nodes, i::AbstractArray, new) = let
        lengths = [length(n.block, n) for n in nodes]
        Base.splice!(nodes, i, [new])
    end
    splice!(nodes, i::Tuple{AbstractArray,AbstractArray}, new) = let 
        Δ = sum(length(n.block, n) for n in p.node[i[1]])
        if p.offset === nothing
            p.offset = -Δ
        else
            p.offset -= Δ
        end

        lengths = [length(n.block, n) for n in nodes]

        Base.splice!(nodes, i[1])
        Base.splice!(nodes, i[2], [new])
    end

    for (i,s) ∈ zip(interval, strand)
        newnode = Node(new, s)
        oldnode = oldnodes(i)

        # oldseq = join(String(sequence(n.block,n)) for n in oldnode)

        splice!(p.node, i, newnode)
        swap!(new, oldnode, newnode)

        # newseq = String(sequence(newnode.block,newnode))

        #=
        if newseq != oldseq
            error("FAIL")
        end
        =#
    end
end

# XXX: store as field in block?
#      conversely we could store a pointer to the Path in each node.
#      this would allow us to access which isolate each node belongs to.
#      for now we keep it global to have memory usage be minimal.
#      important to consider
# XXX: strongly type iterator for more specificity?
function count_isolates(paths)
    blocks  = Dict{Block, Counter}()

    for path in paths
        for node in path.node
            if node.block ∈ keys(blocks)
                add!(blocks[node.block], path.name)
            else
                blocks[node.block] = Counter([path.name => 1])
            end
        end
    end

    return blocks
end

end
