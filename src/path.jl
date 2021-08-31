module Paths

import Base:
    length, show

using Infiltrator

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

function sequence(p::Path; shift=true)
    seq = join(String(sequence(n.block, n)) for n ∈ p.node)
    if p.offset !== nothing && shift
        seq = String(circshift(Array{UInt8,1}(seq), p.offset))
    end
    return seq
end

# used for block merging
function Base.replace!(p::Path, old::Block, new::Array{Block}, orientation::Bool)
    indices = Int[]
    inserts = Array{Node{Block}}[]

    for (i, n₁) in enumerate(p.node)
        n₁.block != old && continue

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
    # ----------------------------
    # internal functions
    
    unzip(a) = map(x->getfield.(a,x), fieldnames(eltype(a)))

    next = p.circular ? (x) -> (mod(x-0,length(p.node)) + 1) : (x) -> (x == length(p.node) ? nothing : x+1)
    prev = p.circular ? (x) -> (mod(x-2,length(p.node)) + 1) : (x) -> (x == 1 ? nothing : x-1)

    oldnodes(i::A, s::Bool)          where A <: AbstractArray = s ? p.node[i] : reverse(p.node[i])
    oldnodes(i::Tuple{A,A}, s::Bool) where A <: AbstractArray = s ? [p.node[i[1]]; p.node[i[2]]] : reverse([p.node[i[1]]; p.node[i[2]]])

    pack(i::A, s) where A <: AbstractArray = [(loci=i, strand=s, oldnode=oldnodes(i,s))]
    pack(i::Tuple{A,A}, s) where A <: AbstractArray = let 
        i₁, i₂ = (minimum(i[1]) < minimum(i[2])) ? (1, 2) : (2, 1)

        Δ = sum(length(n.block, n) for n in p.node[i[i₂]])
        if p.offset === nothing
            p.offset = -Δ
        else
            p.offset -= Δ
        end

        return [(loci=i[i₁], strand=s, oldnode=oldnodes(i,s)), (loci=i[i₂], strand=nothing, oldnode=nothing)]
    end

    # ----------------------------
    # function body
    
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
                    error("bad interval match")
                end
                if p.node[x].strand != (parity ? s : ~s)
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

                # @infiltrate
                # error("REVERSE")

                return (interval=(stop:length(p.node), 1:start), strand=false)
            end
        end
    )

    data = sort([x for (i,s) in zip(interval, strand) for x in pack(i,s)]; by=(x)->minimum(x.loci), rev=true)

    for datum ∈ data
        if datum.oldnode !== nothing
            newnode = Node(new, datum.strand)

            splice!(p.node, datum.loci, [newnode])
            swap!(new, datum.oldnode, newnode)
        else
            splice!(p.node, datum.loci, [])
        end

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
