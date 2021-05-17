module Paths

using ..Nodes
using ..Blocks

import ..Graphs: 
    pair, sequence,
    Counter, add!

export Path
export count_isolates

struct Path
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
Base.show(io::IO, p::Path) = Base.show(io, (name=p.name, blocks=p.node))

function sequence(p::Path)
    seq = join(sequence(node.block, node) for node ∈ p.node)
    if p.offset !== nothing
        seq = String(circshift(Array{UInt8,1}(seq), p.offset))
    end
    return seq
end

# used for block merging
function Base.replace!(p::Path, old::Block, new::Array{Block})
    indices = Int[]
    inserts = Array{Node{Block}}[]
    for (i, n₁) in enumerate(p.node)
        if n₁.block != old
            continue
        end

        push!(indices, i)

        nodes = (n₁.strand ? [Node{Block}(nb) for nb in new] 
                           : [Node{Block}(nb;strand=false) for nb in reverse(new)])
        for (blk,n₂) in zip(n₁.strand ? new : reverse(new), nodes)
            swap!(blk, n₁, n₂)
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
                    # XXX: should never hit this block
                    error("bad interval match")
                    return (interval=nothing,strand=nothing)
                end
                blk != p.node[x].block && return (interval=nothing,strand=nothing)
                stop, x = x, advance(x)
            end

            if parity
                stop ≥ start && return (interval=start:stop, strand=true)           # simple case: | -(--)- |
                !p.circular  && error("invalid circular interval on linear path")   # broken case: | -)--(- |

                #= DEBUG
                print("forward: [")
                for i ∈ stop:length(p.node)
                    print(stderr, p.node[i].block, ", ")
                end
                for i ∈ 1:start
                    print(stderr, p.node[i].block, ", ")
                end
                println("]")
                =#

                return (interval=(start:length(p.node), 1:stop), strand=true)
            else
                stop ≤ start && return (interval=stop:start, strand=false)          # simple case: | -(--)- |
                !p.circular  && error("invalid circular interval on linear path")   # broken case: | -)--(- |

                #= DEBUG
                print("reverse: [")
                for i ∈ start:-1:1
                    print(stderr, p.node[i], ", ")
                end
                for i ∈ length(p.node):-1:stop
                    print(stderr, p.node[i], ", ")
                end
                println("]")
                =#

                return (interval=(1:start, stop:length(p.node)), strand=false)
            end
        end
    )

    swap!(b, i::AbstractArray, new)                       = Blocks.swap!(b, p.node[i], new)
    swap!(b, i::Tuple{AbstractArray,AbstractArray}, new)  = Blocks.swap!(b, [p.node[i[1]]; p.node[i[2]]], new)

    splice!(nodes, i::AbstractArray, new)                       = Base.splice!(nodes, i, [new])
    splice!(nodes, i::Tuple{AbstractArray,AbstractArray}, new)  = let 
        Base.splice!(nodes, i[1], [new])
        Base.splice!(nodes, i[2])
    end

    for (i,s) ∈ zip(interval, strand)
        n = Node(new, s)

        swap!(new, i, n)
        splice!(p.node, i, n)
    end

    #=
    OLD CODE:
    start, stop = old[1].block, old[end].block
    i₁ₛ = findall((n) -> n.block == start, p.node)
    i₂ₛ = findall((n) -> n.block == stop,  p.node)

    @assert length(i₁ₛ) == length(i₂ₛ)

    # XXX: the ordering here is wrong for circular genomes...
    iₛ = p.circular ? intervals(i₁ₛ, i₂ₛ, length(old), length(p.node)) : intervals(i₁ₛ, i₂ₛ)
    sₛ = map(i₁ₛ) do i
        p.node[i].strand == old[1].strand
    end


    # NOTE: each section is reversed so that we modify the larger indices first
    for (i,s) ∈ zip(reverse(iₛ), reverse(sₛ))
        node = Node(new, s)

        @show i,s
        swap!(new, i, node)
        splice!(p.node, i, [node])
    end
    =#
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
