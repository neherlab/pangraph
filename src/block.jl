using Match, FStrings
using Random, StatsBase

# ------------------------------------------------------------------------
# Tag data structure

struct Tag
    isolate::String
    number::Int
end

# ------------------------------------------------------------------------
# Block data structure

struct Block
    uuid::String
    sequence::Array{Char}
    mutation::Dict{Tag,Dict{Int,Char}}
    indel::Dict{Tag,Dict{Int,Union{Array{Char},Int}}}
end

# ---------------------------
# constructors

# simple helpers
Block() = Block(id(), "", {}, {}, {})
Block(name::String, sequence::Array{Char}) = Block(id(),sequence, {Tag(name,0):{}}, {Tag(name,0):{}})
Block(sequence,mutation,indel) = Block(id(),sequence,mutation,indel)

# serial concatenate list of blocks
function Block(bs::Block...)
    @assert all([isolates(bs[1]) == isolates(b) for b in bs[2:end]])

    sequence = join([b.sequence for b in bs])
    mutation = bs[1].mutation
    indel    = bs[1].insertion

    # TODO: finish the concatenation of mutations/etc dictionaries
    #       need to thread degeneracies through the list correctly
    δ = length(bs[1])
    for b in bs[2:end]
        δ += length(b)
    end


    return Block(sequence,mutation,indel)
end

# returns a subslice of block b
function Block(b::Block, range::UnitRange{Int})
    @assert range.start >= 1 && range.stop <= length(b)
    translate(dict) = Dict(iso => Dict(x-range.start => val for (x,val) in d) for (iso,d) in dict)
    sequence = b.sequence[range]
    mutation = translate(b.mutation)
    indel    = translate(b.indel)

    return Block(sequence,mutation,indel)
end

# ---------------------------
# operations

# simple operations
length(b::Block) = length(b.sequence)
depth(b::Block)  = length(b.mutation)

# complex operations

# returns the count of isolates within the block
function isolates(b::Block)
    count = Dict{String,Int}()
    for tag in b.mutation
        if tag.isolate in keys(count)
            count[tag.isolate] += 1
        else
            count[tag.isolate]  = 1
        end
    end

    return count
end

# returns the sequence WITH mutations and indels applied to the consensus for a given tag 
function sequence(b::Block, tag::Tag; gaps=false)
    seq = copy(b.sequence)
    # mutations
    for (locus, snp) in b.mutation[tag]
        seq[locus] = snp
    end
    # indels
    # TODO: deal with case of gaps=true and return explicit gap characters!
    for locus in sort(keys(b.indel),rev=true)
        indel = b.indel[locus]
        @match typeof(indel) begin
            Array{Char} => seq = cat(seq[1:locus], indel, seq[locus+1:end])
            Int         => seq = cat(seq[1:locus], seq[locus+indel:end])
            _ => error(f"unrecognized type of indel = {typeof(indel)}") 
        end
    end

    return seq
end

function test()
    println(id())
end
