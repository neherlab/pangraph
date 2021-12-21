module Nodes

import Base: show, length

export Node

# NOTE: we define an abstract type here to allow for mutually recursive types
#       it is expected that T = block for all purposes within the code
#       this is solely to deal with a limitation of julia
# NOTE: we make this mutable so that independent calls to the constructor
#       with equivalent parameters ensures the return of a brand new object
"""
	mutable struct Node{T}
		block  :: T
		strand :: Bool
	end

Node represents a portion of a sequence path that passes through a single `block`.
`strand` stores whether we pass along the forward strand of `block` (if true) or reverse (if false).
"""
mutable struct Node{T}
    block  :: T
    strand :: Bool
end

# --------------------------------
# constructors

"""
	Node{T}(b::T; strand=true)

Create a Node that passed through block `b`.
Default to forward strand orientation.
"""
Node(b::T; strand=true) where T = Node{T}(b,strand)

# --------------------------------
# operators

# Base.show(io::IO, n::Node) = Base.show(io, UInt64(pointer_from_objref(n)))
Base.show(io::IO, n::Node) = Base.show(io,
    (block = n.block, strand = n.strand)
)

"""
	length(n::Node) = length(n.block, n)

Return the length of sequence stored within node `n`
"""
length(n::Node) = length(n.block, n)

end
