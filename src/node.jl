module Nodes

import Base: show

export Node

# NOTE: we define an abstract type here to allow for mutually recursive types
#       it is expected that T = block for all purposes within the code
#       this is solely to deal with a limitation of julia
# NOTE: we make this mutable so that independent calls to the constructor
#       with equivalent parameters ensures the return of a brand new object
mutable struct Node{T}
    block  :: T
    strand :: Bool
end

# --------------------------------
# constructors

Node{T}(b::T; strand=true) where T = Node{T}(b,strand)

# --------------------------------
# operators

show(io::IO, n::Node) = Base.show(io, UInt64(pointer_from_objref(n)))

end
