module Paths

using ..Nodes
using ..Blocks

import ..Graphs: pair

export Path

struct Path
    name::String
    node::Array{Node{Block}}
    offset::Union{Int,Nothing}
    circular::Bool
end

# --------------------------------
# constructors

Path(name::String,node::Node{Block};circular::Bool=false) = Path(name,[node],circular ? 0 : nothing,circular)

# --------------------------------
# operators

pair(p::Path) = p.name => p

end
