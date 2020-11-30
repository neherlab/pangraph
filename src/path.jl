module Paths

using ..Blocks

export Node, Path

struct Node
    block::Block
    number::Int
    strand::Bool
end

Node(b) = Node(b,0,true)

struct Path
    name::String
    nodes::Array{Node}
    offset::Integer
    circular::Bool
end

Path(name::String,node::Node,circular::Bool) = Path(name,[node],0,circular)
Path(name::String,node::Node) = Path(name,[node],0,false)

end
