struct Node
    block::Block
    number::Int
    strand::Bool
end

struct Path
    isolate::String
    nodes::Array{Node}
    offset::Integer
    circular::Bool
end
