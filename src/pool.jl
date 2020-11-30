module Pool

include("util.jl")
using .Utility: random_id

export Fifo, delete, path
export pool, shutdown

# ------------------------------------------------------------------------
# Named Pipes
# NOTE: assumes a posix environment!
# TODO: use pycall in the case for windows?
struct Fifo
    root::String
    base::String
    mode::Integer

    function Fifo(base, mode=0o666)
        root = tempdir()
        path = joinpath(root,base)
        if ispath(path)
            return nothing
        end

        err  = ccall(:mkfifo, Cint, (Cstring, Cuint), path, mode)
        systemerror("failed to make fifo at $(repr(path))", err != 0)

        return new(root, base, mode)
    end
end

path(f::Fifo) = joinpath(f.root, f.base)

function delete(f::Fifo)
    err = ccall(:remove, Cint, (Cstring,), path(f))
    systemerror("failed to delete fifo at $(repr(f.base))", err != 0)
end

Base.open(f::Fifo, mode::AbstractString; lock = true) = Base.open(path(f), mode; lock=lock)

# --------------------------------
# worker queues

function pool(size)
    chan = Channel{Fifo}(size)
    for i in 1:size
        f = Fifo(random_id())
        while isnothing(f)
            f = Fifo(random_id())
        end
        put!(chan, f)
    end

    return chan
end

function shutdown(fifos::Channel{Fifo})
    close(fifos)
    for f in fifos
        delete(f)
    end
end

end
