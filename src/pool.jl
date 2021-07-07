module Pool

using ..Utility: random_id

import Base:
    open, isopen, close,
    write, unsafe_write

export FIFO, delete, path
export pool, shutdown

# ------------------------------------------------------------------------
# Named Pipes
# NOTE: assumes a posix environment!
# TODO: use pycall in the case for windows?
struct FIFO
    root::String
    base::String
    mode::Integer

    function FIFO(base, mode=0o666)
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

FIFO() = FIFO(random_id())

path(f::FIFO) = joinpath(f.root, f.base)

function delete(f::FIFO)
    err = ccall(:remove, Cint, (Cstring,), path(f))
    systemerror("failed to delete fifo at $(repr(f.base))", err != 0)
end

struct FIFOStream <: IO
    fd :: Int
end

const RDONLY = 0o0
const WRONLY = 0o1
const RDWR   = 0o2

function open(f::FIFO;
    read     :: Union{Bool,Nothing} = nothing,
    write    :: Union{Bool,Nothing} = nothing,
    create   :: Union{Bool,Nothing} = nothing,
    truncate :: Union{Bool,Nothing} = nothing,
    append   :: Union{Bool,Nothing} = nothing,
)
    # TODO: compute correct open flags
    fd = ccall(:open, Cint, (Cstring, Cint), path(f), WRONLY)
    fd < 0 && error("failed to open FIFO '$(path(f))'")

    return FIFOStream(fd)
end

function open(f::FIFO, mode::AbstractString)::FIFOStream
    mode == "r"  ? Base.open(f, read = true)                  :
    mode == "r+" ? Base.open(f, read = true, write = true)    :
    mode == "w"  ? Base.open(f, truncate = true)              :
    mode == "w+" ? Base.open(f, truncate = true, read = true) :
    mode == "a"  ? Base.open(f, append = true)                :
    mode == "a+" ? Base.open(f, append = true, read = true)   :
    throw(ArgumentError("invalid open mode: $mode"))
end

function open(func::Function, f::FIFO, mode::AbstractString)
    io = open(f, mode)
    try
        func(io)
    finally
        close(io)
    end
end

function close(io::FIFOStream)
    err = ccall(:close, Cint, (Cint,), io.fd)
    err != 0 && error("failed to close FIFOStream $(io.fd)")

    return
end

# write a single byte
function write(io::FIFOStream, x::UInt8)
    n = ccall(:write, Cssize_t, (Cint,Ptr{Cvoid},Csize_t), io.fd, Ref{UInt8}(x), 1)
    n != 1 && error("failed to write to FIFOStream $(io.fd), errno = $(Base.Libc.errno())")
    return 1
end

# write many bytes
function unsafe_write(io::FIFOStream, p::Ptr{UInt8}, n::UInt)
    m = ccall(:write, Cssize_t, (Cint,Ptr{Cvoid},Csize_t), io.fd, p, n)
    m != n && error("failed to write to FIFOStream $(io.fd), errno = $(Base.Libc.errno())")
    return Int(n)
end

# --------------------------------
# worker queues

function pool(size)
    chan = Channel{FIFO}(size)
    for i in 1:size
        f = FIFO()
        while isnothing(f)
            f = FIFO()
        end
        put!(chan, f)
    end

    return chan
end

function shutdown(fifos::Channel{FIFO})
    close(fifos)
    for f in fifos
        delete(f)
    end
end

end
