module Command

using Random
using Base.Iterators

sequence(l::Int) = String(rand(['A','C','G','T'], l))
function line(s; col=80)
    n = length(s) ÷ col
    l =  ( s[(i-1)*col+1:(i*col)] for i in 1:n )

    length(s) % col == 0 && return join(l,"\n")
    return join(l,"\n") * "\n" * s[(n*col)+1:end]
end

function usage()
    println("usage: julia script/ancestors.jl -n [int] -l [int]")
    exit(2)
end

function main(args)
    N,L = 100, Int(1e5)
    i = 1
    while i ≤ length(args)
        if args[i] == "-N"
            N = parse(Int,args[i+1])
            i += 2
            continue
        end
        if args[i] == "-L"
            L = parse(Int,args[i+1])
            i += 2
            continue
        end
        break
    end

    for i in 1:N
        s = sequence(L)
        name = "isolate_$(i)"
        write(stdout, ">$(name)", '\n')
        write(stdout, line(s), '\n')
    end
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end
