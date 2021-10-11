module Shell

export minimap2, mash, mafft

# command line execution
function execute(cmd::Cmd; now=true)
    out = Pipe()
    err = Pipe()

    proc = run(pipeline(cmd, stdout=out, stderr=err); wait=now)

    close(out.in)
    close(err.in)

    stdout = @async String(read(out))
    stderr = @async String(read(err))

    if now
        return (
            out  = fetch(stdout),
            err  = fetch(stderr),
            code = proc.exitcode, #err,
        )
    else
        return (
            out  = stdout,
            err  = stderr,
            proc = proc,
        )
    end
end

function mash(input)
    result = execute(`mash triangle $input`; now=false)
    stdout = IOBuffer(fetch(result.out))

    N     = parse(Int64,readline(stdout))
    dist  = zeros(N,N)
    names = Array{String}(undef, N)
    for (i, line) in enumerate(eachline(stdout))
        elt = split(strip(line))
        names[i] = elt[1]
        dist[i,1:(i-1)] = [parse(Float64,x) for x in elt[2:end]]
    end

    dist = dist + dist';

    return dist, names
end

function minimap2(qry::String, ref::String)
    return execute(`minimap2 -x asm10 -m 10 -n 1 -s 30 -D -c $ref $qry`; now=false)
end

function mafft(path::AbstractString)
    return execute(`mafft --auto --quiet --nuc $path`; now=false)
end

end
