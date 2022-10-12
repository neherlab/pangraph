using JLD2

if abspath(PROGRAM_FILE) == @__FILE__
    out = ARGS[1]
    data = ARGS[2:end]

    jldopen(out, "w") do database
        for d in data
            Δ = load(d)
            for (k, v) in Δ
                database[k] = v
            end
        end
    end
end
