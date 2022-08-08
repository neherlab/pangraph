using JLD2
using JSON


if abspath(PROGRAM_FILE) == @__FILE__

    # filenames
    jld2_file = ARGS[1]
    json_file = ARGS[2]

    # load data
    data = load(jld2_file)

    # collect one entry per simulation
    jdata = Dict()
    for (k, v) in data
        # parse key of data as simulation parameters
        hgt, snps, trial, kind = split(k, "/")
        hgt = parse(Float64, hgt)
        snps = parse(Float64, snps)
        trial = parse(Int, trial)

        # create new dictionary key
        key = (hgt = hgt, snps = snps, trial = trial)
        # if key not present, add it
        ~(haskey(jdata, key)) && (jdata[key] = Dict{String,Any}())
        # append measurement
        jdata[key][kind] = v
    end

    # format as list of nested dictionaries, one per simulation
    jdata = [
        Dict("hgt" => k.hgt, "snps" => k.snps, "trial" => k.trial, "values" => v) for
        (k, v) in jdata
    ]

    #output as json file
    open(json_file, "w") do f
        JSON.print(f, jdata)
    end

end
