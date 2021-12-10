using Plots

cdfplot(x; kwargs...)  = plot(sort(x),range(0,1,length=length(x)); kwargs...)
cdfplot!(x; kwargs...) = plot!(sort(x),range(0,1,length=length(x)); kwargs...)

Base.zero(x::Type{Array{Float64,1}}) = Float64[]

function unpack(key)
    elt = split(key,'/')
    return (
       hgt=parse(Float64,elt[1]),
       snp=parse(Float64,elt[2]),
       nit=parse(Int64,elt[3]),
       elt=elt[4],
    )
end

function collectentry(data, type, name, accumulate!, post!)
    hgt, snp = params(data)
    val = [ zero(type) for i in 1:length(hgt), j in 1:length(snp) ]
    num = zeros(Int,length(hgt),length(snp))
    for key in keys(data)
        param = unpack(key)
        param.elt == name || continue

        i = findfirst(hgt .== param.hgt)
        j = findfirst(snp .== param.snp)

        data[key] !== nothing || continue
        accumulate!(val,i,j,data[key])
        num[i,j] += 1
    end
    return post!(val, num), hgt, snp
end

# specific getters
entropy(data)    = collectentry(data,Float64,"tiles",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)
accuracy(data)   = collectentry(data,Array{Float64,1},"costs",(arr,i,j,x)->append!(arr[i,j],x),(arr,num)->arr)
diversity(data)  = collectentry(data,Float64,"dists",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)
complexity(data) = collectentry(data,Float64,"nblks",(arr,i,j,x)->arr[i,j]+=x, (arr,num)->arr./num)
