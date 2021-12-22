cdfplot(x; kwargs...)  = plot(sort(x),range(0,1,length=length(x)); kwargs...)
cdfplot!(x; kwargs...) = plot!(sort(x),range(0,1,length=length(x)); kwargs...)
