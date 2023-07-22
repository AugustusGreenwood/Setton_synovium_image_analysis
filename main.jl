using Images, DataStructures, MAT, BenchmarkTools




include("ImageIO.jl")
include("ImageProperties.jl")





function get_level_entropy(img::AbstractArray{Gray{T}})::Vector{Float64} where T
    ent = Vector{Float64}(undef, 256)
    g = glcm(img, [0 1]; symmetric=true)[[0,1]]
    g_norm = g / sum(g)
    ent[end] = entropy(g_norm)
    for i in 1:255
        g = glcm(img, [0 1]; levels=i, symmetric=true)[[0,1]]
        g_norm = g / sum(g)
        ent[i] = entropy(g_norm)
    end
    return ent
end

function get_entropy_error(img::AbstractArray{Gray{T}})::Vector{Float64} where T
    ent = get_level_entropy(img)

    error = Vector{Float64}(undef, 256)
    for (i, e) in enumerate(ent)
        error[i] = abs(ent[end] - e)
    end
    return error
end

function get_average_level_error_for_stack(stack::Array{Gray{T}, 3})::Vector{Float64} where T
    total_error = zeros(256)
    Threads.@threads for img in eachslice(stack, dims=3)
        total_error += get_entropy_error(img)
    end
    return total_error ./ size(image)[3]
end

function get_entropy_for_stack(stack::Array{Gray{T}, 3})::Vector{Float64} where T
    ent = Vector{Float64}(undef, size(stack)[3])
    for (i, img) in enumerate(eachslice(stack, dims=3))
        g = glcm(img, [0 8; 8 8; 8 0; 8 -8])
        ent[i] = entropy(g)
    end
    return ent
end
