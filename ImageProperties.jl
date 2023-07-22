using DataStructures

struct ImageProperties
    contrast::Float64
    correlation::Float64
    energy::Float64
    homogeneity::Float64
    entropy::Float64
    glcm_norm::Matrix{Float64}
    function ImageProperties(glcm::Matrix{Int64})
        rows, cols = size(glcm)
        glcm_norm = glcm / sum(glcm)
        cont = contrast(glcm_norm, rows, cols)
        corr = correlation(glcm_norm, rows, cols)
        energ = energy(glcm_norm)
        homogen = homogeneity(glcm_norm, rows, cols)
        entropy = 0
        return new(cont, corr, energ, homogen, entropy, glcm_norm)
    end
end



# This is based upon how matlab scaled their gray scale images
function scale_image(image::AbstractArray{Gray{T}}, levels::Int64)::Matrix{Int64} where T
    (mini, maxi) = minimum(image), maximum(image)
    slope = levels / (maxi - mini)
    eq(x) = slope * x + (1 - slope * mini)
    scaled_image =  @. floor(eq(image))
    scaled_image[scaled_image .> levels] .= levels
    scaled_image[scaled_image .< 1] .= 1
    return Int64.(scaled_image)
end


function __get_row_col_range(rows::Int64, cols::Int64, drow::Int64, dcol::Int64)::Tuple{UnitRange{Int64}, UnitRange{Int64}}
    if drow < 0
        row_range = 1-drow:rows
    else
        row_range = 1:rows-drow
    end
    
    if dcol < 0
        col_range = 1-dcol:cols
    else
        col_range = 1:cols-dcol
    end

    return row_range, col_range
end


# GLCM calls are multithreaded, but this may not be best depending on the situation, I'm leaving it because
# as of right now, the benefits when it helps are much greater than the negatives when it doesn't
function asymmetric_glcm(scaled_image::Matrix{Int64}, offset::AbstractVector{Int64}, levels::Int64)::Matrix{Int64}
    drow, dcol = offset
    n, m = size(scaled_image)
    row_range, col_range = __get_row_col_range(n, m, drow, dcol)

    glcm_mat = zeros(Int64, levels, levels)
    Threads.@threads for i in row_range
        for j in col_range
            pix = @views scaled_image[i,j]
            nei = @views scaled_image[i+drow,j+dcol]
            glcm_mat[pix, nei] += 1
        end
    end
    return glcm_mat
end


# I think that a symmetric_glcm can be calculated by adding the transpose? Don't know, but this is 
# fast enough so I probably wont change it
function symmetric_glcm(scaled_image::Matrix{Int64}, offset::AbstractVector{Int64}, levels::Int64)::Matrix{Int64}
    drow, dcol = offset
    n, m = size(scaled_image)
    row_range, col_range = __get_row_col_range(n, m, drow, dcol)

    glcm_mat = zeros(Int64, levels, levels)
    Threads.@threads for i in row_range
        for j in col_range
            pix = @views scaled_image[i,j]
            nei = @views scaled_image[i+drow,j+dcol]
            glcm_mat[pix, nei] += 1
            glcm_mat[nei, pix] += 1
        end
    end
    return glcm_mat
end


# Call this function, it will take care of the other stuff that asymmetric_glcm/symmetric_glcm dont
function glcm(image::AbstractArray{Gray{T}}, offsets::AbstractArray{Int64}; levels::Int64=0, symmetric::Bool=true)::Dict{Vector{Int64}, Matrix{Int64}} where T
    scaled_image = scale_image(image, levels)
    
    all_glcms = Dict{Vector{Int64}, Matrix{Int64}}()
    if symmetric == true
        for offset in eachrow(offsets)
            all_glcms[offset] = @views symmetric_glcm(scaled_image, offset, levels)
        end
    else
        for offset in eachrow(offsets)
            all_glcms[offset] = @views asymmetric_glcm(scaled_image, offset, levels)
        end
    end
    return all_glcms
end



# These aren't 'properties' so they are prefaced with "calculate". Does it make sense? No. Do I care? Yes. 
# Will I change it? Probably not.
function calculate_mu(glcm_norm::Matrix{T}, rows::Int64, cols::Int64)::Tuple{Float64, Float64} where T
    μi = 0
    μj = 0
    for i in 1:rows
        for j in 1:cols
            μi += i * glcm_norm[i,j]
            μj += j * glcm_norm[i,j]
        end   
    end
    return μi, μj
end

function calculate_sigma(glcm_norm::Matrix{T}, rows::Int64, cols::Int64, μi::Float64, μj::Float64)::Tuple{Float64, Float64} where T
    sigmai = 0
    sigmaj = 0
    for i in 1:rows
        for j in 1:cols
            sigmai += (i - μi)^2 * glcm_norm[i,j]
            sigmaj += (j - μj)^2 * glcm_norm[i,j]
        end
    end
    return sqrt(sigmai), sqrt(sigmaj)
end


# All below are properties which have been validated against matlab outputed values. EXCEPT FOR ENTROPY
function contrast(glcm_norm::Matrix{T}, rows::Int64, cols::Int64)::Float64 where T
    cont = 0
    for i in 1:rows
        for j in 1:cols
            cont += (i - j)^2 * glcm_norm[i,j]
        end
    end
    return cont
end


function correlation(glcm_norm::Matrix{T}, rows::Int64, cols::Int64)::Float64 where T
    μi, μj = calculate_mu(glcm_norm, rows, cols)
    sigmai, sigmaj = calculate_sigma(glcm_norm, rows, cols, μi, μj)
    
    corr = 0
    for i in 1:rows
        for j in 1:cols
            corr += (i - μi) * (j - μj) * glcm_norm[i,j] / (sigmai * sigmaj)
        end
    end
    return corr
end


function energy(glcm_norm::Matrix{T})::Float64 where T
    return sum(glcm_norm.^2)
end


function homogeneity(glcm_norm::Matrix{T}, rows::Int64, cols::Int64)::Float64 where T
    homogen = 0
    for i in 1:rows
        for j in 1:cols
            homogen += glcm_norm[i,j] / (1 + abs(i - j))
        end
    end
    return homogen
end



## NOT VALIDATED
function ent(glcm_norm::Matrix{T})::Float64 where T
    prob_dict = DefaultDict{T, Int64}(1)
    for pixel in glcm_norm
        prob_dict[pixel] += 1
    end
    p = values(prob_dict) ./ length(prob_dict)
    return -sum(p .* log.(p))
end








