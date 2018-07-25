module ChebyshevNodes

export chebTzeros, chebUzeros, chebVzeros, chebWzeros,
       chebTzeros01, chebUzeros01, chebVzeros01, chebWzeros01,
       chebTextrema, chebUextrema, chebVextrema, chebWextrema,
       chebTextrema01, chebUextrema01, chebVextrema01, chebWextrema01
       chebUUextrema, chebVVextrema, chebWWextrema,
       chebUUextrema01, chebVVextrema01, chebWWextrema01




#   roots of Chebyshev polynomials (T,U,V,W)
#   within [-1,1] and shifted within [0,+1]

# zeros of T(x)

chebTzero(n, k) = chebTzero(Float64, n, k)
chebTzero(::Type{T}, n, k) where {T} = cospi(T(n-k + 1/2) / n)

shift_chebTzero(n, k) = shift_chebTzero(Float64, n, k)
shift_chebTzero(::Type{T}, n, k) where {T} = (chebTzero(n,k) + 1) / 2

chebTzeros(n) = chebTzeros(Float64, n)
function chebTzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebTzero(T,n,k)
    end
    return result
end

chebTzeros01(n) = chebTzeros01(Float64, n)
function chebTzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebTzero(T,n,k)
    end
    return result
end

# zeros of U(x)

chebUzero(n, k) = chebUzero(Float64, n, k)
chebUzero(::Type{T}, n, k) where {T} = cospi(T(n-k+1) / T(n+1))

shift_chebUzero(n, k) = shift_chebUzero(Float64, n, k)
shift_chebUzero(::Type{T}, n, k) where {T} = (chebUzero(T,n,k) + 1) / 2

chebUzeros(n) = chebUzeros(Float64, n)
function chebUzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebUzero(T,n,k)
    end
    return result
end

chebUzeros01(n) = chebUzeros01(Float64, n)
function chebUzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebUzero(T,n,k)
    end
    return result
end

# zeros of V(x)

chebVzero(n, k) = chebVzero(Float64, n, k)
chebVzero(::Type{T}, n, k) where {T} = cospi(T(n-k+1/2) / T(n+1/2))

shift_chebVzero(n, k) = shift_chebVzero(Float64, n, k)
shift_chebVzero(::Type{T}, n, k) where {T} = (chebVzero(T,n,k) + 1) / 2

chebVzeros(n) = chebVzeros01(Float64, n)
function chebVzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebVzero(T,n,k)
    end
    return result
end

chebVzeros01(n) = chebVzeros01(Float64, n)
function chebVzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebVzero(T,n,k)
    end
    return result
end

# zeros of W(x)

chebWzero(n, k) = chebWzero(Float64, n, k)
chebWzero(::Type{T}, n, k) where {T} = cospi(T(n-k+1) / T(n+1/2))

shift_chebWzero(n, k) = shift_chebWzero(Float64, n, k)
shift_chebWzero(::Type{T}, n, k) where {T} = (chebWzero(T,n,k) + 1) / 2

chebWzeros(n) = chebWzeros(Float64, n)
function chebWzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebWzero(T,n,k)
    end
    return result
end

chebWzeros01(n) = chebWzeros01(Float64, n)
function chebWzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWzero(T,n,k)
    end
    return result
end

#   extrema of weighted Chebyshev polynomials (T,wtU,wtV,wtW)
#   within [-1,1] and shifted within [0,+1] and
#      unique(sort([extrema_shifted..., (1 .- extrema_shifted)...,]))

# extrema of T(x)

chebTextremum(n, k) = chebTextremum(Float64, n, k)
chebTextremum(::Type{T}, n, k) where {T} = cospi(T(n-k) / n)

shift_chebTextremum(n, k) = shift_chebTextrema(Float64, n, k) 
shift_chebTextremum(::Type{T}, n, k) where {T} = (chebTextremum(T, n, k) + 1)/2

chebTextrema(n) = chebTextrema(Float64, n)
function chebTextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebTextremum(T,n,k)
    end
    return result
end

chebTextrema01(n) = chebTextrema01(Float64, n)
function chebTextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebTextremum(T,n,k)
    end
    return result
end

# extrema of sqrt(1-x^2) * U(x)

chebUextremum(n, k) = chebUextremum(Float64, n, k) 
chebUextremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*(n+1)))

shift_chebUextremum(n, k) = shift_chebUextremum(Float64, n, k) 
shift_chebUextremum(::Type{T}, n, k) where {T} = (chebUextremum(T, n, k) + 1)/2

chebUextrema(n) = chebUextrema(Float64, n)
function chebUextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebUextremum(T,n,k)
    end
    return result
end

chebUextrema01(n) = chebUextrema01(Float64, n)
function chebUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebUextremum(T,n,k)
    end
    return result
end

chebUUextrema01(n, k) = chebUUextrema01(Float64, n, k)
function chebUUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebUextrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1+x) * V(x)

chebVextremum(n, k) = chebVextrema(Float64, n, k) 
chebVextremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)) / T(2*n+1))

shift_chebVextremum(n, k) = shift_chebVextremum(Float64, n, k) 
shift_chebVextremum(::Type{T}, n, k) where {T} = (chebVextremum(T, n, k) + 1)/2

chebVextrema(n) = chebVextrema(Float64, n)
function chebVextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebVextremum(T,n,k)
    end
    return result
end

chebVextrema01(n) = chebVextrema01(Float64, n)
function chebVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebVextremum(T,n,k)
    end
    return result
end

chebVVextrema01(n) = chebVVextrema01(Float64, n)
function chebVVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebVextrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1-x) * W(x)

chebWextremum(n, k) = chebWextrema(Float64, n, k) 
chebWextremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*n+1))

shift_chebWextremum(n, k) = shift_chebWextremum(Float64, n, k) 
shift_chebWextremum(::Type{T}, n, k) where {T} = (chebWextremum(T, n, k) + 1)/2

chebWextrema(n) = chebWextrema(Float64, n)
function chebWextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebWextremum(T,n,k)
    end
    return result
end

chebWextrema01(n) = chebWextrema01(Float64, n)
function chebWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWextremum(T,n,k)
    end
    return result
end

chebWWextrema01(n) = chebWWextrema01(Float64, n)
function chebWWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebWextrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end



end # ChebyshevNodes
