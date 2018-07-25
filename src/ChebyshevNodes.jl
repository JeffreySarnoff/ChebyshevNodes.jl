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

chebTzeros(n, k) = chebTzeros(Float64, n, k)
function chebTzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebTzero(T,n,k)
    end
    return result
end

chebTzeros01(n, k) = chebTzeros01(Float64, n, k)
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

chebUzeros(n, k) = chebUzeros(Float64, n, k)
function chebUzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebUzero(T,n,k)
    end
    return result
end

chebUzeros01(n, k) = chebUzeros01(Float64, n, k)
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

chebVzeros(n, k) = chebVzeros01(Float64, n, k)
function chebVzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebVzero(T,n,k)
    end
    return result
end

chebVzeros01(n, k) = chebVzeros01(Float64, n, k)
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

chebWzeros(n, k) = chebWzeros(Float64, n, k)
function chebWzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebWzero(T,n,k)
    end
    return result
end

chebWzeros01(n, k) = chebWzeros01(Float64, n, k)
function chebWzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWzero(T,n,k)
    end
    return result
end

# extrema of _weighted_ Chebyshev polynomials (types T,U,V,W or 1,2,3,4)
# includes reverse order 1-extrema with extrema 

#   extrema of weighted Chebyshev polynomials (T,wtU,wtV,wtW)
#   within [-1,1] and shifted within [0,+1] and
#      unique(sort([extrema_shifted..., (1 .- extrema_shifted)...,]))


# extrema of T(x)

chebTextrema(n, k) = chebTextrema(Float64, n, k)
chebTextrema(::Type{T}, n, k) where {T} = cospi(T(n-k) / n)

shift_chebTextrema(n, k) = shift_chebTextrema(Float64, n, k) 
shift_chebTextrema(::Type{T}, n, k) where {T} = (chebTextrema(T, n, k) + 1)/2

chebTextrema(n, k) = chebTextrema(Float64, n, k)
function chebTextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebTextrema(T,n,k)
    end
    return result
end

chebTextrema01(n, k) = chebTextrema01(Float64, n, k)
function chebTextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebTextrema(T,n,k)
    end
    return result
end



# extrema of sqrt(1-x^2) * U(x)

chebUextrema(n, k) = chebUextrema(Float64, n, k) 
chebUextrema(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*(n+1)))

shift_chebUextrema(n, k) = shift_chebUextrema(Float64, n, k) 
shift_chebUextrema(::Type{T}, n, k) where {T} = (chebUextrema(T, n, k) + 1)/2

chebUextrema(n, k) = chebUextrema(Float64, n, k)
function chebUextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebUextrema(T,n,k)
    end
    return result
end

chebUextrema01(n, k) = chebUextrema01(Float64, n, k)
function chebUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebUextrema(T,n,k)
    end
    return result
end

chebUUextrema01(n, k) = chebUUextrema01(Float64, n, k)
function chebUUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebUextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end


# extrema of sqrt(1+x) * V(x)

chebVextrema(n, k) = chebVextrema(Float64, n, k) 
chebVextrema(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)) / T(2*n+1))

shift_chebVextrema(n, k) = shift_chebVextrema(Float64, n, k) 
shift_chebVextrema(::Type{T}, n, k) where {T} = (chebVextrema(T, n, k) + 1)/2

chebVextrema(n, k) = chebVextrema(Float64, n, k)
function chebVextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebVextrema(T,n,k)
    end
    return result
end

chebVextrema01(n, k) = chebVextrema01(Float64, n, k)
function chebVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebVextrema(T,n,k)
    end
    return result
end

chebVVextrema01(n, k) = chebVVextrema01(Float64, n, k)
function chebVVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebVextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1-x) * W(x)

chebWextrema(n, k) = chebWextrema(Float64, n, k) 
chebWextrema(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*n+1))

shift_chebWextrema(n, k) = shift_chebWextrema(Float64, n, k) 
shift_chebWextrema(::Type{T}, n, k) where {T} = (chebWextrema(T, n, k) + 1)/2

chebWextrema(n, k) = chebWextrema(Float64, n, k)
function chebWextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebWextrema(T,n,k)
    end
    return result
end

chebWextrema01(n, k) = chebWextrema01(Float64, n, k)
function chebWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWextrema(T,n,k)
    end
    return result
end

chebWWextrema01(n, k) = chebWWextrema01(Float64, n, k)
function chebWWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebWextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end



end # ChebyshevNodes
