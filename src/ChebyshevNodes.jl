module ChebyshevNodes

export chebt_zeros, chebu_zeros, chebv_zeros, chebw_zeros,
       chebt_zeros01, chebu_zeros01, chebv_zeros01, chebw_zeros01,
       chebtt_zeros01, chebuu_zeros01, chebvv_zeros01, chebww_zeros01,
       chebt_extrema, chebu_extrema, chebv_extrema, chebw_extrema,
       chebt_extrema01, chebu_extrema01, chebv_extrema01, chebw_extrema01
       chebtt_extrema, chebuu_extrema, chebvv_extrema, chebww_extrema,
       chebtt_extrema01, chebuu_extrema01, chebvv_extrema01, chebww_extrema01,
       fromab, intoab


# mapping from [a, b] to [-1, 1]
@inline fromab(a,b,x) =  (2*x - (a+b))/(b - a)
# mapping from [-1, 1] to [a, b]
@inline intoab(a,b,y) = (y * (b - a) + (a + b))/2
# mapping from [0, 1] to [-1, 1]
@inline from01(x) = 2*x - 1
# mapping from [-1, 1] to [0, 1]
@inline into01(y) = (y + 1)/2


# k=1..n

# -1..1
T0(k,n)  = cospi((2(n-k+1)-1)/2n)
Tᴱ0(k,n) = cospi((2(k-n-1)+1)/2n)/cospi(inv(2n))   # extended to -1..1
Tᴬ0(k,n) = k==1 ? -1 : (k==n ? 1 : T0(k-1,n-2))   # augmented with -1,+1

# shifted 0..1
𝑇0(k,n)  = into01(T0(k,n))
𝑇ᴱ0(k,n) = into01(Tᴱ0(k,n))   # extended to -1..1
𝑇ᴬ0(k,n) = into01(Tᴬ0(k,n))   # augmented with -1,+1

 function chebT(n,x)
     if -1 <= x <= 1
         Θ = acos(x)
         Θ *= n
         cos(Θ)
     elseif 1 < x <= realmax(x)
         Θ = acosh(x)
         Θ *= n
         cosh(Θ)
     elseif -realmax(x) <= x < -1
         Θ = acosh(-x)
         Θ *= n
         iseven(n) ? cosh(Θ) : -cosh(Θ)
     else
         throw(DomainError("$x is outside of the domain"))
     end
end

# -1..1
U0(k,n)  = cospi((n-k+1)/(n+1))
Uᴱ0(k,n) = cospi((n-k+1)/(n+1))/cospi(1/(n+1))    # extended to -1..1
Uᴬ0(k,n) = k==1 ? -1 : (k==n ? 1 : U0(k-1,n-2))   # augmented with -1,+1

# shifted 0..1
𝑈0(k,n)  = into01(𝑈0(k,n))
𝑈ᴱ0(k,n) = into01(𝑈ᴱ0(k,n))   # extended to -1..1
𝑈ᴬ0(k,n) = into01(𝑈ᴬ0(k,n))   # augmented with -1,+1

TGLextrema(k,n) = -cospi((k-1)/(n-1)) # Gauss-Lobatto Chebyshev nodes


#   roots of Chebyshev polynomials (T,U,V,W)
#   within [-1,1] and shifted within [0,+1]

# zeros of T(x)

chebt_zero(n, k) = chebt_zero(Float64, n, k)
chebt_zero(::Type{T}, n, k) where {T} = cospi(T(n-k + 1/2) / n)

chebtstar_zero(n, k) = chebtstar_zero(Float64, n, k)
chebtstar_zero(::Type{T}, n, k) where {T} = (chebt_zero(n,k) + 1) / 2

chebt_zeros(n) = chebt_zeros(Float64, n)
function chebt_zeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebt_zero(T,n,k)
    end
    return result
end

chebt_zeros01(n) = chebt_zeros01(Float64, n)
function chebt_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebtstar_zero(T,n,k)
    end
    return result
end

chebtt_zeros01(n) = chebtt_zeros01(Float64, n)
function chebtt_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebt_zeros(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# zeros of U(x)

chebu_zero(n, k) = chebu_zero(Float64, n, k)
chebu_zero(::Type{T}, n, k) where {T} = cospi(T(n-k+1) / T(n+1))

chebustar_zero(n, k) = chebustar_zero(Float64, n, k)
chebustar_zero(::Type{T}, n, k) where {T} = (chebu_zero(T,n,k) + 1) / 2

chebu_zeros(n) = chebu_zeros(Float64, n)
function chebu_zeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebustar_zero(T,n,k)
    end
    return result
end

chebu_zeros01(n) = chebu_zeros01(Float64, n)
function chebu_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebustar_zero(T,n,k)
    end
    return result
end

chebuu_zeros01(n) = chebuu_zeros01(Float64, n)
function chebuu_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebu_zeros(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end


# zeros of V(x)

chebv_zero(n, k) = chebv_zero(Float64, n, k)
chebv_zero(::Type{T}, n, k) where {T} = cospi(T(n-k+1/2) / T(n+1/2))

chebvstar_zero(n, k) = chebvstar_zero(Float64, n, k)
chebvstar_zero(::Type{T}, n, k) where {T} = (chebv_zero(T,n,k) + 1) / 2

chebv_zeros(n) = chebv_zeros01(Float64, n)
function chebv_zeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebvstar_zero(T,n,k)
    end
    return result
end

chebv_zeros01(n) = chebv_zeros01(Float64, n)
function chebv_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebvstar_zero(T,n,k)
    end
    return result
end

chebvv_zeros01(n) = chebvv_zeros01(Float64, n)
function chebv_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebv_zeros(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# zeros of W(x)

chebw_zero(n, k) = chebw_zero(Float64, n, k)
chebw_zero(::Type{T}, n, k) where {T} = cospi(T(n-k+1) / T(n+1/2))

chebwstar_zero(n, k) = chebwstar_zero(Float64, n, k)
chebwstar_zero(::Type{T}, n, k) where {T} = (chebw_zero(T,n,k) + 1) / 2

chebw_zeros(n) = chebw_zeros(Float64, n)
function chebw_zeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = chebwstar_zero(T,n,k)
    end
    return result
end

chebw_zeros01(n) = chebw_zeros01(Float64, n)
function chebw_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebwstar_zero(T,n,k)
    end
    return result
end

chebww_zeros01(n) = chebww_zeros01(Float64, n)
function chebww_zeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebw_zeros(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

#   extrema of weighted Chebyshev polynomials (T,wtU,wtV,wtW)
#   within [-1,1] and shifted within [0,+1] and
#      unique(sort([extrema_shifted..., (1 .- extrema_shifted)...,]))

# extrema of T(x)

chebt_extremum(n, k) = chebt_extremum(Float64, n, k)
chebt_extremum(::Type{T}, n, k) where {T} = cospi(T(n-k) / n)

chebtstar_extremum(n, k) = shift_chebt_extrema(Float64, n, k)
chebtstar_extremum(::Type{T}, n, k) where {T} = (chebt_extremum(T, n, k) + 1)/2

chebt_extrema(n) = chebt_extrema(Float64, n)
function chebt_extrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebt_extremum(T,n,k)
    end
    return result
end

chebt_extrema01(n) = chebt_extrema01(Float64, n)
function chebt_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebtstar_extremum(T,n,k)
    end
    return result
end

chebtt_extrema01(n) = chebtt_extrema01(Float64, n)
function chebtt_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebt_extrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1-x^2) * U(x)

chebu_extremum(n, k) = chebu_extremum(Float64, n, k)
chebu_extremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*(n+1)))

chebustar_extremum(n, k) = chebustar_extremum(Float64, n, k)
chebustar_extremum(::Type{T}, n, k) where {T} = (chebu_extremum(T, n, k) + 1)/2

chebu_extrema(n) = chebu_extrema(Float64, n)
function chebu_extrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 0:n
        result[k] = chebu_extremum(T,n,k)
    end
    return result
end

chebu_extrema01(n) = chebu_extrema01(Float64, n)
function chebu_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebustar_extremum(T,n,k)
    end
    return result
end

chebuu_extrema01(n) = chebuu_extrema01(Float64, n)
function chebuu_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebu_extrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1+x) * V(x)

chebv_extremum(n, k) = chebv_extrema(Float64, n, k)
chebv_extremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)) / T(2*n+1))

chebvstar_extremum(n, k) = chebvstar_extremum(Float64, n, k)
chebvstar_extremum(::Type{T}, n, k) where {T} = (chebv_extremum(T, n, k) + 1)/2

chebv_extrema(n) = chebv_extrema(Float64, n)
function chebv_extrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebv_extremum(T,n,k)
    end
    return result
end

chebv_extrema01(n) = chebv_extrema01(Float64, n)
function chebv_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebvstar_extremum(T,n,k)
    end
    return result
end

chebvv_extrema01(n) = chebvv_extrema01(Float64, n)
function chebvv_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebv_extrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

# extrema of sqrt(1-x) * W(x)

chebw_extremum(n, k) = chebw_extrema(Float64, n, k)
chebw_extremum(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*n+1))

chebwstar_extremum(n, k) = chebwstar_extremum(Float64, n, k)
chebwstar_extremum(::Type{T}, n, k) where {T} = (chebw_extremum(T, n, k) + 1)/2

chebw_extrema(n) = chebw_extrema(Float64, n)
function chebw_extrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebw_extremum(T,n,k)
    end
    return result
end

chebw_extrema01(n) = chebw_extrema01(Float64, n)
function chebw_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = chebwstar_extremum(T,n,k)
    end
    return result
end

chebww_extrema01(n) = chebww_extrema01(Float64, n)
function chebww_extrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebw_extrema(T, n)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    unique!(sort!(result))
    return result
end

end # ChebyshevNodes
