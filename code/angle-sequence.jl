# A pure Julia imppementation of the angle sequence algorithm for Laurent polynomials (Polynomials.jl).
# Designed to work with the LaurentPolynomial type from Polynomials.jl which DOES NOT require definite parity.
# Based on the paper https://arxiv.org/pdf/2003.02831
# Adapted from the Python implementation in 
#   https://github.com/alibaba-edu/angle-sequence
# and
#   https://github.com/ichuang/pyqsp


function aligned(p::LaurentPolynomial, lo::Int, hi::Int, nonzero::Bool=false)
    # Returns coefficients of p aligned from lo to hi
    dmin, dmax = [firstindex(p),lastindex(p)]
    @assert lo <= dmin &&  dmax <= hi "check bounds or use trunc()"
    lower = dmin - round(Int, (dmin-lo)/ 2, RoundDown) * 2
    upper = dmax + div(hi-dmax, 2) * 2
    el = eltype(p)
    delta = (nonzero ? 2 : 1)
    return vcat(
        zeros(el, round(Int, (dmin-lo) / 2, RoundDown) * 2),
        p.coeffs,
        zeros(el, div(hi-dmax, 2) * 2)
    )[1:delta:end], lower:delta:upper
end

function trunc(p::LaurentPolynomial, lo::Int, hi::Int, nonzero::Bool=false)
    dmin, dmax = [firstindex(p),lastindex(p)]
    lower = min(lo, dmin)
    upper = max(hi, dmax)
    aligned_coeffs, idx = aligned(p, lower, upper +2, nonzero)
    coeff_dict = Dict(zip(idx, aligned_coeffs))
    return LaurentPolynomial([get(coeff_dict,k, 0.0) for k in lo:hi],lo)
end

function expand_coeffs(coeffs::Vector)
    z = zeros(eltype(coeffs), 2*length(coeffs)-1)
    z[1:2:end] = coeffs
    return z
end

function linear_system(IPoly::LaurentPolynomial, XPoly::LaurentPolynomial, ldeg::Int)
    deg = max(maximum(abs.([firstindex(IPoly),lastindex(IPoly)])),maximum(abs.([firstindex(XPoly),lastindex(XPoly)]))) 

    aligned_icoefs, iidx = aligned(IPoly, -deg, deg, true)
    aligned_xcoefs, xidx = aligned(XPoly, -deg, deg, true)

    function vec_to_mat(vec)
        full_vec = vcat(vec, zeros(ldeg))
        col = vcat([vec[1]], zeros(ldeg))
        return Toeplitz(full_vec, col)
    end

    M = vcat(
        hcat(vec_to_mat(aligned_icoefs), vec_to_mat(-reverse(aligned_xcoefs))),
        hcat(vec_to_mat(aligned_xcoefs), vec_to_mat(reverse(aligned_icoefs)))
    )

    a = ones(ldeg + 1)
    b = zeros(ldeg + 1)

    m = vcat(
        permutedims(vcat(a, b)),
        permutedims(vcat(b, a)),
        M[1:ldeg, :],
        M[deg+2 : deg + 2 * ldeg + 1, :], 
        M[end-ldeg+1:end, :]
    )
    s = zeros(size(m, 1))
    s[1] = 1.0

    return m, s, deg
end

function decompose(IPoly::LaurentPolynomial, XPoly::LaurentPolynomial, ldeg::Int)
    m, s, deg = linear_system(IPoly, XPoly, ldeg)
    v = m \ s
    
    li = expand_coeffs(v[1:ldeg+1])
    lx = expand_coeffs(v[end-ldeg:end])

    l_IPoly = LaurentPolynomial(li, -ldeg)
    l_XPoly = LaurentPolynomial(lx, -ldeg)

    # g is represented as (IPoly, XPoly), l as (l_IPoly, l_XPoly)
    # We compute r = l * g, then truncate
    l = [l_IPoly, l_XPoly]
    #g = [IPoly, XPoly]
    inv_IPoly = LaurentPolynomial(reverse(conj(IPoly.coeffs)), -lastindex(IPoly))
    inv_XPoly = LaurentPolynomial(reverse(conj(XPoly.coeffs)), -lastindex(XPoly))
    
    r = [(l_IPoly * IPoly - l_XPoly * inv_XPoly), 
         (l_IPoly * XPoly + l_XPoly * inv_IPoly)]
    r = [trunc(r[1], -deg + ldeg, deg - ldeg), trunc(r[2], -deg + ldeg, deg-ldeg)]
    r = [LaurentPolynomial(r[1], -deg + ldeg), LaurentPolynomial(r[2], -deg + ldeg)]

    #inv_l = ~l = [~l_IPoly, -l_XPoly]
    inv_l = [LaurentPolynomial(reverse(conj(li)), -ldeg), -l_XPoly]
    return inv_l, r
end

function left_and_right_angles(IPoly::LaurentPolynomial, XPoly::LaurentPolynomial)
    deg = max(maximum(abs.([firstindex(IPoly),lastindex(IPoly)])),maximum(abs.([firstindex(XPoly),lastindex(XPoly)]))) 
    @assert deg == 1
    x = IPoly(1) + im*XPoly(1)
    j = eltype(IPoly) <: Complex ? big(1.0im) : 1.0im
    y = IPoly(j) - im*XPoly(j)
    s = atan(x.im, x.re)
    d = atan(y.im, y.re) - π/2
    res = [(s+d)/2, (s-d)/2]
    return res
end

function angseq(IPoly::LaurentPolynomial, XPoly::LaurentPolynomial)
    deg = max(maximum(abs.([firstindex(IPoly),lastindex(IPoly)])),maximum(abs.([firstindex(XPoly),lastindex(XPoly)]))) 
    if deg == 1
        return left_and_right_angles(IPoly, XPoly)
    else
        l, r = decompose(IPoly, XPoly, div(deg, 2)) 
        a = angseq(l...)
        b = angseq(r...)
        return vcat(a[1:end-1], [a[end] + b[1]], b[2:end])
    end
end