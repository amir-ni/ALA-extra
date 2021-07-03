using LinearAlgebra
using Base: require_one_based_indexing


"""
    systemSolver(A::AbstractMatrix, b::Vector, x_initial::Vector) -> Vector

    Compute the orthogonal projection of a given point to the 
solution set of a system of linear equations. Also a direct 
method for solving systems of linear equations. The output 
is either the projection or inconsistency of the system.

# Arguments
 - `A::AbstractMatrix{Real}`: coefficients of linear equations
 - `b::AbstractVector{Real}`: intended results of equations
 - `x0::AbstractVector{Real}`: initial point for projection
 - `β::Real`: set in order to obtain comparable summands in norm

# References
> Ján Plesník, Finding the orthogonal projection of a point onto 
> an affine subspace, Linear Algebra and its Applications, Volume 422, 
> Issues 2–3, 2007, Pages 455-470
"""
function systemSolver(A::AbstractMatrix{<:Real},    # in: equations matrix (m x n) rows = coefficients
                    b::AbstractVector{<:Real},      # in: equations results
                    x0::AbstractVector{<:Real},     # in: initial point
                    β::Real=1.0)                    # in: optional normalizer
    
    require_one_based_indexing(A, b, x0)
    m, n = size(A)
    if m != length(b) || n != length(x0)
        throw(DimensionMismatch("dimensions of `A`, $m * $n, does not match one of the lengths of `x0`, $(length(x0)), or `b`, $(length(b))"))
    end
    for i in 1:m
        iszero(A[i, :]) && throw(ArgumentError("`A` rows must be nonzero"))
    end
    iszero(β) && throw(ArgumentError("`β` must be nonzero"))

    if m == 0
        return x0
    end

    k = 1
    p = 1
    V = hcat(A[1, :])
    x = x0 + (b[1] - dot(A[1, :], x0))/(dot(A[1, :], V[:,1])) * V[:,1]
    
    while k != m
        k = k + 1
        p = p + 1
        if dot((A[k, :] - A[1, :]),V[:, k-1]) == 0
            A[k, :] = - A[k, :]
            b[k] = - b[k]
        end

        w = x0 + β * A[k, :]
        for j in 2:k
            w = w + (b[j]-b[1]-dot(A[j, :]-A[1, :],w))/dot(A[j, :]-A[1, :],V[:,j-1]) * V[:,j-1]
        end

        u = x + (b[k] - dot(A[k, :], x))/dot(A[k, :]-A[1, :], V[:,k-1]) * V[:,k-1]
        if u != w
            V = hcat(V, w - u)
            x = u + (b[1] - dot(A[1, :], u))/dot(A[1, :], V[:, k]) * V[:, k]
        elseif dot(A[k, :],x) == b[k]
            if k != m
                A = A[setdiff(1:end, k), :]
                b = b[setdiff(1:end, k)]
                k = k - 1
                m = m - 1
            end
        else
            throw(DomainError("the ($p)th equation is inconsistent"))
        end
    end
    return x
end