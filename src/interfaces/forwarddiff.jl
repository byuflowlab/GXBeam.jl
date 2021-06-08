struct ForwardDiffFactorization{T,D,N,M,F} <: Factorization{D}
    value::M
    partials::NTuple{N,M}
    factorization::F
end

ForwardDiffFactorization(A) = ForwardDiffFactorization(factorize, A)

function ForwardDiffFactorization(f, A::AbstractMatrix{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    Av = ForwardDiff.value.(A)
    Ap = ntuple(i->ForwardDiff.partials.(A, i), N)
    Afac = f(Av)
    return ForwardDiffFactorization{T,ForwardDiff.Dual{T,V,N}}(Av, Ap, Afac)
end

function ForwardDiffFactorization{T,D}(Av::M, Ap::NTuple{N,M}, Afac::F) where {T,D,N,M,F}

    return ForwardDiffFactorization{T,D,N,M,F}(Av, Ap, Afac)
end

function Base.adjoint(A::ForwardDiffFactorization{T,D,N,M,F}) where {T,D,N,M,F}
    # value, partials, and factorization of A
    Av = adjoint(A.value)
    Ap = adjoint.(A.partials)
    Afac = adjoint(A.factorization)
    return ForwardDiffFactorization{T,D,N,M,F}(Av, Ap, Afac)
end

function Base.transpose(A::ForwardDiffFactorization{T,D,N,M,F}) where {T,D,N,M,F}
    # value, partials, and factorization of A
    Av = transpose(A.value)
    Ap = transpose.(A.partials)
    Afac = transpose(A.factorization)
    return ForwardDiffFactorization{T,D,N,M,F}(Av, Ap, Afac)
end

function LinearAlgebra.ldiv!(A::ForwardDiffFactorization{T,D,N,M,F}, B) where {T,D,N,M,F}

    # value, partials, and factorization of A
    Av = A.value
    Ap = A.partials
    Afac = A.factorization

    # values and partials of B
    Bv = ForwardDiff.value.(B)
    Bp = ntuple(i->ForwardDiff.partials.(B, i), N)

    # calculate output values, store in Bv
    ldiv!(Afac, Bv)

    # calculate output partials, store in Bp
    ldiv!.(Ref(Afac), mul!.(Bp, Ap, Ref(Bv), -1, 1))

    # overwrite B with result
    for j in eachindex(B)
        B[j] = ForwardDiff.Dual{T}(Bv[j], ntuple(i->Bp[i][j], N)...)
    end

    return B
end

safe_lu(A; check::Bool = true) = lu(A; check)

function safe_lu(A::SparseMatrixCSC{<:ForwardDiff.Dual, <:Any}; check::Bool = true)
    return ForwardDiffFactorization(lu, A)
end
