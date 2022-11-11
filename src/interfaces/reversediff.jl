# NOTE: The following code engages in type piracy and therefore should be relocated

# overloads to resolve method ambiguities
Base.:+(A::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M}, B::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}) = Array(A) + B
Base.:-(A::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M}, B::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}) = Array(A) - B
Base.:*(A::StaticArraysCore.StaticArray{Tuple{N, M}, T, 2} where {N, M, T}, B::ReverseDiff.TrackedArray{V, D, 1} where {V, D}) = Array(A)*B
Base.:*(A::StaticArraysCore.StaticArray{Tuple{N, M}, T, 2} where {N, M, T}, B::ReverseDiff.TrackedArray{V, D, 2} where {V, D}) = Array(A)*B

Base.:+(A::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}, B::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M})  = A + Array(B)
Base.:-(A::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}, B::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M}) = A - Array(B)
Base.:*(A::ReverseDiff.TrackedArray{V, D, 2} where {V, D}, B::StaticArraysCore.StaticArray{S, T, 1} where {S<:Tuple, T}) = A*Array(B)
Base.:*(A::ReverseDiff.TrackedArray{V, D, 2} where {V, D}, B::StaticArraysCore.StaticArray{S, T, 2} where {S<:Tuple, T}) = A*Array(B)

# handle promotion for NLsolve.SolverResults
# function NLsolve.SolverResults(method::String, initial_x::I, zero::Z, residual_norm::rT, 
#     iterations::Int64, x_converged::Bool, xtol::rT, f_converged::Bool, ftol::rT, 
#     trace::NLsolve.SolverTrace, f_calls::Int64, g_calls::Int64
#     ) where {rT<:ReverseDiff.TrackedReal, T<:Union{rT, Complex{rT}}, I<:(AbstractArray{T}), Z<:(AbstractArray{T})}
    
#     # real type
#     new_rT = promote_type(real(eltype(initial_x)), real(eltype(zero)), typeof(residual_norm), typeof(xtol), typeof(ftol))
    
#     # real/complex type
#     if promote_type(eltype(initial_x), eltype(zero)) <: Complex
#         new_T = Complex{new_rT}
#     else
#         new_T = new_rT
#     end

#     # correct guess element type (if necessary)
#     if !(eltype(initial_x) <: new_T)
#         initial_x = new_T.(initial_x)
#     end

#     # guess type
#     new_I = typeof(initial_x)

#     # correct zero element type (if necessary)
#     if !(eltype(zero) <: new_T)
#         zero = new_T.(zero)
#     end

#     # zero type
#     new_Z = typeof(initial_x)

#     return NLsolve.SolverResults{new_rT, new_T, new_I, new_Z}(method, initial_x, zero, residual_norm, 
#         iterations, x_converged, xtol, f_converged, ftol, trace, f_calls, g_calls)
# end

# Base.:*(A::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M}, B::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}) = Array(A)*B
# Base.:*(A::ReverseDiff.TrackedArray{V, D, N} where {V, D, N}, B::StaticArraysCore.StaticArray{S, T, M} where {S<:Tuple, T, M}) = A*Array(B)
