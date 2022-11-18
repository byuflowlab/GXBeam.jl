# For StaticArrays
ChainRulesCore.rrule(T::Type{<:SArray}, xs::Number...) = T(xs...), dv->(NoTangent(), dv...)
ChainRulesCore.rrule(T::Type{<:SArray}, x::Tuple) = T(x), dv->(NoTangent(), tuple(dv...))
ChainRulesCore.rrule(T::Type{<:SArray}, x::AbstractArray) = T(x), dv->(NoTangent(), dv)

# chain rule definitions for `_rotation_parameter_scaling(θ1, θ2, θ3)`

# function ChainRulesCore.frule((NoTangent(), Δθ1, Δθ2, Δθ3),
#     ::typeof(_rotation_parameter_scaling), θ1::Number, θ2::Number, θ3::Number)

#     println("HERE: frule")
#     scaling = _rotation_parameter_scaling(θ1, θ2, θ3)
#     scaling_θ = ifelse(isone(scaling), zero(θ), (1 - scaling)*θ/(θ'*θ))

#     return scaling, ((scaling_θ[1] * Δθ1 + scaling_θ[2] * Δθ2 + scaling_θ[3] * Δθ3),)
# end

# function ChainRulesCore.rrule(::typeof(_rotation_parameter_scaling), θ1, θ2, θ3)

#     println("HERE: rrule")
#     θ = SVector(θ1, θ2, θ3)
#     scaling = rotation_parameter_scaling(θ)
#     scaling_θ = ifelse(isone(scaling), zero(θ), (1.0 - scaling)*θ/(θ'*θ))

#     return scaling, dv -> (NoTangent(), scaling_θ[1] * dv, scaling_θ[2] * dv, scaling_θ[3] * dv)
# end