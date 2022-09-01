ChainRulesCore.rrule(T::Type{<:SArray}, xs::Number...) = T(xs...), dv->(NoTangent(), dv...)
ChainRulesCore.rrule(T::Type{<:SArray}, x::Tuple) = T(x), dv->(NoTangent(), tuple(dv...))
ChainRulesCore.rrule(T::Type{<:SArray}, x::AbstractArray) = T(x), dv->(NoTangent(), dv)