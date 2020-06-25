struct PrescribedCondition{T}
    pt::Int
    val::NTuple{6, T}
    force::NTuple{6, Bool}
    follower::NTuple{6, Bool}
end
