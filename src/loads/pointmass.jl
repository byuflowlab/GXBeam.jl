"""
    PointMass{T}

Type which contains the aggregated inertial properties of one or more point masses which 
are rigidly attached to the center of an element.

# Fields:
 - `mass`: Mass matrix corresponding to the point masses.
"""
struct PointMass{T}
    mass::SMatrix{6,6,T,36}
end
Base.eltype(::PointMass{T}) where T = T

function PointMass{T}(p::PointMass) where T
    PointMass{T}(p.mass)
end

"""
    PointMass(mass)

Define a point mass given its mass matrix
"""
PointMass(mass::AbstractMatrix{T}) where T = PointMass{T}(SMatrix{6,6,T,36}(mass))

"""
    PointMass(m, p, I)

Define a point mass given its mass `m`, offset `p`, and inertia matrix `I`
"""
PointMass(m, p, I) = PointMass([m*I3 -m*tilde(p); m*tilde(p) I-m*tilde(p)*tilde(p)])

"""
    combine_masses(masses)

Combine the point masses in the iterable collection `masses`
"""
function combine_masses(masses)
    M = @SMatrix zeros(6,6)
    for mass in masses
        M += mass.M
    end
    return PointMass(M)
end