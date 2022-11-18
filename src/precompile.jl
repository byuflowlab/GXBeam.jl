
let 
    
    # set default input types
    jacob=SparseMatrixCSC{Float64, Int}
    resid=Vector{Float64}
    x=Vector{Float64}
    indices=SystemIndices 
    two_dimensional=Bool 
    force_scaling=Bool
    assembly=Assembly{Float64, Vector{SVector{3,Float64}}, Vector{Int}, Vector{Element{Float64}}}
    pcond=Dict{Int, PrescribedConditions{Float64}}
    dload=Dict{Int, DistributedLoads{Float64}}
    pmass=Dict{Int, PointMass{Float64}}
    gvec=SVector{3,Float64}

    ub_p = SVector{3,Float64} 
    θb_p = SVector{3,Float64} 
    vb_p = SVector{3,Float64} 
    ωb_p = SVector{3,Float64} 
    ab_p = SVector{3,Float64} 
    αb_p = SVector{3,Float64}

    V0 = Vector{SVector{3,Float64}}
    Omega0 = Vector{SVector{3,Float64}}
    Vdot0 = Vector{SVector{3,Float64}}
    Omegadot0 = Vector{SVector{3,Float64}}

    ubdot = SVector{3,Float64}
    θbdot = SVector{3,Float64}
    vbdot = SVector{3,Float64}
    ωbdot = SVector{3,Float64}

    udot = Vector{SVector{3,Float64}}
    θdot = Vector{SVector{3,Float64}}
    Vdot = Vector{SVector{3,Float64}}
    Ωdot = Vector{SVector{3,Float64}}

    dt = Float64

    # precompile static methods

    precompile(static_system_residual!, (resid, x, indices, two_dimensional, force_scaling, 
        assembly, pcond, dload, pmass, gvec))

    precompile(static_system_jacobian!, (jacob, x, indices, two_dimensional, force_scaling, 
        assembly, pcond, dload, pmass, gvec))

    # precompile steady state methods

    precompile(steady_system_residual!(similar(x), x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    precompile(steady_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    # precompile initial condition methods

    precompile(initial_system_residual!(similar(x), x, indices, rate_vars1, 
        rate_vars2, two_dimensional, force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p,
        u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    precompile(initial_system_jacobian!(J, x, indices, rate_vars1, rate_vars2, two_dimensional,
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p,
        u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    # precompile newmark time integration methods

    precompile(newmark_system_residual!(similar(x), x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

    precompile(newmark_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

    # precompile dynamic methods

    precompile(dynamic_system_residual!(similar(x), dx, x, indices, two_dimensional,
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    dynamic_system_jacobian!(J, dx, x, indices, force_scaling, two_dimensional, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    precompile(system_mass_matrix!(M, x, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    # precompile steady constant mass matrix methods

    precompile(expanded_steady_system_residual!(similar(x), x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    precompile(expanded_steady_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    # precompile dynamic constant mass matrix methods

    precompile(expanded_dynamic_system_residual!(similar(dx), dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    precompile(expanded_dynamic_system_jacobian!(J, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    precompile(expanded_system_mass_matrix!(M, indices, two_dimensional, force_scaling, assembly, 
        pcond, pmass)

end