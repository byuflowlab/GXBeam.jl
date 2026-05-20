"""
    examples_benchmark.jl

Benchmark suite covering the 11 workloads under test/examples/. Each example's
assembly setup is constructed once outside `@benchmark`; only the solve call(s)
are timed. Loops with many parameter steps (load sweeps, rpm sweeps, time-domain
runs) are trimmed where necessary to keep individual workloads in the ~30-120s
budget — see the inventory in DI_MIGRATION_PLAN.md for budgets.

Some examples are split into multiple benchmark entries when they exercise
distinctly different solve paths (e.g. linear/nonlinear/follower for
static-joined-wing).

Run from the GXBeam root directory:
    julia --project=test test/examples_benchmark.jl
"""

using GXBeam, LinearAlgebra, BenchmarkTools, Statistics, Printf, Elliptic

# ============================================================================
# Helpers
# ============================================================================

ms(t) = t / 1e6

function pretty_table(rows)
    println("| Workload | Min (ms) | Mean (ms) | Std (ms) | Med (ms) | Max (ms) | Allocs | Mem (MB) |")
    println("|---|---|---|---|---|---|---|---|")
    for (name, b) in rows
        ts = b.times
        @printf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %d | %.2f |\n",
            name, ms(minimum(ts)), ms(mean(ts)), ms(std(ts)), ms(median(ts)), ms(maximum(ts)),
            b.allocs, b.memory/1e6)
    end
end

# ============================================================================
# 1. cantilever — linear static_analysis
# ============================================================================

function setup_cantilever()
    nelem = 12
    a = 0.3; b = 0.7; L = 1.0
    n1 = n3 = div(nelem, 3); n2 = nelem - n1 - n3
    x1 = range(0, a, length=n1+1); x2 = range(a, b, length=n2+1); x3 = range(b, L, length=n3+1)
    x = vcat(x1, x2[2:end], x3[2:end])
    y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]
    start = 1:nelem; stop = 2:nelem+1
    EI = 1e9
    stiffness = fill(Diagonal([0.0, 0, 0, 0, EI, 0]), nelem)
    assembly = Assembly(points, start, stop, stiffness=stiffness)
    prescribed_conditions = Dict(
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    q = 1000.0
    distributed_loads = Dict()
    for ielem in n1+1:n1+n2
        distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz=(s)->q)
    end
    return (; assembly, prescribed_conditions, distributed_loads)
end

function solve_cantilever(setup)
    static_analysis(setup.assembly,
        prescribed_conditions=setup.prescribed_conditions,
        distributed_loads=setup.distributed_loads, linear=true)
end

# ============================================================================
# 2. curved — nonlinear static_analysis
# ============================================================================

function setup_curved()
    R = 100; L = R*pi/4; h = w = 1
    E = 1e7; ν = 0.0; G = E/(2*(1+ν))
    r = [0.0, 0, 0]
    frame = [0 -1 0; 1 0 0; 0 0 1]
    curvature = [0.0, 0, -1/R]
    A = h*w; Ay = A; Az = A
    Iyy = w*h^3/12; Izz = w^3*h/12; J = Iyy + Izz
    nelem = 16
    ΔL, xp, xm, Cab = discretize_beam(L, r, nelem; frame=frame, curvature=curvature)
    P = 600.0
    start = 1:nelem; stop = 2:nelem+1
    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)
    assembly = Assembly(xp, start, stop, compliance=compliance, frames=Cab, lengths=ΔL, midpoints=xm)
    prescribed_conditions = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        nelem+1 => PrescribedConditions(Fz=P))
    return (; assembly, prescribed_conditions)
end

function solve_curved(setup)
    static_analysis(setup.assembly, prescribed_conditions=setup.prescribed_conditions)
end

# ============================================================================
# 3. overdetermined — linear static_analysis
# ============================================================================

function setup_overdetermined()
    nelem = 16; L = 1
    x = range(0, L, length=nelem+1); y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]
    start = 1:nelem; stop = 2:nelem+1
    EI = 1e7
    compliance = fill(Diagonal([0.0, 0, 0, 0, 1/EI, 0]), nelem)
    assembly = Assembly(points, start, stop, compliance=compliance)
    prescribed_conditions = Dict(
        1 => PrescribedConditions(uz=0),
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    qmax = 1000.0
    distributed_loads = Dict()
    for i = 1:nelem
        distributed_loads[i] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fz=(s)->qmax*s)
    end
    return (; assembly, prescribed_conditions, distributed_loads)
end

function solve_overdetermined(setup)
    static_analysis(setup.assembly,
        prescribed_conditions=setup.prescribed_conditions,
        distributed_loads=setup.distributed_loads, linear=true)
end

# ============================================================================
# 4. tipforce — nonlinear static_analysis! in a load-step loop
# ============================================================================

function setup_tipforce()
    L = 1; EI = 1e6
    λ = 0:0.5:16
    p = EI/L^2
    P = λ*p
    nelem = 16
    x = range(0, L, length=nelem+1); y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:nelem+1]
    start = 1:nelem; stop = 2:nelem+1
    compliance = fill(Diagonal([0.0, 0, 0, 0, 1/EI, 0]), nelem)
    assembly = Assembly(points, start, stop, compliance=compliance)
    system = StaticSystem(assembly)
    pcs = [
        Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            nelem+1 => PrescribedConditions(Fz = P[i])
        ) for i = 1:length(P)]
    return (; assembly, system, pcs)
end

function solve_tipforce(setup)
    reset_state!(setup.system)
    for i = 1:length(setup.pcs)
        static_analysis!(setup.system, setup.assembly, prescribed_conditions=setup.pcs[i])
    end
end

# ============================================================================
# 5. tipmoment — nonlinear static_analysis! in a moment-step loop
# ============================================================================

function setup_tipmoment()
    L = 12; h = w = 1; E = 30e6
    A = h*w; Iyy = w*h^3/12; Izz = w^3*h/12
    λ = [0.0, 0.4, 0.8, 1.2, 1.6, 1.8, 2.0]
    m = pi*E*Iyy/L
    M = λ*m
    nelem = 16
    x = range(0, L, length=nelem+1); y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]
    start = 1:nelem; stop = 2:nelem+1
    compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)
    assembly = Assembly(points, start, stop, compliance=compliance)
    system = StaticSystem(assembly)
    pcs = [
        Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            nelem+1 => PrescribedConditions(Mz = M[i])
        ) for i = 1:length(M)]
    return (; assembly, system, pcs)
end

function solve_tipmoment(setup)
    reset_state!(setup.system)
    for i = 1:length(setup.pcs)
        static_analysis!(setup.system, setup.assembly, prescribed_conditions=setup.pcs[i])
    end
end

# ============================================================================
# 6. static-joined-wing — split into linear / nonlinear / nonlinear-follower
#
# Original uses 141 load steps; trimmed to 21 here for benchmark feasibility
# (consistent across all phases).
# ============================================================================

function setup_static_joined_wing()
    # Same geometry as test/examples/static-joined-wing.jl
    p1 = [-7.1726, -12, -3.21539]; p2 = [-5.37945, -9, -2.41154]
    p3 = [-3.5863, -6, -1.6077];   p4 = [-1.79315, -3, -0.803848]
    p5 = [0.0, 0, 0];              p6 = [7.1726, -12, 3.21539]

    tmp1 = sqrt(p1[1]^2 + p1[2]^2)
    c1, s1 = -p1[1]/tmp1, -p1[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]
    tmp2 = sqrt(p1[1]^2 + p1[2]^2 + p1[3]^2)
    c2, s2 = tmp1/tmp2, -p1[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]
    Cab_1 = rot1*rot2

    tmp1 = sqrt(p6[1]^2 + p6[2]^2)
    c1, s1 = p6[1]/tmp1, p6[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]
    tmp2 = sqrt(p6[1]^2 + p6[2]^2 + p6[3]^2)
    c2, s2 = tmp1/tmp2, p6[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]
    Cab_2 = rot1*rot2

    L_b1 = norm(p2-p1); r_b1 = p1; nelem_b1 = 5
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1, frame=Cab_1)
    compliance_b1 = fill(Diagonal([1.05204e-9, 3.19659e-9, 2.13106e-8, 1.15475e-7, 1.52885e-7, 7.1672e-9]), nelem_b1)

    L_b2 = norm(p3-p2); r_b2 = p2; nelem_b2 = 5
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=Cab_1)
    compliance_b2 = fill(Diagonal([1.24467e-9, 3.77682e-9, 2.51788e-8, 1.90461e-7, 2.55034e-7, 1.18646e-8]), nelem_b2)

    L_b3 = norm(p4-p3); r_b3 = p3; nelem_b3 = 5
    lengths_b3, xp_b3, xm_b3, Cab_b3 = discretize_beam(L_b3, r_b3, nelem_b3, frame=Cab_1)
    compliance_b3 = fill(Diagonal([1.60806e-9, 4.86724e-9, 3.24482e-8, 4.07637e-7, 5.57611e-7, 2.55684e-8]), nelem_b3)

    L_b4 = norm(p5-p4); r_b4 = p4; nelem_b4 = 5
    lengths_b4, xp_b4, xm_b4, Cab_b4 = discretize_beam(L_b4, r_b4, nelem_b4, frame=Cab_1)
    compliance_b4 = fill(Diagonal([2.56482e-9, 7.60456e-9, 5.67609e-8, 1.92171e-6, 2.8757e-6, 1.02718e-7]), nelem_b4)

    L_b5 = norm(p6-p5); r_b5 = p5; nelem_b5 = 20
    lengths_b5, xp_b5, xm_b5, Cab_b5 = discretize_beam(L_b5, r_b5, nelem_b5, frame=Cab_2)
    compliance_b5 = fill(Diagonal([2.77393e-9, 7.60456e-9, 1.52091e-7, 1.27757e-5, 2.7835e-5, 1.26026e-7]), nelem_b5)

    nelem = nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + nelem_b5
    points = vcat(xp_b1, xp_b2[2:end], xp_b3[2:end], xp_b4[2:end], xp_b5[2:end])
    start = 1:nelem; stop = 2:nelem+1
    lengths = vcat(lengths_b1, lengths_b2, lengths_b3, lengths_b4, lengths_b5)
    midpoints = vcat(xm_b1, xm_b2, xm_b3, xm_b4, xm_b5)
    Cab = vcat(Cab_b1, Cab_b2, Cab_b3, Cab_b4, Cab_b5)
    compliance = vcat(compliance_b1, compliance_b2, compliance_b3, compliance_b4, compliance_b5)

    assembly = Assembly(points, start, stop, compliance=compliance,
        frames=Cab, lengths=lengths, midpoints=midpoints)

    # Trimmed from length=141 to length=21 for benchmark feasibility
    Fz = range(0, 70e3, length=21)
    system = StaticSystem(assembly)

    point4 = nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1

    pcs_force = [
        Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            point4 => PrescribedConditions(Fz=Fz[i]),
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        ) for i = 1:length(Fz)]

    pcs_follower = [
        Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            point4 => PrescribedConditions(Fz_follower=Fz[i]),
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        ) for i = 1:length(Fz)]

    return (; assembly, system, pcs_force, pcs_follower)
end

function solve_static_joined_wing_linear(setup)
    reset_state!(setup.system)
    for i = 1:length(setup.pcs_force)
        static_analysis!(setup.system, setup.assembly, prescribed_conditions=setup.pcs_force[i], linear=true)
    end
end

function solve_static_joined_wing_nonlinear(setup)
    reset_state!(setup.system)
    for i = 1:length(setup.pcs_force)
        static_analysis!(setup.system, setup.assembly, prescribed_conditions=setup.pcs_force[i], reset_state=false)
    end
end

function solve_static_joined_wing_follower(setup)
    reset_state!(setup.system)
    for i = 1:length(setup.pcs_follower)
        static_analysis!(setup.system, setup.assembly, prescribed_conditions=setup.pcs_follower[i], reset_state=false)
    end
end

# ============================================================================
# 7. rotating — split into steady-state sweep + eigenvalue sweep
#
# Steady-state: trimmed from 31 rpm values to 7 (nonlinear only).
# Eigenvalue:   trimmed from 19×3 to 4×2 = 8 analyses.
# ============================================================================

function setup_rotating_assembly_short()
    # Short beam (16 elements) for steady-state sweep
    L_b1 = 31.5; r_b1 = [2.5, 0, 0]; nelem_b1 = 13
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)
    L_b2 = 6.0; r_b2 = [34.0, 0, 0]; nelem_b2 = 3
    sweep0 = 45 * pi/180
    cs, ss = cos(sweep0), sin(sweep0)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem; stop = 2:nelem+1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    w = 1.0; h = 0.063
    E = 1.06e7; ν = 0.325; ρ = 2.51e-4
    ky = 1.2000001839588001; kz = 14.625127919304001; kt = 65.85255016982444
    A = h*w; Iyy = w*h^3/12; Izz = w^3*h/12; J = Iyy + Izz
    Ay = A/ky; Az = A/kz; Jx = J/kt
    G = E/(2*(1+ν))
    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)
    mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)
    return assembly, E, ν, ρ, A, Iyy, Izz, J, Ay, Az, Jx, G
end

function setup_rotating_assembly_long(sweep)
    # Long beam (40 elements) for eigenvalue sweep
    L_b1 = 31.5; r_b1 = [2.5, 0, 0]; nelem_b1 = 20
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)
    L_b2 = 6.0; r_b2 = [34.0, 0, 0]; nelem_b2 = 20
    cs, ss = cos(sweep), sin(sweep)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem; stop = 2:nelem+1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    w = 1.0; h = 0.063
    E = 1.06e7; ν = 0.325; ρ = 2.51e-4
    ky = 1.2000001839588001; kz = 14.625127919304001; kt = 65.85255016982444
    A = h*w; Iyy = w*h^3/12; Izz = w^3*h/12; J = Iyy + Izz
    Ay = A/ky; Az = A/kz; Jx = J/kt
    G = E/(2*(1+ν))
    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)
    mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    return Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)
end

function setup_rotating_steady_state()
    assembly, _, _, _, _, _, _, _, _, _, _, _ = setup_rotating_assembly_short()
    prescribed_conditions = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    # Trimmed: 7 rpm values
    rpm_vals = 0:125:750
    return (; assembly, prescribed_conditions, rpm_vals)
end

function solve_rotating_steady_state(setup)
    for rpm in setup.rpm_vals
        w0 = [0.0, 0.0, rpm*(2*pi)/60]
        steady_state_analysis(setup.assembly,
            angular_velocity=w0,
            prescribed_conditions=setup.prescribed_conditions)
    end
end

function setup_rotating_eigenvalue()
    prescribed_conditions = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    # Trimmed: 4 sweep angles × 2 rpm values
    sweep_vals = (0:15:45) * pi/180
    rpm_vals = [0, 750]
    nev = 30
    # Pre-build one assembly per sweep
    assemblies = [setup_rotating_assembly_long(s) for s in sweep_vals]
    return (; assemblies, prescribed_conditions, rpm_vals, nev)
end

function solve_rotating_eigenvalue(setup)
    for assembly in setup.assemblies
        system = DynamicSystem(assembly)
        for rpm in setup.rpm_vals
            w0 = [0.0, 0.0, rpm*(2*pi)/60]
            eigenvalue_analysis!(system, assembly,
                angular_velocity=w0,
                prescribed_conditions=setup.prescribed_conditions,
                nev=setup.nev)
        end
    end
end

# ============================================================================
# 8. wind-turbine-blade — time_domain_analysis (split: undamped / damped)
# ============================================================================

function setup_wind_turbine_blade()
    L = 60; nelem = 5
    x = range(0, L, length=nelem+1); y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]
    start = 1:nelem; stop = 2:nelem+1
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8], nelem)
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418], nelem)
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)
    tvec = 0:0.001:0.1
    prescribed_conditions = (t) -> Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t)))
    return (; assembly, tvec, prescribed_conditions)
end

function solve_wind_turbine_blade_undamped(setup)
    time_domain_analysis(setup.assembly, setup.tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=false)
end

function solve_wind_turbine_blade_damped(setup)
    time_domain_analysis(setup.assembly, setup.tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=true)
end

# ============================================================================
# 9. dynamic-joined-wing — time_domain_analysis (split: undamped / damped)
# ============================================================================

function setup_dynamic_joined_wing()
    p1 = [0.0, 0, 0]; p2 = [-7.1726, -12, -3.21539]; p3 = [7.1726, -12,  3.21539]
    Cab_1 = [0.5 0.866025 0.0; 0.836516 -0.482963 0.258819; 0.224144 -0.12941 -0.965926]
    Cab_2 = [0.5 0.866025 0.0; -0.836516 0.482963 0.258819; 0.224144 -0.12941 0.965926]

    L_b1 = norm(p1-p2); r_b1 = p2; nelem_b1 = 4
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1, frame=Cab_1)
    L_b2 = norm(p3-p1); r_b2 = p1; nelem_b2 = 4
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=Cab_2)

    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem; stop = 2:nelem+1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)
    compliance = fill(Diagonal([2.93944738387698e-10, 8.42991725049126e-10, 3.38313996669689e-08,
        4.69246721094557e-08, 6.79584100559513e-08, 1.37068861370898e-09]), nelem)
    mass = fill(Diagonal([4.86e-2, 4.86e-2, 4.86e-2,
        1.0632465e-2, 2.10195e-4, 1.042227e-2]), nelem)
    assembly = Assembly(points, start, stop; compliance=compliance, mass=mass,
        frames=Cab, lengths=lengths, midpoints=midpoints)

    tvec = range(0, 0.01, length=201)

    F_L = (t) -> begin
        if 0.0 <= t < 0.01
            1e6*t
        elseif 0.01 <= t < 0.02
            -1e6*(t-0.02)
        else
            zero(t)
        end
    end
    F_S = (t) -> begin
        if 0.0 <= t < 0.02
            5e3*(1-cos(pi*t/0.02))
        else
            1e4
        end
    end
    prescribed_conditions = (t) -> Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        nelem_b1 + 1 => PrescribedConditions(Fx=F_L(t), Fy=F_L(t), Fz=F_S(t)),
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

    return (; assembly, tvec, prescribed_conditions)
end

function solve_dynamic_joined_wing_undamped(setup)
    time_domain_analysis(setup.assembly, setup.tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=false)
end

function solve_dynamic_joined_wing_damped(setup)
    time_domain_analysis(setup.assembly, setup.tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=true)
end

# ============================================================================
# 10. impliciteuler — GXBeam.simulate
# ============================================================================

function setup_impliciteuler()
    # Same assembly as wind-turbine-blade
    return setup_wind_turbine_blade()
end

function solve_impliciteuler(setup)
    GXBeam.simulate(setup.assembly, setup.tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=true)
end

# ============================================================================
# 11. stepsystem — manual initialize_system! + step_system! loop
# ============================================================================

function setup_stepsystem()
    return setup_wind_turbine_blade()
end

function solve_stepsystem(setup)
    system = GXBeam.DynamicSystem(setup.assembly)
    tvec = setup.tvec
    nt = length(tvec)
    system, gxstate, constants, paug, xgx, _ = GXBeam.initialize_system!(
        system, setup.assembly, tvec;
        prescribed_conditions=setup.prescribed_conditions,
        structural_damping=true, reset_state=true)
    prev = gxstate[1]
    for i in 2:nt
        system, prev, constants, paug, xgx, _ = GXBeam.step_system!(
            system, paug, xgx, constants, prev, setup.assembly, tvec, i;
            prescribed_conditions=setup.prescribed_conditions,
            structural_damping=true)
    end
end

# ============================================================================
# Driver
# ============================================================================

const WORKLOADS = [
    ("cantilever",                   setup_cantilever,                solve_cantilever,                30, 30),
    ("curved",                       setup_curved,                    solve_curved,                    20, 30),
    ("overdetermined",               setup_overdetermined,            solve_overdetermined,            30, 30),
    ("tipforce",                     setup_tipforce,                  solve_tipforce,                  10, 60),
    ("tipmoment",                    setup_tipmoment,                 solve_tipmoment,                 10, 60),
    ("static-joined-wing-linear",    setup_static_joined_wing,        solve_static_joined_wing_linear,     5, 90),
    ("static-joined-wing-nonlinear", setup_static_joined_wing,        solve_static_joined_wing_nonlinear,  5, 90),
    ("static-joined-wing-follower",  setup_static_joined_wing,        solve_static_joined_wing_follower,   5, 90),
    ("rotating-steady-state",        setup_rotating_steady_state,     solve_rotating_steady_state,         5, 90),
    ("rotating-eigenvalue",          setup_rotating_eigenvalue,       solve_rotating_eigenvalue,           3, 120),
    ("wind-turbine-blade-undamped",  setup_wind_turbine_blade,        solve_wind_turbine_blade_undamped,   3, 120),
    ("wind-turbine-blade-damped",    setup_wind_turbine_blade,        solve_wind_turbine_blade_damped,     3, 120),
    ("dynamic-joined-wing-undamped", setup_dynamic_joined_wing,       solve_dynamic_joined_wing_undamped,  3, 120),
    ("dynamic-joined-wing-damped",   setup_dynamic_joined_wing,       solve_dynamic_joined_wing_damped,    3, 120),
    ("impliciteuler",                setup_impliciteuler,             solve_impliciteuler,                 3, 120),
    ("stepsystem",                   setup_stepsystem,                solve_stepsystem,                    3, 120),
]

println("EXAMPLES BENCHMARK SUITE")
println("="^70)
println()

# Warmup pass: one untimed call per workload to JIT-compile
println("Warming up (JIT compilation)...")
for (name, setup_fn, solve_fn, _, _) in WORKLOADS
    print("  $name... ")
    s = setup_fn()
    solve_fn(s)
    println("ok")
end
println()

# Timed pass
results = Tuple{String, BenchmarkTools.Trial}[]
for (name, setup_fn, solve_fn, n_samples, n_seconds) in WORKLOADS
    println("--- $name (samples=$n_samples, seconds=$n_seconds) ---")
    s = setup_fn()
    bm = @benchmarkable $solve_fn($s)
    bm.params.samples = n_samples
    bm.params.seconds = n_seconds
    bm.params.evals = 1
    b = run(bm)
    display(b); println()
    push!(results, (name, b))
end

# Summary
println()
println("="^110)
println("SUMMARY")
println("="^110)
println()
pretty_table(results)
println()
