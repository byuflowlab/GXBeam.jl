@testset "DifferentialEquations" begin

    L = 60 # m

    # create points
    nelem = 10
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix for each beam element
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix for each beam element
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    # prescribed conditions
    prescribed_conditions = (t) -> begin
        Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t))
            )
    end

    # define simulation time
    tspan = (0.0, 0.1)

    # define named tuple with parameters
    p = (; prescribed_conditions)

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions)

    # construct ODEProblem
    prob = ODEProblem(system, assembly, tspan, p; 
        constant_mass_matrix = false,
        structural_damping = false)

    # solve ODEProblem
    sol = solve(prob, Rodas4())

    # test solution convergence
    @test sol.t[end] == 0.1

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions,
        constant_mass_matrix = false,
        structural_damping = true)

    # construct ODEProblem
    prob = ODEProblem(system, assembly, tspan, p; 
        constant_mass_matrix = false,
        structural_damping = true)

    # solve ODEProblem
    sol = solve(prob, Rodas4())

    # test solution convergence
    @test sol.t[end] == 0.1

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions,
        constant_mass_matrix = true,
        structural_damping = false)

    # construct ODEProblem
    prob = ODEProblem(system, assembly, tspan, p; 
        constant_mass_matrix = true,
        structural_damping = false)

    # solve ODEProblem
    sol = solve(prob, Rodas4())

    # test solution convergence
    @test sol.t[end] == 0.1

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions,
        constant_mass_matrix = true,
        structural_damping = false)

    # construct ODEProblem
    prob = ODEProblem(system, assembly, tspan, p; 
        constant_mass_matrix = true,
        structural_damping = true)

    # solve ODEProblem
    sol = solve(prob, Rodas4())

    # test solution convergence
    @test sol.t[end] == 0.1

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions,
        structural_damping = false)

    # construct DAEProblem
    prob = DAEProblem(system, assembly, tspan, p; 
        structural_damping = false)

    # solve DAEProblem
    sol = solve(prob, DABDF2())

    # test solution convergence
    @test sol.t[end] == 0.1

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; 
        prescribed_conditions = prescribed_conditions,
        structural_damping = false)

    # construct DAEProblem
    prob = DAEProblem(system, assembly, tspan, p; 
        structural_damping = true)

    # solve DAEProblem
    sol = solve(prob, DABDF2())

    # test solution convergence
    @test sol.t[end] == 0.1
end