using GXBeam, LinearAlgebra, ForwardDiff, Test

@testset "ForwardDiff" begin

    # Linear Analysis of a Beam Under a Linear Distributed Load

    function linear_analysis_test_with_AD(length) # this should affect just about everything

        nelem = 16

        # create points
        L = length[1]
        x = collect(range(0, L, length=nelem+1))
        y = zero(x)
        z = zero(x)

        points = [[x[i],y[i],z[i]] for i = 1:size(x,1)]

        # index of endpoints for each beam element
        start = 1:nelem
        stop = 2:nelem+1

        # create compliance matrix for each beam element
        EI = 1e7
        compliance = fill(Diagonal([0, 0, 0, 0, 1/EI, 0]), nelem)

        # create assembly
        assembly = Assembly(points, start, stop, compliance=compliance)

        # set prescribed conditions
        prescribed_conditions = Dict(
            # simply supported left endpoint
            1 => PrescribedConditions(uz=0),
            # clamped right endpoint
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

        # create distributed load
        qmax = 1000
        distributed_loads = Dict()
        for i = 1:nelem
            distributed_loads[i] = DistributedLoads(assembly, i; s1=x[i],
                s2=x[i+1], fz = (s) -> qmax*s)
        end

        # solve system
        system, state, converged = static_analysis(assembly;
            prescribed_conditions = prescribed_conditions,
            distributed_loads = distributed_loads,
            linear = true)

        return system.x
    end

    # run FrowardDiff - no specific test, just make sure it runs fine
    J = ForwardDiff.jacobian(linear_analysis_test_with_AD, [1.0]) #length=1
end