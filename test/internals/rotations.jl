using GXBeam, Random, ForwardDiff, Test

@testset "Rotation Parameters" begin
    
    RNG = MersenneTwister(1234)

    θ = 1e3*rand(RNG, 3)
    Δθ = 1e3*rand(RNG, 3)

    # get_C_θ
    C_θ1, C_θ2, C_θ3 = GXBeam.get_C_θ(θ)
    @test isapprox(C_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_C([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(C_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_C([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(C_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_C([θ[1], θ[2], θ3]), θ[3]))

    # get_Q_θ
    Q_θ1, Q_θ2, Q_θ3 = GXBeam.get_Q_θ(θ)
    @test isapprox(Q_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_Q([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(Q_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_Q([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(Q_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_Q([θ[1], θ[2], θ3]), θ[3]))

    # get_Qinv_θ
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = GXBeam.get_Qinv_θ(θ)
    @test isapprox(Qinv_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_Qinv([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(Qinv_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_Qinv([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(Qinv_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_Qinv([θ[1], θ[2], θ3]), θ[3]))

    # get_ΔQ
    ΔQ = GXBeam.get_ΔQ(θ, Δθ)
    @test isapprox(ΔQ, GXBeam.mul3(Q_θ1, Q_θ2, Q_θ3, Δθ))

    # get_ΔQ_θ
    ΔQ_θ1, ΔQ_θ2, ΔQ_θ3 = GXBeam.get_ΔQ_θ(θ, Δθ)
    @test isapprox(ΔQ_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_ΔQ([θ1, θ[2], θ[3]], Δθ), θ[1]))
    @test isapprox(ΔQ_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_ΔQ([θ[1], θ2, θ[3]], Δθ), θ[2]))
    @test isapprox(ΔQ_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_ΔQ([θ[1], θ[2], θ3], Δθ), θ[3]))

end