#=
Type-stability check for the DI-migrated AD path.

Standalone — not included in runtests.jl. Run with:

    julia --project=test test/type_stability.jl
    julia --project=test test/type_stability.jl --verbose

Three layers, increasing tolerance for `Any`:

  1. Function-barrier interiors (`_autodiff_jacobian_inner!`, `_build_linearmap`):
     hard `@inferred` assertions — these MUST be concretely inferable, that's
     the entire point of the barriers. We construct the cached NamedTuples
     directly here (using DI on a toy residual) so the test doesn't depend on
     which analysis path GXBeam happens to dispatch into.

  2. Outer DI wrappers via `@code_warntype` summary. Expected to have some
     `::Any` because `system.prep_jacobian` is `::Any` by design; counts are
     printed so a regression bumps the number visibly.

  3. Top-level solve entries via `@code_warntype` summary. These have many
     kwargs and Dict-typed params, so meaningful `::Any` counts are normal;
     numbers exist for trend tracking, not assertion.

Exit code 0 if layer-1 assertions pass; nonzero otherwise. Layers 2 and 3
are informational.
=#

using GXBeam
using LinearAlgebra
using Test
using InteractiveUtils
using DifferentiationInterface
const DI = DifferentiationInterface

# const VERBOSE = "--verbose" in ARGS || "-v" in ARGS
const VERBOSE = false 

# ---------------------------------------------------------------------------
# Build assembly so we have a real Assembly type for Layer 3 warntype scans
# ---------------------------------------------------------------------------

function build_assembly()
    nelem = 6
    L = 1.0
    x = collect(range(0, L; length=nelem+1))
    y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i in eachindex(x)]
    EI = 1e9
    stiffness = fill(Diagonal([0, 0, 0, 0, EI, 0]), nelem)
    assembly = Assembly(points, 1:nelem, 2:nelem+1; stiffness=stiffness)
    pc = Dict(nelem + 1 => PrescribedConditions(
        ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    return assembly, pc
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function warntype_string(f, argtypes)
    io = IOBuffer()
    ctx = IOContext(io, :color => false)
    InteractiveUtils.code_warntype(ctx, f, argtypes)
    return String(take!(io))
end

function summarize_warntype(label, text)
    # NOTE: when `code_warntype` writes to a no-color IOContext it marks
    # unstable types as ::ANY (uppercase) instead of using ANSI red on
    # ::Any. We match both. Also strip any leftover ANSI escape codes so
    # a colored capture matches the same way.
    cleaned    = replace(text, r"\e\[[0-9;]*m" => "")
    any_count  = count(s -> occursin("::ANY", s) || occursin("::Any", s),
                       split(cleaned, '\n'))
    body_any   = occursin(r"Body::(ANY|Any)\b", cleaned) ? 1 : 0
    body_union = count(s -> occursin(r"^Body::Union\{", s), split(cleaned, '\n'))
    println(rpad(label, 50), "  any=", any_count,
            "  body_any=", body_any, "  body_union=", body_union)
    if VERBOSE
        println("--- $label ---")
        println(cleaned)
        println()
    end
    return (; any_count, body_any, body_union)
end

# ---------------------------------------------------------------------------
# Build a `cached` NamedTuple of the exact shape stored in
# `system.prep_jacobian`. We mirror `autodiff_jacobian!` (src/analyses.jl
# around line 3850-3862) so the type matches what production code sees.
# ---------------------------------------------------------------------------

function build_autodiff_cache(::Type{T}, n::Int) where {T}
    # Toy residual mirroring the (r, x, p, c) shape used in GXBeam:
    # accepts a residual buffer, state x, params p, captures c. The
    # function-barrier interior only depends on the cached NamedTuple's
    # type, so the inside of the residual doesn't matter.
    residual! = (r, x, p, c) -> (r .= x; r)

    functor = GXBeam.ResidualWithCapture(residual!, nothing, nothing)
    backend = DI.AutoSparse(
        DI.AutoForwardDiff();
        sparsity_detector  = DI.DenseSparsityDetector(DI.AutoForwardDiff(); atol=1e-5),
        coloring_algorithm = GXBeam.SparseMatrixColorings.GreedyColoringAlgorithm(),
    )
    x_ref = zeros(T, n)
    r_buf = similar(x_ref)
    prep  = DI.prepare_jacobian(functor, r_buf, backend, x_ref)
    key = (typeof(residual!), Nothing, Nothing, n)
    return (; key=key, functor=functor, prep=prep,
              backend=backend, r_buf=r_buf)
end

function build_matrixfree_cache(::Type{T}, n::Int) where {T}
    residual! = (r, x, p, c) -> (r .= x; r)
    functor = GXBeam.ResidualWithCapture(residual!, nothing, nothing)
    backend = DI.AutoForwardDiff()

    x_ref = zeros(T, n)
    v0  = zeros(T, n); w0  = zeros(T, n)
    dy0 = zeros(T, n); dx0 = zeros(T, n)
    prep_pf = DI.prepare_pushforward(functor, dy0, backend, x_ref, (v0,))
    prep_pb = DI.prepare_pullback(functor,  dx0, backend, x_ref, (w0,))

    key = (typeof(residual!), Nothing, Nothing, n)
    return (; key=key, functor=functor, backend=backend,
              prep_pf=prep_pf, prep_pb=prep_pb, n=n,
              v_buf=zeros(T, n), w_buf=zeros(T, n),
              dy_buf=zeros(T, n), dx_buf=zeros(T, n))
end

# ---------------------------------------------------------------------------
# Layer 1: hard @inferred assertions on the function-barrier interiors
# ---------------------------------------------------------------------------

println("\n=== Layer 1: function-barrier interiors (hard inference assertions) ===")

@testset "_autodiff_jacobian_inner! is inferable" begin
    n = 8
    cached = build_autodiff_cache(Float64, n)
    x     = zeros(Float64, n)
    jacob = zeros(Float64, n, n)
    # @inferred throws if inference can't pin down a concrete return type.
    @test (@inferred GXBeam._autodiff_jacobian_inner!(jacob, x, cached)) === jacob
end

@testset "_build_linearmap is inferable" begin
    n = 8
    cached = build_matrixfree_cache(Float64, n)
    x = zeros(Float64, n)
    L = @inferred GXBeam._build_linearmap(x, cached)
    @test size(L) == (n, n)
end

# ---------------------------------------------------------------------------
# Layer 1b: warm-cache check
#
# The Layer 1 tests above use synthetic caches built right here in the script.
# This block verifies the same `@inferred` guarantee against caches produced
# by *production* code paths — `autodiff_jacobian!` and `matrixfree_jacobian`
# populate `system.prep_jacobian` / `system.prep_jvp` on first call, and we
# then assert the function-barrier interiors stay inferable when handed the
# real stored cache.
#
# This catches a regression where the cache's *runtime* shape stops matching
# what `_autodiff_jacobian_inner!` / `_build_linearmap` expect, even if a
# synthetic mock still passes.
# ---------------------------------------------------------------------------

println("\n=== Layer 1b: warm-cache @inferred checks ===")

@testset "warm autodiff_jacobian! cache stays inferable" begin
    nelem = 6; L_len = 1.0
    xp = collect(range(0, L_len; length=nelem+1))
    yp = zero(xp); zp = zero(xp)
    pts = [[xp[i], yp[i], zp[i]] for i in eachindex(xp)]
    stiff = fill(Diagonal([0, 0, 0, 0, 1e9, 0]), nelem)
    asm = Assembly(pts, 1:nelem, 2:nelem+1; stiffness=stiff)
    pc = Dict(nelem + 1 => PrescribedConditions(
        ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    sys, _, _ = static_analysis(asm; prescribed_conditions=pc, linear=true)

    constants = (system=sys,)
    residual! = (r, x, p, c) -> (r .= x; r)
    n = length(sys.x)
    x  = zeros(Float64, n)
    K  = zeros(Float64, n, n)

    # Cold call: populates sys.prep_jacobian.
    GXBeam.autodiff_jacobian!(K, residual!, x, nothing, constants)
    @test sys.prep_jacobian isa NamedTuple

    # Now feed the warm cache through the function barrier and assert
    # inference still gives a concrete return type.
    cached_warm = sys.prep_jacobian
    @test (@inferred GXBeam._autodiff_jacobian_inner!(K, x, cached_warm)) === K
end

@testset "warm matrixfree_jacobian cache stays inferable" begin
    nelem = 6; L_len = 1.0
    xp = collect(range(0, L_len; length=nelem+1))
    yp = zero(xp); zp = zero(xp)
    pts = [[xp[i], yp[i], zp[i]] for i in eachindex(xp)]
    stiff = fill(Diagonal([0, 0, 0, 0, 1e9, 0]), nelem)
    asm = Assembly(pts, 1:nelem, 2:nelem+1; stiffness=stiff)
    pc = Dict(nelem + 1 => PrescribedConditions(
        ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    sys, _, _ = static_analysis(asm; prescribed_conditions=pc, linear=true)

    constants = (system=sys,)
    residual! = (r, x, p, c) -> (r .= x; r)
    x = zeros(Float64, length(sys.x))

    # Cold call: populates sys.prep_jvp and returns a LinearMap.
    L = GXBeam.matrixfree_jacobian(residual!, x, nothing, constants)
    @test sys.prep_jvp isa NamedTuple
    @test size(L, 1) == length(x)

    # Warm function-barrier check.
    cached_warm = sys.prep_jvp
    L2 = @inferred GXBeam._build_linearmap(x, cached_warm)
    @test size(L2) == size(L)
end

# ---------------------------------------------------------------------------
# Layer 2: outer wrapper warntype summary (informational)
# ---------------------------------------------------------------------------

println("\n=== Layer 2: outer wrapper warntype (informational) ===")
println("    `Any` from `system.prep_jacobian::Any` field access is expected;")
println("    function-barrier helpers should remain clean.")
println()

let n = 8
    cached_aj = build_autodiff_cache(Float64, n)
    x  = zeros(Float64, n)
    J  = zeros(Float64, n, n)
    text = warntype_string(GXBeam._autodiff_jacobian_inner!,
                           Tuple{typeof(J), typeof(x), typeof(cached_aj)})
    summarize_warntype("_autodiff_jacobian_inner!", text)

    # @code_warntype GXBeam._autodiff_jacobian_inner!(J, x, cached_aj) #checking julia's flags for good typing. 

    cached_mf = build_matrixfree_cache(Float64, n)
    text = warntype_string(GXBeam._build_linearmap,
                           Tuple{typeof(x), typeof(cached_mf)})
    summarize_warntype("_build_linearmap", text)

    # @code_warntype GXBeam._build_linearmap(x, cached_mf) #checking julia's flags for good typing.
end


# The outer wrappers take `residual!`, `x`, `p`, `constants` where `constants`
# is a NamedTuple specific to each analysis. We use a stub `constants` with a
# StaticSystem (matching the real production type) so `system.prep_jacobian`
# field access is honest.
let n = 8
    assembly, pc = build_assembly()
    sys, _, _ = static_analysis(assembly; prescribed_conditions=pc, linear=true)
    constants = (system=sys,)
    residual! = (r, x, p, c) -> (r .= x; r)
    x = zeros(Float64, n)
    J = zeros(Float64, n, n)

    text = warntype_string(GXBeam.autodiff_jacobian!,
        Tuple{typeof(J), typeof(residual!), typeof(x), Nothing, typeof(constants)})
    summarize_warntype("autodiff_jacobian! (outer wrapper)", text)

    # @code_warntype GXBeam.autodiff_jacobian!(J, residual!, x, nothing, constants) #checking julia's flags for good typing.

    text = warntype_string(GXBeam.matrixfree_jacobian,
        Tuple{typeof(residual!), typeof(x), Nothing, typeof(constants)})
    summarize_warntype("matrixfree_jacobian (outer wrapper)", text)

    # @code_warntype GXBeam.matrixfree_jacobian(residual!, x, nothing, constants) #checking julia's flags for good typing. #Claude says that the 18 Anys here are expected. The idea is that the outer wrapper is allowed to be type unstable, but the inner wrappers are stable. #Claude says there is a type instability here (inside `_build_linearmap`), and we can make a helper function to return a concretely typed `LinearMap`
end

# ---------------------------------------------------------------------------
# Layer 3: top-level solve entries (informational)
# ---------------------------------------------------------------------------

println("\n=== Layer 3: top-level solve entries (informational) ===")
println("    Heavy kwargs/Dict params here — nonzero `Any` counts are normal.")
println()

let
    assembly, pc = build_assembly()
    sys, _, _ = static_analysis(assembly; prescribed_conditions=pc, linear=true)

    text = warntype_string(static_analysis, Tuple{typeof(assembly)})
    summarize_warntype("static_analysis(Assembly)", text)

    # @code_warntype static_analysis(assembly) #checking julia's flags for good typing.

    text = warntype_string(static_analysis!, Tuple{typeof(sys), typeof(assembly)})
    summarize_warntype("static_analysis!(StaticSystem, Assembly)", text)

    # @code_warntype static_analysis!(sys, assembly) #checking julia's flags for good typing. 

    text = warntype_string(steady_state_analysis, Tuple{typeof(assembly)})
    summarize_warntype("steady_state_analysis(Assembly)", text)

    text = warntype_string(eigenvalue_analysis, Tuple{typeof(assembly)})
    summarize_warntype("eigenvalue_analysis(Assembly)", text)
end

println("\ntype_stability.jl done. Layer 1 failures (if any) are the only ones",
        "\nthat should block a release.")
