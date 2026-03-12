"""
    jacobian_benchmark.jl

Benchmark script comparing SparseDiffTools and DifferentiationInterface
for sparse Jacobian computation on a cantilever beam problem.

Run from the GXBeam root directory:
    julia --project=test test/jacobian_benchmark.jl
"""

using GXBeam, LinearAlgebra, Printf, Statistics, BenchmarkTools, ForwardDiff, SparseDiffTools
using DifferentiationInterface, SparseMatrixColorings, StaticArrays

# ============================================================================
# Setup: Create cantilever beam assembly
# ============================================================================

function create_cantilever_assembly(nelem=12)
    a = 0.3; b = 0.7; L = 1.0
    n1 = n3 = div(nelem, 3)
    n2 = nelem - n1 - n3
    x1 = range(0, a, length=n1+1)
    x2 = range(a, b, length=n2+1)
    x3 = range(b, L, length=n3+1)
    x = vcat(x1, x2[2:end], x3[2:end])
    y = zero(x); z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]
    start = 1:nelem; stop = 2:nelem+1
    EI = 1e9
    stiffness = fill(Diagonal([0, 0, 0, 0, EI, 0]), nelem)
    assembly = Assembly(points, start, stop, stiffness=stiffness)
    prescribed_conditions = Dict(
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )
    q = 1000
    distributed_loads = Dict()
    for ielem in n1+1:n1+n2
        distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz=(s)->q)
    end
    return assembly, prescribed_conditions, distributed_loads
end

assembly, prescribed_conditions, distributed_loads = create_cantilever_assembly(12)

system, state, converged = static_analysis(
    assembly,
    prescribed_conditions = prescribed_conditions,
    distributed_loads = distributed_loads,
    linear = true
)

# Extract state vector and parameters
x0 = copy(system.x)
nstate = length(x0)
p = nothing

# Build constants namedtuple matching what static_analysis! constructs internally.
# static_residual! only unpacks: indices, two_dimensional, force_scaling, assembly,
# prescribed_conditions, distributed_loads, point_masses, gravity, xpfunc, pfunc, t
constants = (;
    assembly,
    indices          = system.indices,
    two_dimensional  = false,
    force_scaling    = system.force_scaling,
    xpfunc           = nothing,
    pfunc            = (p, t) -> (;),
    t                = 0.0,
    prescribed_conditions,
    distributed_loads,
    point_masses     = Dict{Int, PointMass{Float64}}(),
    gravity          = (@SVector zeros(3)),
)

# Preallocated dense Jacobian storage (for SparseDiffTools)
jacob = zeros(nstate, nstate)

println("Benchmark Configuration:")
println("  Problem size (state variables): $nstate")
println("  Assembly: cantilever with 12 elements")
println()

# ============================================================================
# One-time setup: compute coloring/prep outside of "compute-only" benchmarks
# ============================================================================

# Out-of-place wrapper (for ForwardDiff.jacobian and sparsity detection)
function residual_oop(x)
    r = zeros(eltype(x), nstate)
    GXBeam.static_residual!(r, x, p, constants)
    return r
end

# In-place wrapper used for ALL DI calls — must be the same object for prep + exec
struct ResidualInPlace{P,C}
    p::P
    c::C
end
(f::ResidualInPlace)(r, x) = GXBeam.static_residual!(r, x, f.p, f.c)

f_di! = ResidualInPlace(p, constants)

# SparseDiffTools: compute sparsity pattern + colors once
J_dense = ForwardDiff.jacobian(residual_oop, x0)
nnz_count = count(!iszero, J_dense)
# Threshold to get structural sparsity pattern
J_pattern = abs.(J_dense) .> (maximum(abs, J_dense) * 1e-12)
colors_st = SparseDiffTools.matrix_colors(J_pattern)

# DifferentiationInterface: prepare sparse backend once
di_backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector  = DenseSparsityDetector(AutoForwardDiff(); atol=1e-5),
    coloring_algorithm = GreedyColoringAlgorithm(),
)
resid_di = zeros(nstate)
prep_di = prepare_jacobian(f_di!, resid_di, di_backend, x0)

# DI matrix-free: prepare pushforward and pullback once
v0 = randn(nstate)   # tangent vector
w0 = randn(nstate)   # cotangent vector
dy0 = zeros(nstate)
dx0 = zeros(nstate)

mf_backend = AutoForwardDiff()

prep_pf = prepare_pushforward(f_di!, dy0, mf_backend, x0, (v0,))
prep_pb = prepare_pullback(f_di!, dx0, mf_backend, x0, (w0,))

println("Sparsity info:")
println("  Jacobian nonzeros: $nnz_count / $(nstate*nstate) ($(round(100*nnz_count/(nstate*nstate),digits=1))%)")
println("  Colors (SparseDiffTools): $(maximum(colors_st)) for $nstate DOFs")
println()

# ============================================================================
# Benchmark functions
# ============================================================================

# --- SparseDiffTools: compute only (pre-computed colors) ---
function bench_st_compute!(jacob, x0, p, constants, colors)
    fill!(jacob, 0)
    GXBeam.autodiff_jacobian!(jacob, GXBeam.static_residual!, x0, p, constants; colors=colors)
    return jacob
end

# --- SparseDiffTools: full pipeline (pattern detection + coloring + compute) ---
function bench_st_full(x0, p, constants)
    J_local = ForwardDiff.jacobian(
        x -> (r = zeros(eltype(x), nstate); GXBeam.static_residual!(r, x, p, constants); r),
        x0
    )
    colors = SparseDiffTools.matrix_colors(J_local)
    J_out = zeros(nstate, nstate)
    GXBeam.autodiff_jacobian!(J_out, GXBeam.static_residual!, x0, p, constants; colors=colors)
    return J_out
end

# --- DifferentiationInterface: compute only (pre-computed prep) ---
# f_di! must be the SAME object used in prepare_jacobian (DI checks type identity)
function bench_di_compute!(f!, resid, prep, backend, x0)
    return jacobian(f!, resid, prep, backend, x0)
end

# --- DifferentiationInterface: full pipeline (sparsity detection + coloring + compute) ---
function bench_di_full(f!, x0)
    backend = AutoSparse(
        AutoForwardDiff();
        sparsity_detector  = DenseSparsityDetector(AutoForwardDiff(); atol=1e-5),
        coloring_algorithm = GreedyColoringAlgorithm(),
    )
    r_local = zeros(nstate)
    prep = prepare_jacobian(f!, r_local, backend, x0) 
    return jacobian(f!, r_local, prep, backend, x0)
end

# --- SparseDiffTools: construct JacVec operator ---
function bench_st_jacvec_construct(x0, p, constants)
    return GXBeam.matrixfree_jacobian(GXBeam.static_residual!, x0, p, constants)
end

# --- SparseDiffTools JacVec: apply J*v via mul! ---
function bench_st_jacvec_apply!(w, JV, v)
    mul!(w, JV, v)
    return w
end

# --- DI pushforward: compute JVP (pre-computed prep) ---
function bench_di_pushforward!(dy, f!, prep, backend, x, v)
    return pushforward(f!, dy, prep, backend, x, (v,))
end

# --- DI pullback: compute VJP (pre-computed prep) ---
function bench_di_pullback!(dx, f!, prep, backend, x, w)
    return pullback(f!, dx, prep, backend, x, (w,))
end

# ============================================================================
# Warmup
# ============================================================================

println("Warming up (JIT compilation)...")
bench_st_compute!(jacob, x0, p, constants, colors_st)
bench_st_full(x0, p, constants)
bench_di_compute!(f_di!, resid_di, prep_di, di_backend, x0)
bench_di_full(f_di!, x0)
JV = bench_st_jacvec_construct(x0, p, constants)
w_buf = zeros(nstate)
bench_st_jacvec_construct(x0, p, constants)
bench_st_jacvec_apply!(w_buf, JV, v0)
bench_di_pushforward!(dy0, f_di!, prep_pf, mf_backend, x0, v0)
bench_di_pullback!(dx0, f_di!, prep_pb, mf_backend, x0, w0)
println("Done.")
println()

# Sanity check: JVP from ST and DI should agree
err = norm(w_buf - dy0)
println("JVP sanity check (ST mul! vs DI pushforward): norm(w_buf - dy0) = $(round(err, sigdigits=3))")
println()

# ============================================================================
# Run benchmarks
# ============================================================================

println("=" ^ 70)
println("JACOBIAN BENCHMARK: SparseDiffTools vs DifferentiationInterface")
println("=" ^ 70)
println()

println("--- Compute-only (sparsity/coloring pre-computed) ---")
println()

println("SparseDiffTools (colored ForwardDiff, compute only)...")
bst_c = @benchmark bench_st_compute!($jacob, $x0, $p, $constants, $colors_st) samples=30 seconds=30
display(bst_c); println()

println("DifferentiationInterface (sparse ForwardDiff, compute only)...")
bdi_c = @benchmark bench_di_compute!($f_di!, $resid_di, $prep_di, $di_backend, $x0) samples=30 seconds=30
display(bdi_c); println()

println("--- Full pipeline (includes sparsity detection + coloring) ---")
println()

println("SparseDiffTools (full pipeline)...")
bst_f = @benchmark bench_st_full($x0, $p, $constants) samples=10 seconds=60
display(bst_f); println()

println("DifferentiationInterface (full pipeline)...")
bdi_f = @benchmark bench_di_full($f_di!, $x0) samples=10 seconds=60
display(bdi_f); println()

println("--- Matrix-free Jacobian-vector products ---")
println()

println("SparseDiffTools JacVec construction...")
bst_mf_construct = @benchmark bench_st_jacvec_construct($x0, $p, $constants) samples=30 seconds=30
display(bst_mf_construct); println()

println("SparseDiffTools JacVec apply (J*v)...")
bst_mf_apply = @benchmark bench_st_jacvec_apply!($w_buf, $JV, $v0) samples=100 seconds=30
display(bst_mf_apply); println()

println("DifferentiationInterface pushforward J*v (pre-computed prep)...")
bdi_pf = @benchmark bench_di_pushforward!($dy0, $f_di!, $prep_pf, $mf_backend, $x0, $v0) samples=100 seconds=30
display(bdi_pf); println()

println("DifferentiationInterface pullback J'*v (pre-computed prep)...")
bdi_pb = @benchmark bench_di_pullback!($dx0, $f_di!, $prep_pb, $mf_backend, $x0, $w0) samples=100 seconds=30
display(bdi_pb); println()

# ============================================================================
# Summary table
# ============================================================================

println("=" ^ 110)
println("SUMMARY")
println("=" ^ 110)
println()

# Helper: extract stats in ms from a BenchmarkTools trial
ms(t) = t / 1e6
function stats_ms(b)
    ts = b.times
    return (
        min  = ms(minimum(ts)),
        mean = ms(mean(ts)),
        std  = ms(std(ts)),
        med  = ms(median(ts)),
        max  = ms(maximum(ts)),
    )
end

rows = [
    ("SparseDiffTools compute-only",      bst_c),
    ("DiffInterface   compute-only",      bdi_c),
    ("SparseDiffTools full pipeline",     bst_f),
    ("DiffInterface   full pipeline",     bdi_f),
    ("SparseDiffTools JacVec construct",  bst_mf_construct),
    ("SparseDiffTools JacVec apply  J*v", bst_mf_apply),
    ("DiffInterface   pushforward   J*v", bdi_pf),
    ("DiffInterface   pullback      J'*v",bdi_pb),
]

println("| Method | Min (ms) | Mean (ms) | Std (ms) | Med (ms) | Max (ms) | Allocs | Mem (MB) |")
println("|---|---|---|---|---|---|---|---|")
for (label, b) in rows
    s = stats_ms(b)
    @printf "| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %d | %.2f |\n" label s.min s.mean s.std s.med s.max b.allocs b.memory/1e6
end
println()

st_c_min = minimum(bst_c.times);       di_c_min  = minimum(bdi_c.times)
st_f_min = minimum(bst_f.times);       di_f_min  = minimum(bdi_f.times)
st_apply_min = minimum(bst_mf_apply.times); pf_min = minimum(bdi_pf.times)
println("**Compute-only speedup (ST/DI):**   $(round(ms(st_c_min)/ms(di_c_min), digits=2))x — $(st_c_min < di_c_min ? "ST faster" : "DI faster")")
println("**Full pipeline speedup (ST/DI):**   $(round(ms(st_f_min)/ms(di_f_min), digits=2))x — $(st_f_min < di_f_min ? "ST faster" : "DI faster")")
println("**JVP speedup (ST apply / DI pf):**  $(round(ms(st_apply_min)/ms(pf_min), digits=2))x — $(st_apply_min < pf_min ? "ST faster" : "DI faster")")
println()
