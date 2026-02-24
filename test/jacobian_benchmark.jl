"""
    jacobian_benchmark.jl

Benchmark script comparing SparseDiffTools and DifferentiationInterface 
for sparse Jacobian computation on a realistic cantilever beam problem.

Includes the cost of sparsity pattern detection and coloring in the total runtime.
"""

using GXBeam, LinearAlgebra, BenchmarkTools, ForwardDiff, SparseDiffTools
using DifferentiationInterface

# ============================================================================
# Setup: Create cantilever beam assembly from test example
# ============================================================================

function create_cantilever_assembly(nelem=12)
    """Create a cantilever beam assembly for benchmarking."""
    
    # create points
    a = 0.3
    b = 0.7
    L = 1.0
    n1 = n3 = div(nelem, 3)
    n2 = nelem - n1 - n3
    x1 = range(0, a, length=n1+1)
    x2 = range(a, b, length=n2+1)
    x3 = range(b, L, length=n3+1)
    x = vcat(x1, x2[2:end], x3[2:end])
    y = zero(x)
    z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]

    # index of endpoints for each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # create compliance matrix for each beam element
    EI = 1e9
    stiffness = fill(Diagonal([0, 0, 0, 0, EI, 0]), nelem)

    # create the assembly
    assembly = Assembly(points, start, stop, stiffness=stiffness)

    # set prescribed conditions (fixed right endpoint)
    prescribed_conditions = Dict(
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

    # create distributed load
    q = 1000
    distributed_loads = Dict()
    for ielem in n1+1:n1+n2
        distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz = (s) -> q)
    end

    return assembly, prescribed_conditions, distributed_loads
end

# ============================================================================
# Get initial state and problem parameters
# ============================================================================

assembly, prescribed_conditions, distributed_loads = create_cantilever_assembly(12)

# Run analysis to get initial state
system, state, converged = static_analysis(
    assembly, 
    prescribed_conditions=prescribed_conditions,
    distributed_loads=distributed_loads, 
    linear=true
)

# Extract state vector and problem parameters
x = GXBeam.get_state_vector(system, state)
nstate = length(x)

# Create parameter dictionary
p = GXBeam.assemble_system(assembly, prescribed_conditions, distributed_loads)

# Residual vector (preallocated)
resid = zeros(nstate)

# Jacobian matrix (preallocated)
jacob = zeros(nstate, nstate)

println("Benchmark Configuration:")
println("  Problem size (state variables): $nstate")
println("  Assembly: cantilever with 12 elements")
println()

# ============================================================================
# Method 1: SparseDiffTools with colored forward differentiation
# ============================================================================

function benchmark_sparsetools_colored()
    """
    Benchmark SparseDiffTools colored forward-mode autodiff.
    Includes the cost of computing sparsity coloring.
    """
    
    # Step 1: Compute sparsity pattern via finite differences
    # (this is similar to jacobian_colors logic, but simplified)
    jacob_sample = similar(jacob)
    ForwardDiff.jacobian!(jacob_sample, (r, x) -> GXBeam.steady_residual!(r, x, p, assembly), resid, x)
    
    # Step 2: Compute colors for sparse jacobian
    colors = SparseDiffTools.matrix_colors(abs.(jacob_sample))
    
    # Step 3: Colored forward-mode jacobian computation
    fill!(jacob, 0)
    GXBeam.autodiff_jacobian!(jacob, GXBeam.steady_residual!, x, p, assembly; colors=colors)
    
    return jacob
end

# ============================================================================
# Method 2: DifferentiationInterface sparse jacobian
# ============================================================================

function benchmark_diffintrface()
    """
    Benchmark DifferentiationInterface sparse jacobian computation.
    Includes detection of sparsity pattern.
    """
    
    # Define residual function with correct signature for DifferentiationInterface
    function residual_wrapper(x)
        resid_local = zeros(nstate)
        GXBeam.steady_residual!(resid_local, x, p, assembly)
        return resid_local
    end
    
    # Compute jacobian using DifferentiationInterface
    # sparse_jacobian automatically detects sparsity pattern
    jacob_di = sparse_jacobian(residual_wrapper, x)
    
    return jacob_di
end

# ============================================================================
# Method 3: SparseDiffTools matrix-free (for comparison)
# ============================================================================

function benchmark_sparsetools_matrixfree()
    """
    Benchmark SparseDiffTools matrix-free jacobian-vector product.
    Note: This doesn't produce a dense matrix but a lazy operator.
    Included for reference on the third jacobian method used in GXBeam.
    """
    
    # Create matrix-free operator
    jvp_op = GXBeam.matrixfree_jacobian(GXBeam.steady_residual!, x, p, assembly)
    
    return jvp_op
end

# ============================================================================
# Run benchmarks
# ============================================================================

println("=" * 70)
println("JACOBIAN BENCHMARK: SparseDiffTools vs DifferentiationInterface")
println("=" * 70)
println()

# Warmup (compile)
println("Warming up (compilation)...")
benchmark_sparsetools_colored()
try
    benchmark_diffintrface()
catch e
    println("Warning: DifferentiationInterface benchmark encountered error: $e")
end
benchmark_sparsetools_matrixfree()
println("Warmup complete.")
println()

# SparseDiffTools colored forward-mode
println("Benchmarking SparseDiffTools (colored forward-mode)...")
bench_st = @benchmark benchmark_sparsetools_colored() samples=10 seconds=30
println("  Time: $(minimum(bench_st.times) / 1e6) ms")
println("  Allocations: $(bench_st.allocs) allocs, $(bench_st.memory / 1e6) MB")
println()

# DifferentiationInterface
println("Benchmarking DifferentiationInterface (sparse jacobian)...")
try
    bench_di = @benchmark benchmark_diffintrface() samples=10 seconds=30
    println("  Time: $(minimum(bench_di.times) / 1e6) ms")
    println("  Allocations: $(bench_di.allocs) allocs, $(bench_di.memory / 1e6) MB")
    println()
catch e
    println("  Error: $e")
    println()
end

# SparseDiffTools matrix-free (reference)
println("Benchmarking SparseDiffTools (matrix-free JacVec)...")
bench_mf = @benchmark benchmark_sparsetools_matrixfree() samples=10 seconds=30
println("  Time: $(minimum(bench_mf.times) / 1e6) ms")
println("  Allocations: $(bench_mf.allocs) allocs, $(bench_mf.memory / 1e6) MB")
println()

println("=" * 70)
println("Summary Statistics")
println("=" * 70)
println()

# Create summary table
st_time = minimum(bench_st.times) / 1e6
st_allocs = bench_st.allocs
st_mem = bench_st.memory / 1e6

mf_time = minimum(bench_mf.times) / 1e6
mf_allocs = bench_mf.allocs
mf_mem = bench_mf.memory / 1e6

println("SparseDiffTools (Colored Forward-Mode):")
println("  Min time: $st_time ms")
println("  Allocations: $st_allocs")
println("  Memory: $st_mem MB")
println()

try
    di_time = minimum(bench_di.times) / 1e6
    di_allocs = bench_di.allocs
    di_mem = bench_di.memory / 1e6
    
    println("DifferentiationInterface (Sparse Jacobian):")
    println("  Min time: $di_time ms")
    println("  Allocations: $di_allocs")
    println("  Memory: $di_mem MB")
    println()
    
    println("Speedup (ST vs DI):")
    println("  Time ratio: $(di_time / st_time)x")
    println()
catch
    println("DifferentiationInterface benchmark failed.")
    println()
end

println("SparseDiffTools (Matrix-Free):")
println("  Min time: $mf_time ms")
println("  Allocations: $mf_allocs")
println("  Memory: $mf_mem MB")
println()
