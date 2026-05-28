#=
Coloring algorithm comparison for a representative GXBeam sparse Jacobian.

Run with:

    julia --project=test test/coloring_comparison.jl
    julia --project=test test/coloring_comparison.jl --nelem=48

What it does:
- Builds a cantilever with `nelem` elements, runs `static_analysis(linear=true)`
  to get a representative sparse stiffness/Jacobian sparsity pattern.
- For each `SparseMatrixColorings` column-ordering, computes a greedy partial-
  distance-2 coloring and reports the number of colors used and the wall-time
  for the coloring step. Fewer colors = fewer columns of FD perturbations per
  Jacobian assembly = faster `autodiff_jacobian!`.
- Also benchmarks `GreedyColoringAlgorithm((ord1, ord2, ...))` which runs each
  ordering and keeps the one with the fewest colors — i.e. the "try all and
  pick best" mode.

This script uses GXBeam's re-exported `SparseMatrixColorings` to avoid bloating
`test/Project.toml`. No new test deps needed.
=#

using GXBeam
using LinearAlgebra
using SparseArrays
using Random
using Printf

const SMC = GXBeam.SparseMatrixColorings

# ---- args ----------------------------------------------------------------

const NELEM = let i = findfirst(s -> startswith(s, "--nelem="), ARGS)
    i === nothing ? 24 : parse(Int, ARGS[i][9:end])
end
const NSAMPLES = 5    # for the coloring time (it's deterministic-ish, this is
                      # mostly to dodge first-call noise)

# ---- representative problem ----------------------------------------------

function build_problem(nelem)
    L = 1.0
    xp = collect(range(0, L; length=nelem+1))
    yp = zero(xp); zp = zero(xp)
    points = [[xp[i], yp[i], zp[i]] for i in eachindex(xp)]
    # Full 6-DOF stiffness so the per-element coupling block is dense — gives
    # the coloring something nontrivial to do (vs. the EI-only diagonal that
    # the cantilever bending tests use).
    stiff = fill(Diagonal(fill(1e9, 6)), nelem)
    assembly = Assembly(points, 1:nelem, 2:nelem+1; stiffness=stiff)
    pc = Dict(nelem + 1 => PrescribedConditions(
        ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    distributed_loads = Dict(
        i => DistributedLoads(assembly, i; fz = s -> 1000.0) for i in 1:nelem)
    sys, _, _ = static_analysis(assembly;
        prescribed_conditions=pc,
        distributed_loads=distributed_loads,
        linear=true)
    return sys
end

# ---- driver --------------------------------------------------------------

sys = build_problem(NELEM)
J = sparse(sys.K)
n = size(J, 1)
println("Problem: nelem=$NELEM, n=$n, nnz=$(nnz(J)), ",
        "density=$(round(100 * nnz(J) / n^2; digits=2))%")
println()

const PROBLEM = SMC.ColoringProblem()

orderings = [
    ("NaturalOrder",         SMC.NaturalOrder()),
    ("LargestFirst",         SMC.LargestFirst()),
    ("SmallestLast",         SMC.SmallestLast()),
    ("IncidenceDegree",      SMC.IncidenceDegree()),
    ("DynamicLargestFirst",  SMC.DynamicLargestFirst()),
    ("RandomOrder(seed=0)",  SMC.RandomOrder(MersenneTwister(0), 0)),
]

"""
    bench_coloring(alg) -> (ncolors::Int, min_ms::Float64, mean_ms::Float64)

Run `coloring` once to warm, then `NSAMPLES` times and report both the min and
the mean wall-time in milliseconds. Min is the cleanest signal (closest to
"no scheduling noise"); mean tracks the cost the resolver would actually see
on average and surfaces outlier behavior like first-call JIT residue.
"""
function bench_coloring(alg)
    SMC.coloring(J, PROBLEM, alg)            # warm-up (JIT + first dispatch)
    samples_ns = Vector{UInt64}(undef, NSAMPLES)
    ncolors = 0
    for i in 1:NSAMPLES
        t0 = time_ns()
        result = SMC.coloring(J, PROBLEM, alg)
        t1 = time_ns()
        samples_ns[i] = t1 - t0
        ncolors = maximum(SMC.column_colors(result))
    end
    min_ms  = minimum(samples_ns) / 1e6
    mean_ms = (sum(samples_ns) / NSAMPLES) / 1e6
    return ncolors, min_ms, mean_ms
end

println("--- single ordering ---")
@printf("%-25s  %8s  %12s  %12s\n",
        "ordering", "ncolors", "min (ms)", "mean (ms)")
println("-" ^ 65)
single_results = Tuple{String, Int, Float64, Float64}[]
for (name, ord) in orderings
    alg = SMC.GreedyColoringAlgorithm(ord)
    nc, min_ms, mean_ms = bench_coloring(alg)
    @printf("%-25s  %8d  %12.3f  %12.3f\n", name, nc, min_ms, mean_ms)
    push!(single_results, (name, nc, min_ms, mean_ms))
end

# Pick the per-ordering best for comparison. Tie-break on min wall-time among
# rows that share the lowest color count.
fewest = minimum(r[2] for r in single_results)
best_single = argmin(r -> r[3], filter(r -> r[2] == fewest, single_results))
println()
println("Best single-ordering: $(best_single[1]) at $(best_single[2]) colors, ",
        "min=$(round(best_single[3]; digits=3)) ms / ",
        "mean=$(round(best_single[4]; digits=3)) ms")

println()
println("--- tuple-of-orderings (try all, keep fewest colors) ---")
@printf("%-25s  %8s  %12s  %12s\n",
        "ordering", "ncolors", "min (ms)", "mean (ms)")
println("-" ^ 65)
all_orderings = Tuple(o for (_, o) in orderings)
alg_tuple = SMC.GreedyColoringAlgorithm(all_orderings)
nc_tup, min_tup, mean_tup = bench_coloring(alg_tuple)
@printf("%-25s  %8d  %12.3f  %12.3f\n", "tuple(all six)",
        nc_tup, min_tup, mean_tup)

println()
if nc_tup < best_single[2]
    println("Tuple mode beats best single by ",
            best_single[2] - nc_tup, " color(s); pays ",
            @sprintf("%.2fx (min) / %.2fx (mean)",
                     min_tup / best_single[3],
                     mean_tup / best_single[4]),
            " more wall time.")
elseif nc_tup == best_single[2]
    println("Tuple mode matches best single ($(nc_tup) colors). Wall-time ratio ",
            @sprintf("%.2fx (min) / %.2fx (mean)",
                     min_tup / best_single[3],
                     mean_tup / best_single[4]),
            " vs the cheapest individual ordering.")
else
    println("Tuple mode landed at $(nc_tup) colors — unexpectedly worse than ",
            "the best single ordering. Investigate.")
end

println()
println("Note: coloring is a one-time prep cost per `System` (cached on ",
        "`prep_jacobian`). The metric that matters in production is ncolors, ",
        "because every Newton step does `ncolors` forward-diff evaluations to ",
        "assemble the Jacobian. Coloring wall-time only matters for cold-start ",
        "or for solves with a very small number of Newton steps.")
