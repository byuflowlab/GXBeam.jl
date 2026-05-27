#!/usr/bin/env julia
"""
    compare_examples_benchmark.jl <file_a> <file_b>

Parses the markdown SUMMARY tables emitted by `test/examples_benchmark.jl`
from two output files and prints a per-workload comparison: min time A,
min time B, ratio (B/A), and a tag (faster / slower / ~same).

Defaults to comparing benchmarks/baseline_examples_sparsediff.txt against
benchmarks/examples_sparsediff_longer.txt. Pass any two paths to override.

Usage:
    julia test/compare_examples_benchmark.jl
    julia test/compare_examples_benchmark.jl <file_a> <file_b>
    julia test/compare_examples_benchmark.jl <file_a> <file_b> --col=mean
"""

using Printf

# ---- argument parsing -----------------------------------------------------

const COL_LABELS = Dict(
    "min" => "Min (ms)", "mean" => "Mean (ms)", "std" => "Std (ms)",
    "med" => "Med (ms)", "max" => "Max (ms)",
)
const COL_INDEX = Dict("min" => 1, "mean" => 2, "std" => 3, "med" => 4, "max" => 5)

function parse_args(argv)
    files = String[]
    col = "min"
    for arg in argv
        if startswith(arg, "--col=")
            col = arg[7:end]
        elseif startswith(arg, "--")
            error("unknown flag: $arg")
        else
            push!(files, arg)
        end
    end
    if isempty(files)
        files = ["benchmarks/baseline_examples_sparsediff.txt",
                "benchmarks/examples_sparsediff_longer.txt"]
    elseif length(files) != 2
        error("expected 2 positional file args, got $(length(files))")
    end
    haskey(COL_INDEX, col) || error("--col must be one of $(keys(COL_INDEX))")
    return (file_a=files[1], file_b=files[2], col=col)
end

# ---- table parsing --------------------------------------------------------

"""
Parse the SUMMARY table out of an examples_benchmark.jl output file.
Returns a Vector of (workload_name, [min, mean, std, med, max]) tuples,
preserving original ordering. Allocs and Mem are dropped.
"""
function parse_summary(path::AbstractString)
    rows = Tuple{String, Vector{Float64}}[]
    in_table = false
    for line in eachline(path)
        line = strip(line)
        # Look for the table header that marks the start.
        if startswith(line, "| Workload |")
            in_table = true
            continue
        end
        # Skip the markdown separator row.
        if in_table && startswith(line, "|---")
            continue
        end
        # Data row: starts with "| " followed by a workload name.
        if in_table && startswith(line, "| ")
            parts = split(line, "|"; keepempty=false)
            # parts: [" workload ", " min ", " mean ", " std ", " med ", " max ", " allocs ", " mem "]
            length(parts) >= 6 || continue
            name = strip(parts[1])
            # The first five numeric columns are the timings.
            timings = Float64[]
            for i in 2:6
                try
                    push!(timings, parse(Float64, strip(parts[i])))
                catch
                    # Non-numeric (e.g. blank line below table) — stop parsing.
                    return rows
                end
            end
            push!(rows, (name, timings))
        elseif in_table && isempty(line)
            # Blank line ends the table.
            return rows
        end
    end
    return rows
end

# ---- comparison + report --------------------------------------------------

function classify(ratio::Float64; threshold::Float64=0.03)
    if abs(ratio - 1) <= threshold
        return "~same"
    elseif ratio < 1
        return "faster"
    else
        return "slower"
    end
end

function compare(a_rows, b_rows; verdict_col::String="min")
    verdict_idx = COL_INDEX[verdict_col]
    min_idx, mean_idx = COL_INDEX["min"], COL_INDEX["mean"]
    a_dict = Dict(name => t for (name, t) in a_rows)
    b_dict = Dict(name => t for (name, t) in b_rows)

    # Use union, preserving insertion order from a then b.
    names = String[]
    seen = Set{String}()
    for (n, _) in a_rows
        push!(seen, n); push!(names, n)
    end
    for (n, _) in b_rows
        n in seen || (push!(seen, n); push!(names, n))
    end

    rows = NamedTuple[]
    for n in names
        a = get(a_dict, n, nothing)
        b = get(b_dict, n, nothing)
        a_min  = isnothing(a) ? missing : a[min_idx]
        b_min  = isnothing(b) ? missing : b[min_idx]
        a_mean = isnothing(a) ? missing : a[mean_idx]
        b_mean = isnothing(b) ? missing : b[mean_idx]
        a_v = isnothing(a) ? missing : a[verdict_idx]
        b_v = isnothing(b) ? missing : b[verdict_idx]
        r_min  = (ismissing(a_min)  || ismissing(b_min))  ? missing : (b_min  / a_min)
        r_mean = (ismissing(a_mean) || ismissing(b_mean)) ? missing : (b_mean / a_mean)
        r_v    = (ismissing(a_v)    || ismissing(b_v))    ? missing : (b_v    / a_v)
        tag = ismissing(r_v) ? "missing" : classify(r_v)
        push!(rows, (; name=n, a_min, b_min, r_min, a_mean, b_mean, r_mean, tag))
    end
    return rows
end

function print_report(file_a, file_b, verdict_col, rows)
    println()
    println("Comparing timings (B vs A)")
    println("  A (baseline):  $file_a")
    println("  B (rerun):     $file_b")
    println("  ratio = B / A;  >1 means B slower, <1 means B faster")
    println("  verdict uses the $verdict_col column")
    println()
    # Header: workload, min A, min B, min B/A, mean A, mean B, mean B/A, verdict
    @printf("%-32s %10s %10s %8s   %10s %10s %8s   %s\n",
            "Workload", "MinA(ms)", "MinB(ms)", "min B/A",
            "MeanA(ms)", "MeanB(ms)", "mean B/A", "verdict")
    println(repeat('-', 110))
    n_slower = 0; n_faster = 0; n_same = 0; n_missing = 0
    worst_name = ""; worst_ratio = -Inf
    for r in rows
        a_min  = ismissing(r.a_min)  ? "    -"   : @sprintf("%10.3f", r.a_min)
        b_min  = ismissing(r.b_min)  ? "    -"   : @sprintf("%10.3f", r.b_min)
        rmin   = ismissing(r.r_min)  ? "    -"   : @sprintf("%8.2fx", r.r_min)
        a_mean = ismissing(r.a_mean) ? "    -"   : @sprintf("%10.3f", r.a_mean)
        b_mean = ismissing(r.b_mean) ? "    -"   : @sprintf("%10.3f", r.b_mean)
        rmean  = ismissing(r.r_mean) ? "    -"   : @sprintf("%8.2fx", r.r_mean)
        @printf("%-32s %s %s %s   %s %s %s   %s\n",
                r.name, a_min, b_min, rmin, a_mean, b_mean, rmean, r.tag)
        if r.tag == "slower"
            n_slower += 1
            verdict_r = verdict_col == "mean" ? r.r_mean : r.r_min
            if !ismissing(verdict_r) && verdict_r > worst_ratio
                worst_ratio = verdict_r; worst_name = r.name
            end
        elseif r.tag == "faster"; n_faster += 1
        elseif r.tag == "~same"; n_same += 1
        else; n_missing += 1
        end
    end
    println()
    println("Summary: $n_faster faster, $n_same ~same, $n_slower slower, $n_missing missing")
    if n_slower > 0 && !isempty(worst_name)
        @printf("Worst slowdown (%s): %s at %.2fx\n", verdict_col, worst_name, worst_ratio)
    end
    return n_slower
end

# ---- main -----------------------------------------------------------------

opts = parse_args(ARGS)
a_rows = parse_summary(opts.file_a)
b_rows = parse_summary(opts.file_b)
isempty(a_rows) && error("could not parse a SUMMARY table from $(opts.file_a)")
isempty(b_rows) && error("could not parse a SUMMARY table from $(opts.file_b)")
rows = compare(a_rows, b_rows; verdict_col=opts.col)
print_report(opts.file_a, opts.file_b, opts.col, rows)
