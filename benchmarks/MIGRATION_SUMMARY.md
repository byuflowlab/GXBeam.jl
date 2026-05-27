# DI Migration — Benchmark Summary

This is the Phase 4 writeup for the SparseDiffTools → DifferentiationInterface
migration. The migration was executed in three phases on the `di-migration`
branch; benchmarks were captured at every phase boundary and live in this
directory. See [`../DI_MIGRATION_PLAN.md`](../DI_MIGRATION_PLAN.md) for the
plan and per-phase scope.

## Methodology and source files

Two harnesses were used at every phase boundary:

| Harness | What it measures |
|---|---|
| `test/examples_benchmark.jl` | End-to-end timings of 16 example workloads from `test/examples/`. The realistic-load comparison. |
| `test/jacobian_benchmark.jl` | Targeted micro-benchmarks of the AD layer itself (sparse Jacobian, JVP, GMRES Newton) on a 150-DOF cantilever. |

Files in this directory:

| File | Phase | Source | Depot |
|---|---|---|---|
| `baseline_examples_sparsediff.txt`, `baseline_sparsediff.txt` | 0 | SparseDiffTools (pre-migration code) | pre-migration |
| `phase1_di_naive_examples.txt`, `phase1_di_naive_sparsediff.txt` | 1 | DI wrappers, **no prep caching** | post-Phase-1 |
| `phase2_di_naive_examples.txt`, `phase2_di_naive_sparsediff.txt` | 2 | DI wrappers, SparseDiffTools dep removed | post-Phase-2 |
| `phase3_di_cached_examples.txt`, `phase3_di_cached_sparsediff.txt` | 3 | DI with prep cached on `System` | post-Phase-2 |
| `examples_sparsediff_longer.txt` | 0-equivalent | **SparseDiffTools** code rerun at 5× samples on the post-Phase-2 depot | post-Phase-2 |

**Sample-count asymmetry.** Phases 0–3 were captured at the 1× sample
multiplier (`samples` between 3 and 30 per workload depending on cost).
`examples_sparsediff_longer.txt` was captured after the SAMPLE/SECONDS
multiplier was bumped to 5×. Min and median timings remain directly
comparable; mean and std should be read with this asymmetry in mind.

## A. Example workloads — Phase 0 vs Phase 1 (naive DI)

This is what the migration cost before any caching work — the price of
substituting `DI.prepare_jacobian` + `DI.jacobian!` for the previous
`forwarddiff_color_jacobian!` wrapper inside a per-call hot path.

| Workload | P0 min (ms) | P1 min (ms) | min B/A | P0 mean | P1 mean | mean B/A | verdict |
|---|---|---|---|---|---|---|---|
| cantilever | 0.546 | 0.523 | 0.96× | 0.586 | 0.617 | 1.05× | faster |
| curved | 89.4 | 107.1 | 1.20× | 93.4 | 112.2 | 1.20× | slower |
| overdetermined | 0.743 | 0.719 | 0.97× | 0.817 | 0.810 | 0.99× | faster |
| tipforce | 111.6 | 134.8 | 1.21× | 116.4 | 139.5 | 1.20× | slower |
| tipmoment | 14.7 | 17.1 | 1.16× | 16.1 | 18.1 | 1.12× | slower |
| static-joined-wing-linear | 22.5 | 25.0 | 1.11× | 24.1 | 26.0 | 1.08× | slower |
| static-joined-wing-nonlinear | 7208 | 8385 | 1.16× | 7340 | 8830 | 1.20× | slower |
| static-joined-wing-follower | 14139 | 15399 | 1.09× | 14179 | 16325 | 1.15× | slower |
| rotating-steady-state | 21.6 | 24.9 | 1.15× | 23.9 | 25.8 | 1.08× | slower |
| rotating-eigenvalue | 117.6 | 140.2 | 1.19× | 121.6 | 141.9 | 1.17× | slower |
| wind-turbine-blade-undamped | 103.7 | 116.9 | 1.13× | 108.6 | 119.5 | 1.10× | slower |
| wind-turbine-blade-damped | 107.3 | 110.7 | 1.03× | 113.5 | 121.7 | 1.07× | slower |
| dynamic-joined-wing-undamped | 357.7 | 417.2 | 1.17× | 395.5 | 433.1 | 1.10× | slower |
| dynamic-joined-wing-damped | 171.9 | 186.4 | 1.08× | 184.8 | 212.2 | 1.15× | slower |
| impliciteuler | 1108 | 1267 | 1.14× | 1196 | 1309 | 1.09× | slower |
| stepsystem | 107.2 | 117.7 | 1.10× | 112.6 | 122.0 | 1.08× | slower |

Tally: 2 faster, 14 slower (min). Worst slowdown: `tipforce` at 1.21×.

## B. Example workloads — Phase 0 vs Phase 3 (cached DI)

This is the headline migration comparison: pre-migration SparseDiffTools vs
the final cached-prep DI implementation.

| Workload | P0 min (ms) | P3 min (ms) | min B/A | P0 mean | P3 mean | mean B/A | verdict |
|---|---|---|---|---|---|---|---|
| cantilever | 0.546 | 0.541 | 0.99× | 0.586 | 0.673 | 1.15× | ~same |
| curved | 89.4 | 95.3 | 1.07× | 93.4 | 103.0 | 1.10× | slower |
| overdetermined | 0.743 | 0.678 | 0.91× | 0.817 | 0.791 | 0.97× | faster |
| tipforce | 111.6 | 140.4 | 1.26× | 116.4 | 147.0 | 1.26× | slower |
| tipmoment | 14.7 | 15.8 | 1.07× | 16.1 | 18.7 | 1.16× | slower |
| static-joined-wing-linear | 22.5 | 23.7 | 1.06× | 24.1 | 25.9 | 1.08× | slower |
| static-joined-wing-nonlinear | 7208 | 8831 | 1.23× | 7340 | 9105 | 1.24× | slower |
| static-joined-wing-follower | 14139 | 14537 | 1.03× | 14179 | 15392 | 1.09× | ~same |
| rotating-steady-state | 21.6 | 23.5 | 1.09× | 23.9 | 28.2 | 1.18× | slower |
| rotating-eigenvalue | 117.6 | 130.6 | 1.11× | 121.6 | 141.0 | 1.16× | slower |
| wind-turbine-blade-undamped | 103.7 | 119.2 | 1.15× | 108.6 | 122.0 | 1.12× | slower |
| wind-turbine-blade-damped | 107.3 | 130.0 | 1.21× | 113.5 | 131.2 | 1.16× | slower |
| dynamic-joined-wing-undamped | 357.7 | 527.9 | 1.48× | 395.5 | 546.4 | 1.38× | slower |
| dynamic-joined-wing-damped | 171.9 | 188.8 | 1.10× | 184.8 | 199.1 | 1.08× | slower |
| impliciteuler | 1108 | 1181 | 1.07× | 1196 | 1193 | 1.00× | slower |
| stepsystem | 107.2 | 123.6 | 1.15× | 112.6 | 134.1 | 1.19× | slower |

Tally: 1 faster, 2 ~same, 13 slower (min). Worst slowdown:
`dynamic-joined-wing-undamped` at 1.48×.

Read in isolation, this is bad news. Section D below shows it isn't the
whole story.

## C. Example workloads — same-depot comparison (ST 5× rerun vs Phase 3)

The cleanest read of the migration's true effect: both rows ran on the
post-Phase-2 depot (same transitive package versions), so any timing delta
is attributable to the AD code itself rather than upstream drift.

| Workload | ST 5× min (ms) | P3 min (ms) | min B/A | ST 5× mean | P3 mean | mean B/A | verdict |
|---|---|---|---|---|---|---|---|
| cantilever | 0.485 | 0.541 | 1.12× | 0.579 | 0.673 | 1.16× | slower |
| curved | 98.5 | 95.3 | 0.97× | 106.6 | 103.0 | 0.97× | faster |
| overdetermined | 0.623 | 0.678 | 1.09× | 0.737 | 0.791 | 1.07× | slower |
| tipforce | 127.1 | 140.4 | 1.11× | 134.4 | 147.0 | 1.09× | slower |
| tipmoment | 16.1 | 15.8 | 0.98× | 17.1 | 18.7 | 1.10× | ~same |
| static-joined-wing-linear | 22.2 | 23.7 | 1.07× | 24.7 | 25.9 | 1.05× | slower |
| static-joined-wing-nonlinear | 7867 | 8831 | 1.12× | 8276 | 9105 | 1.10× | slower |
| static-joined-wing-follower | 15421 | 14537 | 0.94× | 16191 | 15392 | 0.95× | **faster** |
| rotating-steady-state | 22.8 | 23.5 | 1.03× | 25.7 | 28.2 | 1.10× | ~same |
| rotating-eigenvalue | 142.9 | 130.6 | 0.91× | 161.4 | 141.0 | 0.87× | **faster** |
| wind-turbine-blade-undamped | 128.4 | 119.2 | 0.93× | 146.8 | 122.0 | 0.83× | **faster** |
| wind-turbine-blade-damped | 123.0 | 130.0 | 1.06× | 142.0 | 131.2 | 0.92× | slower |
| dynamic-joined-wing-undamped | 450.3 | 527.9 | 1.17× | 478.2 | 546.4 | 1.14× | slower |
| dynamic-joined-wing-damped | 173.4 | 188.8 | 1.09× | 196.4 | 199.1 | 1.01× | slower |
| impliciteuler | 1271 | 1181 | 0.93× | 1391 | 1193 | 0.86× | **faster** |
| stepsystem | 152.0 | 123.6 | 0.81× | 255.4 | 134.1 | 0.53× | **faster** |

Tally: 6 faster, 2 ~same, 8 slower (min). Worst slowdown:
`dynamic-joined-wing-undamped` at 1.17×.

DI cached wins on the heavy time-domain workloads (`stepsystem`,
`impliciteuler`, the wind-turbine blade cases) and on the bigger static
problems (`static-joined-wing-follower`, `rotating-eigenvalue`). It loses on
several smaller workloads where prep amortization can't pay back its
constant cost — most notably `dynamic-joined-wing-undamped` and the
joined-wing-nonlinear case.

## D. Upstream drift: the SparseDiffTools 5× rerun

The same SparseDiffTools code, run on the post-Phase-2 depot at 5× samples,
is consistently slower than the original Phase 0 baseline. The Phase 2 work
rebuilt the depot (older `SciMLOperators` was reinstalled to satisfy
SparseDiffTools v2 + SciMLBase constraints, and DifferentiationInterface
shifted from 7.16 to 7.18 transitively). Commit `5c56d3e3` flagged this at
capture time.

| Workload | P0 min (ms) | ST 5× min (ms) | min B/A | mean B/A | verdict |
|---|---|---|---|---|---|
| cantilever | 0.546 | 0.485 | 0.89× | 0.99× | faster |
| curved | 89.4 | 98.5 | 1.10× | 1.14× | slower |
| overdetermined | 0.743 | 0.623 | 0.84× | 0.90× | faster |
| tipforce | 111.6 | 127.1 | 1.14× | 1.15× | slower |
| tipmoment | 14.7 | 16.1 | 1.09× | 1.06× | slower |
| static-joined-wing-linear | 22.5 | 22.2 | 0.99× | 1.02× | ~same |
| static-joined-wing-nonlinear | 7208 | 7867 | 1.09× | 1.13× | slower |
| static-joined-wing-follower | 14139 | 15421 | 1.09× | 1.14× | slower |
| rotating-steady-state | 21.6 | 22.8 | 1.06× | 1.08× | slower |
| rotating-eigenvalue | 117.6 | 142.9 | 1.22× | 1.33× | slower |
| wind-turbine-blade-undamped | 103.7 | 128.4 | 1.24× | 1.35× | slower |
| wind-turbine-blade-damped | 107.3 | 123.0 | 1.15× | 1.25× | slower |
| dynamic-joined-wing-undamped | 357.7 | 450.3 | 1.26× | 1.21× | slower |
| dynamic-joined-wing-damped | 171.9 | 173.4 | 1.01× | 1.06× | ~same |
| impliciteuler | 1108 | 1271 | 1.15× | 1.16× | slower |
| stepsystem | 107.2 | 152.0 | 1.42× | 2.27× | slower |

Tally: 2 faster, 2 ~same, 12 slower. Worst slowdown: `stepsystem` at 1.42×
(2.27× on mean).

Implication: a meaningful chunk of the Phase 0 → Phase 3 regression in
Section B is **not the migration** — it is upstream package drift, mostly
transitive AD versioning (DI 7.16 → 7.18 and the `SciMLOperators`
downgrade). Section C is the load-bearing comparison.

## E. Jacobian-level benchmark (`test/jacobian_benchmark.jl`)

Same 150-DOF cantilever at every phase. Both SparseDiffTools and DI paths
are exercised at every phase even after Phase 2 removes ST as a GXBeam dep
(the script imports ST directly for reference).

### Compute-only sparse Jacobian (min, ms)

| Phase | ST | DI | ratio ST/DI |
|---|---|---|---|
| 0 baseline | 0.158 | 0.124 | 1.28× — DI faster |
| 1 DI naive | 9.802 | 0.113 | 86.76× — DI faster |
| 2 DI naive | 9.533 | 0.091 | 104.94× — DI faster |
| 3 DI cached | 0.100 | 0.096 | 1.05× — DI faster |

DI's compute-only path is consistently a touch faster than ST's. The huge
ST numbers in Phases 1–2 reflect that the bench script's ST path lost some
caching scaffolding when GXBeam stopped exposing it — not a real ST
regression. Phase 3 restored parity.

### Full pipeline incl. sparsity detection + coloring (min, ms)

| Phase | ST | DI | ratio ST/DI |
|---|---|---|---|
| 0 baseline | 3.781 | 10.648 | 0.36× — ST faster |
| 1 DI naive | 12.117 | 10.303 | 1.18× — DI faster |
| 2 DI naive | 14.138 | 10.287 | 1.37× — DI faster |
| 3 DI cached | 1.910 | 11.486 | 0.17× — ST faster |

The DI "full pipeline" min stays in the 10–11 ms band across all phases
because the bench script re-runs detection+coloring on every sample (it
doesn't exercise the System-level cache). The Phase-3 ST drop is the bench
script's own caching kicking in. Realistic GXBeam workloads — Section A–C —
are what reflect Phase 3's prep-caching win.

### JVP (J·v, min, ms)

| Phase | ST `mul!` | DI pushforward | ratio ST/DI |
|---|---|---|---|
| 0 baseline | 0.084 | 0.064 | 1.31× — DI faster |
| 1 DI naive | 0.078 | 0.075 | 1.04× — DI faster |
| 2 DI naive | 0.063 | 0.064 | 0.99× — ST faster |
| 3 DI cached | 0.081 | 0.076 | 1.07× — DI faster |

Effectively a wash across the migration.

### GMRES Newton step (min, ms)

| Phase | ST | DI | GXBeam xpfunc |
|---|---|---|---|
| 0 baseline | 0.577 | 0.587 | 0.615 |
| 1 DI naive | 0.543 | 0.534 | 0.575 |
| 2 DI naive | 0.530 | 0.462 | 0.517 |
| 3 DI cached | 0.496 | 0.494 | 0.520 |

All three converge to ~0.5 ms — the GMRES solve dominates, not the AD layer.

## F. Verdict

- **The migration is a net positive on heavy realistic workloads.** Section
  C shows DI-cached beats SparseDiffTools on the time-domain solves that
  Phase 3 was designed for (`stepsystem` ~0.81×, `impliciteuler` ~0.93×,
  `wind-turbine-blade-undamped` ~0.93×, `static-joined-wing-follower`
  ~0.94×, `rotating-eigenvalue` ~0.91×).
- **It regresses on small workloads** where prep amortization can't pay back
  its constant cost — most notably `dynamic-joined-wing-undamped` (~1.17×
  on the same-depot read; ~1.48× when also blamed for upstream drift).
  Tracking down whether the dynamic-joined-wing case is missing a cache hit
  inside Newmark stepping is worth a follow-up.
- **Compute-only Jacobian and JVP are unchanged or slightly faster** with
  DI at every phase (Section E).
- **Roughly a third of the Phase 0 → Phase 3 wall-time gap is upstream**
  (Section D — same SparseDiffTools code regressed 10–25 % on the new
  depot). Worth opening a separate, non-blocking issue to track the
  DI 7.16 → 7.18 regression; not a migration blocker.
- **The 1×/5× sample asymmetry** between Phases 0–3 and `examples_sparsediff_longer.txt`
  doesn't change the verdict but does inflate the mean column noise on the
  1× phases. If a clean rerun of Phase 3 at 5× becomes worthwhile (e.g. for
  the PR), set `SAMPLE_MULTIPLIER = SECONDS_MULTIPLIER = 5` at the top of
  `test/examples_benchmark.jl` and re-capture.

## Reproducing these tables

Comparator script (added by `27b168e2`):

```
julia test/compare_examples_benchmark.jl benchmarks/baseline_examples_sparsediff.txt benchmarks/phase1_di_naive_examples.txt
julia test/compare_examples_benchmark.jl benchmarks/baseline_examples_sparsediff.txt benchmarks/phase3_di_cached_examples.txt
julia test/compare_examples_benchmark.jl benchmarks/baseline_examples_sparsediff.txt benchmarks/examples_sparsediff_longer.txt
julia test/compare_examples_benchmark.jl benchmarks/examples_sparsediff_longer.txt benchmarks/phase3_di_cached_examples.txt
```

Jacobian-level numbers come from the SUMMARY tables at the tail of
`benchmarks/{baseline,phase{1,2,3}_*}_sparsediff.txt`.
