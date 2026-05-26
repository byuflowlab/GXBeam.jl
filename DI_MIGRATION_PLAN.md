# GXBeam: Migrate AD Layer from SparseDiffTools to DifferentiationInterface

## Context

GXBeam currently uses SparseDiffTools (ST) for three AD operations: sparse Jacobian
coloring, colored Jacobian assembly, and matrix-free JVPs in the `xpfunc` GMRES solve
path. Switching to DifferentiationInterface (DI) consolidates these onto one modern,
backend-agnostic frontend and — critically — lets us cache the sparsity pattern and
coloring once and reuse it across every Newton step and every time step of a dynamic
solve, instead of redetecting per call as ST effectively forces.

Out of scope for this migration: replacing `ImplicitAD` with `ImplicitDifferentiation.jl`.
ImplicitAD keeps working unchanged as long as the new `drdy` return value satisfies its
contract.

This is a multi-session plan intended to be executed over several work blocks (including
overnight agent runs). Each session should leave the repo in a clean buildable state.

---

## Session workflow rules

These apply to every session that touches this migration:

1. **Branch hygiene.**
   - At session start: `git checkout jacobian-comparison && git pull` (if applicable),
     then `git checkout -B di-migration` (or fast-forward if already exists) so the
     migration branch is anchored at the latest `jacobian-comparison`.
   - All migration work happens on `di-migration`.
   - At session end (especially for overnight runs): commit any in-progress work on
     `di-migration`, then `git checkout jacobian-comparison` so Adam can work on other
     tasks during the day without the migration branch checked out.
2. **NEVER edit `Project.toml` or `Manifest.toml` directly.** These files are
   managed exclusively by Julia's `Pkg`. All dependency adds / removes / compat
   changes go through `Pkg.add`, `Pkg.rm`, `Pkg.update`, `Pkg.resolve`,
   `Pkg.compat` (or REPL `] add`/`rm`/etc.). User approval to "add a dep" is
   approval to run the Pkg command, not to hand-edit the file. Direct edits
   desync Project.toml from Manifest.toml and have caused incidents. Also: ask
   the user before running any `Pkg.*` command. See `.claude/rules/julia.md`.
3. **Run tests after each phase**, not at the end. The phases below are sized to be
   independently testable.
4. **No `--no-verify`, no force pushes, no amending shared commits.**

---

## Phase 0 — Baseline benchmark capture (no code changes)

**Goal:** Lock in pre-migration numbers so post-migration comparisons are credible.

1. On `jacobian-comparison`, run `julia --project=test test/jacobian_benchmark.jl` and
   tee the stdout to `benchmarks/baseline_sparsediff.txt`. Create `benchmarks/` if it
   doesn't exist.
2. Capture system info alongside it: Julia version, `Pkg.status()` output, CPU,
   `Threads.nthreads()`. Save to `benchmarks/baseline_env.txt`.
3. Commit both files on `jacobian-comparison` (separate commit from any migration work).

**Verification:** Files exist and contain the full markdown summary table from the
benchmark script.

---

## Phase 1 — Wrapper-level DI swap (no caching yet)

**Goal:** Make the three ST wrappers DI-backed while keeping their existing call
signatures, so downstream code in `analyses.jl` doesn't change yet.

**Files to modify:**
- [src/analyses.jl](src/analyses.jl) — the three wrappers at lines 3797-3818.
- [src/GXBeam.jl](src/GXBeam.jl) — line 17 import.

**Changes:**

1. **`jacobian_colors`** ([analyses.jl:3797](src/analyses.jl#L3797))
   - Replace `SparseDiffTools.matrix_colors(jacob)` with
     `SparseMatrixColorings.coloring(...)` using `GreedyColoringAlgorithm()` against the
     pattern derived from the |J1|+|J2|+|J3| sum (same trick currently used).
   - Return the colors vector in the same shape consumers expect (a `Vector{Int}` of
     column colors).

2. **`autodiff_jacobian!`** ([analyses.jl:3808](src/analyses.jl#L3808))
   - Replace `SparseDiffTools.forwarddiff_color_jacobian!` with a DI-based
     implementation. For this phase, build a fresh `AutoSparse(AutoForwardDiff(); ...)`
     backend and call `prepare_jacobian` + `jacobian!` per call. Phase 3 will cache the
     prep; this phase prioritizes correctness over speed.
   - Keep the `colors` keyword for API stability; if passed, it can be ignored (DI
     handles coloring internally) or used to seed a `PrecomputedColoring` if the
     `SparseMatrixColorings` API allows it.
   - Call sites: [analyses.jl:922](src/analyses.jl#L922), [:925](src/analyses.jl#L925),
     [:1363](src/analyses.jl#L1363), [:1366](src/analyses.jl#L1366) — these should all
     continue to work without modification.

3. **`matrixfree_jacobian`** ([analyses.jl:3816](src/analyses.jl#L3816))
   - Return a `LinearMap{Float64}(n, n; ismutating=true)` that wraps a prepared DI
     `pushforward`. The pattern from `newton_solve_di_gmres!` in
     [test/jacobian_benchmark.jl:252-278](test/jacobian_benchmark.jl#L252-L278) is the
     reference.
   - **Critical: must also support the adjoint action** (`mul!(out, A', v)`).
     `ImplicitAD.implicit` uses `lsolve(A', ybar)` in the reverse path. Wrap a prepared
     `pullback` for the adjoint matvec. `LinearMaps.jl` supports this via the `fc`
     argument or `LinearMap(f, fc, n)`.
   - Verify ImplicitAD's default `lsolve` works with a `LinearMap` (it does GMRES under
     the hood; if it tries `A \ b` directly, that will fail and we'll need to pass a
     custom `lsolve`). Test this with one xpfunc call before declaring the phase done.
   - Call sites: all 7 `ImplicitAD.implicit(...; drdy=matrixfree_jacobian)` lines plus
     the `coupled_jacobian=matrixfree_jacobian` defaults in `matrixfree_lsolve!`
     ([:3527](src/analyses.jl#L3527)) and `matrixfree_nlsolve!` ([:3669](src/analyses.jl#L3669)).
     Also the line-search calls inside `matrixfree_nlsolve!` at
     [:3724](src/analyses.jl#L3724), [:3732](src/analyses.jl#L3732) that do
     `mul!(xtmp, jacob, dx)`.

4. **`GXBeam.jl` imports**: replace `import SparseDiffTools` with
   `import DifferentiationInterface as DI`, `import SparseMatrixColorings`,
   `import LinearMaps`. Adding the corresponding deps to `Project.toml` is **not**
   a text edit — use `Pkg.add` via Julia (see workflow rule 2 at the top of this
   file). Each Pkg command requires explicit user approval.

**Verification:**
- Run `julia --project test/runtests.jl` — all tests pass.
- Re-run `test/jacobian_benchmark.jl`, save output to `benchmarks/phase1_di_naive.txt`.
- Compare GMRES Newton solve times to Phase 0 baseline. Expect roughly comparable
  numbers (Phase 1 doesn't yet exploit prep caching).

---

## Phase 2 — Dependency cleanup

**Goal:** Remove SparseDiffTools from `Project.toml` once nothing imports it.

**Requires explicit user approval before running any `Pkg` commands.**

1. Verify no `SparseDiffTools.` references remain in `src/` (grep).
2. **Do not edit Project.toml by hand.** Use Pkg:
   `julia --project -e 'using Pkg; Pkg.rm("SparseDiffTools"); Pkg.add(["DifferentiationInterface", "SparseMatrixColorings"])'`
   (LinearMaps is already a dep — confirm before adding). Pin compat with
   `Pkg.compat("DifferentiationInterface", "0.7")` etc. Each Pkg command needs
   explicit user approval.
3. Inspect the resulting Project.toml + Manifest.toml diff with `git diff` to
   confirm Pkg made the expected edits.
4. Run full test suite.

**Verification:** Tests pass; `Pkg.status()` shows no SparseDiffTools.

---

## Phase 3 — Cache DI prep on System objects

**Goal:** Pay sparsity detection + coloring once per `System` lifetime instead of per
call. This is the main performance win, especially for time-domain solves where the
same residual structure is re-Jacobian'd thousands of times.

**Files to modify:**
- [src/system.jl](src/system.jl) — add cache fields to `StaticSystem` (line 215),
  `DynamicSystem` (line 267), `ExpandedSystem` (line 323).
- [src/analyses.jl](src/analyses.jl) — thread the cached prep through `autodiff_jacobian!`
  and `matrixfree_jacobian`, and through the `constants` namedtuple where appropriate.

**Design choices to make during implementation:**

- **Where the cache lives.** Options:
  (a) Add a `prep_jacobian::Any` field to each System struct (loses type stability but
      simple). (b) Make System parametric on a prep type (type-stable but invasive).
      (c) Store prep on the `constants` namedtuple at solve setup time (no struct
      changes; cache is per-solve, not per-System). **Lean toward (c) for the first
      cut** — minimum surface area, and the lifetime of `constants` already spans the
      whole solve including all Newton steps.
- **One prep per residual variant.** Static, steady, expanded-steady, initial, and
  newmark each have their own residual function — they each need their own prep. The
  cache should be a small named-tuple of preps keyed by residual variant, built lazily
  on first use within a solve.
- **Pushforward/pullback prep for the matrix-free path** is similarly per-residual.

**Verification:**
- Tests pass.
- Re-run benchmark, save to `benchmarks/phase3_di_cached.txt`. Compare:
  - GMRES Newton solve (DI) vs Phase 1 — expect speedup on multi-Newton-step solves.
  - Compute-only Jacobian (DI) numbers stay the same (already pre-prepped in benchmark).
  - Look at the dynamic / newmark examples in `test/examples/` for the biggest wins.

---

## Phase 4 — Final comparison + writeup

**Goal:** Produce a self-contained comparison document.

1. Re-run the full benchmark suite on `di-migration` HEAD. Save to
   `benchmarks/final_di.txt`.
2. Add a short markdown summary `benchmarks/MIGRATION_SUMMARY.md` with side-by-side
   tables: Phase 0 baseline vs Phase 1 naive DI vs Phase 3 cached DI. Include the
   absolute times *and* the ratios.
3. Run the example notebooks / scripts in `test/examples/` end-to-end to confirm no
   regressions on real workloads (especially the rotating beam and joined-wing cases
   that exercise time-domain solves heavily).
4. Open a PR from `di-migration` → `jacobian-comparison` (or `master`, user's choice)
   with the summary file as the PR body.

---

## Critical files (quick reference)

- [src/analyses.jl:3797-3818](src/analyses.jl#L3797-L3818) — the three ST wrappers
- [src/analyses.jl:922,925,1363,1366](src/analyses.jl#L922) — `autodiff_jacobian!` call sites
- [src/analyses.jl:171,505,511,1910,2164,2768,3116](src/analyses.jl#L171) — `ImplicitAD.implicit(...; drdy=...)` sites
- [src/analyses.jl:3527,3669](src/analyses.jl#L3527) — `matrixfree_lsolve!` / `matrixfree_nlsolve!` defaults
- [src/system.jl:215,267,323](src/system.jl#L215) — System struct definitions
- [src/GXBeam.jl:17](src/GXBeam.jl#L17) — `import SparseDiffTools`
- [Project.toml:24,46](Project.toml#L24) — SparseDiffTools dep + compat
- [test/jacobian_benchmark.jl](test/jacobian_benchmark.jl) — benchmark harness (already
  exercises both backends; the `newton_solve_di_gmres!` pattern is the reference for
  the new `matrixfree_jacobian`)
- [test/runtests.jl](test/runtests.jl) — primary test entry point
- [test/jacobians.jl](test/jacobians.jl) — Jacobian-specific tests

## Reusable code already in the codebase

- `ResidualInPlace` functor pattern from
  [test/jacobian_benchmark.jl:94-100](test/jacobian_benchmark.jl#L94-L100) —
  use the same shape for cached DI residual closures (DI requires type-identical
  function objects between prep and exec).
- `newton_solve_di_gmres!`
  [test/jacobian_benchmark.jl:252-278](test/jacobian_benchmark.jl#L252-L278) —
  reference implementation for `LinearMap`-wrapped pushforward, including the
  `v_buf`/`dy_buf` buffer trick that handles GMRES's SubArray views.
- `prepare_jacobian` + `AutoSparse(AutoForwardDiff(); ...)` setup
  [test/jacobian_benchmark.jl:110-116](test/jacobian_benchmark.jl#L110-L116).

---

## Confirmed decisions

- **Plan file location:** `DI_MIGRATION_PLAN.md` at the repo root.
- **Migration branch name:** `di-migration`, branched from `jacobian-comparison`.
- **Cut strategy:** Hard cut. No flag-gated dual implementation.
