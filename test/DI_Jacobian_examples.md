======================================================================
JACOBIAN BENCHMARK: SparseDiffTools vs DifferentiationInterface
======================================================================

--- Compute-only (sparsity/coloring pre-computed) ---

SparseDiffTools (colored ForwardDiff, compute only)...
BenchmarkTools.Trial: 30 samples with 1 evaluation per sample.
 Range (min … max):  108.677 μs … 283.227 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     135.328 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   143.186 μs ±  34.373 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▂       ▂█                                                     
  █▅█▁▅▁▅█████▅▁▅▁▁▁▁▁█▁▁▁▁▅▁▁▅▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅ ▁
  109 μs           Histogram: frequency by time          283 μs <

 Memory estimate: 337.38 KiB, allocs estimate: 1322.

DifferentiationInterface (sparse ForwardDiff, compute only)...
BenchmarkTools.Trial: 30 samples with 1 evaluation per sample.
 Range (min … max):  101.214 μs … 157.379 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     105.363 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   110.301 μs ±  13.059 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▁▁█                                                            
  ███▇▇▇▄▇▁▁▇▁▁▇▄▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▄ ▁
  101 μs           Histogram: frequency by time          157 μs <

 Memory estimate: 344.44 KiB, allocs estimate: 1324.

--- Full pipeline (includes sparsity detection + coloring) ---

SparseDiffTools (full pipeline)...
BenchmarkTools.Trial: 10 samples with 1 evaluation per sample.
 Range (min … max):  2.104 ms … 9.256 ms  ┊ GC (min … max):  0.00% … 75.05%
 Time  (median):     2.393 ms             ┊ GC (median):     0.00%
 Time  (mean ± σ):   3.927 ms ± 2.904 ms  ┊ GC (mean ± σ):  35.58% ± 31.83%

  █▃                                                      ▃  
  ██▇▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  2.1 ms         Histogram: frequency by time       9.26 ms <

 Memory estimate: 8.26 MiB, allocs estimate: 18810.

DifferentiationInterface (full pipeline)...
BenchmarkTools.Trial: 10 samples with 1 evaluation per sample.
 Range (min … max):  13.350 ms … 28.873 ms  ┊ GC (min … max):  0.00% … 30.49%
 Time  (median):     18.878 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   19.129 ms ±  4.844 ms  ┊ GC (mean ± σ):  12.18% ± 15.25%

  █       ▁  ▁        ▁▁     ▁▁            ▁                ▁  
  █▁▁▁▁▁▁▁█▁▁█▁▁▁▁▁▁▁▁██▁▁▁▁▁██▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  13.4 ms         Histogram: frequency by time        28.9 ms <

 Memory estimate: 21.88 MiB, allocs estimate: 199166.

--- Matrix-free Jacobian-vector products ---

SparseDiffTools JacVec construction...
BenchmarkTools.Trial: 30 samples with 9 evaluations per sample.
 Range (min … max):  2.425 μs …   2.783 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.465 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.502 μs ± 104.024 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▂▂   █ ▅                                                    
  ████▅▅█▅█▅▁▁▅▁▁▁▅▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▁▁▁▅▁▅▁▁▅ ▁
  2.43 μs         Histogram: frequency by time        2.78 μs <

 Memory estimate: 3.00 KiB, allocs estimate: 21.

SparseDiffTools JacVec apply (J*v)...
BenchmarkTools.Trial: 100 samples with 1 evaluation per sample.
 Range (min … max):  69.668 μs … 80.895 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     73.285 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   73.533 μs ±  2.066 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                     █    ▃                                    
  ▄▄▄▁▁▁▇▅▄▇▄▇▇▄▅█▅▁▇█▅▅▅██▇▇▅▅█▇▁▄▇▁▁▅▁▄▄▅█▅▁▇▄▁▄▇▁▁▁▁▁▄▁▁▁▄ ▄
  69.7 μs         Histogram: frequency by time        78.6 μs <

 Memory estimate: 146.38 KiB, allocs estimate: 1318.

DifferentiationInterface pushforward J*v (pre-computed prep)...
BenchmarkTools.Trial: 100 samples with 1 evaluation per sample.
 Range (min … max):  69.057 μs … 118.916 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     74.125 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   75.373 μs ±   6.828 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

     ▁▁▄█▁▃                                                     
  ▄▃▅██████▇▆▅▄▁▃▃▁▃▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃ ▃
  69.1 μs         Histogram: frequency by time          114 μs <

 Memory estimate: 146.38 KiB, allocs estimate: 1318.

DifferentiationInterface pullback J'*v (pre-computed prep)...
BenchmarkTools.Trial: 100 samples with 1 evaluation per sample.
 Range (min … max):  10.519 ms … 19.278 ms  ┊ GC (min … max):  0.00% … 26.68%
 Time  (median):     13.612 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   14.030 ms ±  2.504 ms  ┊ GC (mean ± σ):  14.41% ± 15.22%

      ▆ █ ▆                            ▂▄  ▄                   
  ▆▄▆██▆████▄▆▁▄▁▆▄▄▄████▁▄▄▄▁▁▄▁▁▄█▄▆▆███▄█▆▁▄▆▁▁▄▁▁▁▁▄▄▄▄▁▆ ▄
  10.5 ms         Histogram: frequency by time        19.2 ms <

 Memory estimate: 21.63 MiB, allocs estimate: 198004.

--- Full system solve (GMRES Newton-Krylov, mirrors GXBeam xpfunc path) ---

SparseDiffTools GMRES Newton (JacVec matvec)...
BenchmarkTools.Trial: 30 samples with 1 evaluation per sample.
 Range (min … max):  510.516 μs …   5.157 ms  ┊ GC (min … max):  0.00% … 83.97%
 Time  (median):     518.356 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   672.898 μs ± 847.015 μs  ┊ GC (mean ± σ):  21.45% ± 15.33%

  █                                                              
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅ ▁
  511 μs        Histogram: log(frequency) by time       5.16 ms <

 Memory estimate: 653.38 KiB, allocs estimate: 2366.

DifferentiationInterface GMRES Newton (pushforward matvec)...
BenchmarkTools.Trial: 30 samples with 1 evaluation per sample.
 Range (min … max):  488.312 μs …   5.068 ms  ┊ GC (min … max):  0.00% … 83.85%
 Time  (median):     501.444 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   654.742 μs ± 833.606 μs  ┊ GC (mean ± σ):  21.63% ± 15.31%

  █                                                              
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅ ▁
  488 μs        Histogram: log(frequency) by time       5.07 ms <

 Memory estimate: 655.27 KiB, allocs estimate: 2361.

GXBeam static_analysis! xpfunc path (JacVec+GMRES baseline)...
BenchmarkTools.Trial: 30 samples with 1 evaluation per sample.
 Range (min … max):  609.356 μs … 981.125 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     617.601 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   646.794 μs ±  79.063 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █▅                                                             
  ██▆▅▅▃▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃ ▁
  609 μs           Histogram: frequency by time          981 μs <

 Memory estimate: 673.00 KiB, allocs estimate: 2615.

==============================================================================================================
SUMMARY
==============================================================================================================

| Method | Min (ms) | Mean (ms) | Std (ms) | Med (ms) | Max (ms) | Allocs | Mem (MB) |
|---|---|---|---|---|---|---|---|
| SparseDiffTools compute-only | 0.109 | 0.143 | 0.034 | 0.135 | 0.283 | 1322 | 0.35 |
| DiffInterface   compute-only | 0.101 | 0.110 | 0.013 | 0.105 | 0.157 | 1324 | 0.35 |
| SparseDiffTools full pipeline | 2.104 | 3.927 | 2.904 | 2.393 | 9.256 | 18810 | 8.66 |
| DiffInterface   full pipeline | 13.350 | 19.129 | 4.844 | 18.878 | 28.873 | 199166 | 22.94 |
| SparseDiffTools JacVec construct | 0.002 | 0.003 | 0.000 | 0.002 | 0.003 | 21 | 0.00 |
| SparseDiffTools JacVec apply  J*v | 0.070 | 0.074 | 0.002 | 0.073 | 0.081 | 1318 | 0.15 |
| DiffInterface   pushforward   J*v | 0.069 | 0.075 | 0.007 | 0.074 | 0.119 | 1318 | 0.15 |
| DiffInterface   pullback      J'*v | 10.519 | 14.030 | 2.504 | 13.612 | 19.278 | 198004 | 22.68 |
| SparseDiffTools GMRES Newton (JacVec) | 0.511 | 0.673 | 0.847 | 0.518 | 5.157 | 2366 | 0.67 |
| DiffInterface   GMRES Newton (pushforward) | 0.488 | 0.655 | 0.834 | 0.501 | 5.068 | 2361 | 0.67 |
| GXBeam xpfunc path (JacVec+GMRES baseline) | 0.609 | 0.647 | 0.079 | 0.618 | 0.981 | 2615 | 0.69 |

**Compute-only speedup (ST/DI):**      1.07x — DI faster
**Full pipeline speedup (ST/DI):**      0.16x — ST faster
**JVP speedup (ST apply / DI pf):**     1.01x — DI faster
**GMRES Newton speedup (ST/DI):**       1.05x — DI faster
**GMRES vs GXBeam xpfunc (ST/GXBeam):** 0.84x — ST faster
**GMRES vs GXBeam xpfunc (DI/GXBeam):** 0.8x — DI faster