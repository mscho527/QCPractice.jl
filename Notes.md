# Benchmark Timings on Test Inputs
```bash
minsik in 🌐 telemachus in ~ via △ v3.11.4 took 1m55s
❯ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.5 (2023-01-08)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.8) pkg> activate QCPractice
  Activating project at `~/QCPractice`

julia> using QCPractice, BenchmarkTools

julia> stdout_placeholder = stdout
Base.TTY(RawFD(13) open, 0 bytes waiting)

julia> redirect_stdout(devnull)
Base.DevNull()

julia> mol = Molecule("", "/home/minsik/QCPractice/src/Tests/H2.dat")
Molecule("", "/home/minsik/QCPractice/src/Tests/H2.dat", 2)

julia> r = RHF(mol)
RHF(Molecule("", "/home/minsik/QCPractice/src/Tests/H2.dat", 2), #undef, #undef, #undef, #undef, #undef, #undef, #undef, 5.0e-324, #undef, 0.0)

julia> ben = @benchmark run_rhf!(r, true);

julia> mol = Molecule("", "/home/minsik/QCPractice/src/Tests/H2O.dat")
Molecule("", "/home/minsik/QCPractice/src/Tests/H2O.dat", 10)

julia> r = RHF(mol)
RHF(Molecule("", "/home/minsik/QCPractice/src/Tests/H2O.dat", 10), #undef, #undef, #undef, #undef, #undef, #undef, #undef, 1.12629331307196e-310, #undef, 1.1262933130783e-310)

julia> ben2 = @benchmark run_rhf!(r, true);

julia> mol = Molecule("", "/home/minsik/QCPractice/src/Tests/LiF.dat")
Molecule("", "/home/minsik/QCPractice/src/Tests/LiF.dat", 12)

julia> r = RHF(mol)
RHF(Molecule("", "/home/minsik/QCPractice/src/Tests/LiF.dat", 12), #undef, #undef, #undef, #undef, #undef, #undef, #undef, 1.12629327823914e-310, #undef, 1.1262932782423e-310)

julia> ben3 = @benchmark run_rhf!(r, true);

julia> redirect_stdout(stdout_placeholder)
Base.TTY(RawFD(13) open, 0 bytes waiting)

julia> ben
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  858.951 ms … 893.538 ms  ┊ GC (min … max): 3.57% … 3.87%
 Time  (median):     880.840 ms               ┊ GC (median):    3.80%
 Time  (mean ± σ):   879.845 ms ±  11.366 ms  ┊ GC (mean ± σ):  3.79% ± 0.59%

                                        █
  ▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  859 ms           Histogram: frequency by time          894 ms <

 Memory estimate: 138.89 MiB, allocs estimate: 4978860.

julia> ben2
BenchmarkTools.Trial: 11 samples with 1 evaluation.
 Range (min … max):  477.579 ms … 497.252 ms  ┊ GC (min … max): 2.17% … 3.57%
 Time  (median):     486.405 ms               ┊ GC (median):    3.41%
 Time  (mean ± σ):   486.097 ms ±   6.250 ms  ┊ GC (mean ± σ):  3.42% ± 1.23%

  ▁  ▁   ▁▁                 ▁▁      █     ▁    ▁              ▁
  █▁▁█▁▁▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁██▁▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  478 ms           Histogram: frequency by time          497 ms <

 Memory estimate: 74.71 MiB, allocs estimate: 2534750.

julia> ben3
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  866.334 ms … 892.017 ms  ┊ GC (min … max): 4.87% … 4.88%
 Time  (median):     876.863 ms               ┊ GC (median):    4.06%
 Time  (mean ± σ):   877.530 ms ±   8.357 ms  ┊ GC (mean ± σ):  4.25% ± 0.50%

  █                  █  █    █ █                              █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁█▁▁▁▁█▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  866 ms           Histogram: frequency by time          892 ms <

 Memory estimate: 139.83 MiB, allocs estimate: 4979128.
```
