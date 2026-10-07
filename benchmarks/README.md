# DEM fill comparison

Run after installing LandSpy and activating the GDAL/PROJ environment:

```bash
python benchmarks/benchmark_fill.py --size 4096
python benchmarks/benchmark_fill.py --size 14480
```

The current script measures Priority-Flood only. It reports time, peak RSS
and an output checksum on deterministic synthetic float32 terrain. Input size
is the uncompressed array size, not the size of a GeoTIFF file.

The following historical measurements were collected before removal of the
reconstruction option and NoData restoration. They are retained as a record
of that earlier comparison, not validation of the current NoData behavior.

Measured with Python 3.12.14, NumPy 2.2.6, scikit-image 0.26.0 and Numba 0.68.0
on the cloud Linux environment:

| Input | NoData | Reconstruction time | Priority-Flood time | Reconstruction peak RSS | Priority-Flood peak RSS |
| --- | --- | --- | --- | --- | --- |
| 64 MiB, 4096 x 4096 | 0% | 9.420 s | 6.377 s | 1531.5 MiB | 598.4 MiB |
| 64 MiB, 4096 x 4096 | 50% | 6.550 s | 2.798 s | 1531.5 MiB | 541.9 MiB |
| 799.83 MiB, 14480 x 14480 | 0% | 132.900 s | 103.986 s | 16620.3 MiB | 4112.5 MiB |

All three pairs produced matching checksums. The large pair's SHA-256 was
`44cc7eb0cc54ddb45b084d56876f9f0a5a8bd9cf0a66517e393f1ee6274e4e96`.
The full library test suite also passed: 107 tests, including comparison with
reconstruction on the three repository DEMs, random terrain, NoData and
eight-neighbour drainage.

Fill times exclude warmup, input generation and output hashing. Warmup is
reported separately; cached float32 Priority-Flood warmup took 0.118-0.123 s.
The first compilation on a new environment can take longer. Peak RSS covers
the whole process, including imports, warmup, input and algorithm arrays.
Results depend on the terrain and machine; these synthetic measurements do
not substitute for a comparison on a user's real DEM.

Both algorithms were called with `inplace=True`. That replaces the DEM's
array after filling. Priority-Flood still needs a working output copy,
visitation array, FIFO and a terrain-dependent heap; it is not an out-of-core
algorithm. Reconstruction here includes the earlier removal of redundant
output copies, so the comparison measures the algorithm change itself.

## Flow construction: temporary-memory changes

Compared with commit `3390f70`, Flow construction now skips unused elevation
differences, frees weights before receiver calculation, uses a boolean cell
marker, avoids redundant sorting indices and filters connections once.
The final `_ix`, `_ixc` and `_zx` representations are unchanged.

Exact array and dtype equality was checked on `small25`, `tunez` and `jebja30`
for all combinations of `filled`, `raw_z` and `auxtopo` (24 configurations).
The derived NoData positions were also identical. The full suite passed:
108 tests, including stable sorting and unconnected-cell regressions.

On a deterministic 2048 x 2048 synthetic float32 DEM (16 MiB), separate
processes reported peak RSS of 657.6 MiB before and 647.3 MiB after. All three
output checksums matched. This small reduction in overall peak memory shows
that other construction phases still dominate this terrain; savings from
individual temporary arrays must not be added together or extrapolated
directly to large DEMs.

The complete comparison was subsequently run on a 14480 x 14480 synthetic
float32 DEM (799.83 MiB), using seed 38, with default Flow settings:

| Revision | Peak process RSS | Flow construction time |
| --- | --- | --- |
| Before the five changes (`3390f70`) | 18134.8 MiB (17.71 GiB) | 470.055 s |
| After the five changes (`7ec8c85`) | 18123.9 MiB (17.70 GiB) | 467.716 s |

The three output checksums matched exactly. The peak reduction was only
10.9 MiB (0.06%). The time difference was about 0.5%; one run per revision
does not establish a repeatable speed improvement. These changes reduce
temporary allocations in individual phases but did not materially lower
the overall construction peak on this large terrain. The phase responsible
for the overall peak needs further profiling before choosing the next change.

Raw results are in `flow_memory_800_results.json`. To reproduce the comparison
with an optimized checkout and local history containing the baseline:

```bash
python benchmarks/benchmark_flow_memory.py --size 14480 --baseline 3390f70
```

The script checks checksums and dtypes in sequential, separate processes.
Both constructors use the same current Priority-Flood implementation.
Timing excludes input generation, fill-kernel warmup and output hashing.
Peak RSS includes the entire worker process and the retained input DEM.
Allow more than 18 GiB of free RAM; this full comparison took about 16 minutes
on the cloud machine. These measurements describe uncompressed synthetic
data, not a user's real GeoTIFF.

## Flow memory by construction phase

The optimized constructor was profiled on the same 799.83 MiB DEM, with
default settings. Linux RSS was sampled every 50 ms by a separate process,
so the sampler continues while C extensions hold the worker's Python GIL.
Phase markers use existing verbose progress messages; only the benchmark
worker wraps the final unconnected-cell calculation to mark its start.

| Phase | Sampled peak process RSS | Time |
| --- | --- | --- |
| Priority-Flood fill | 4.01 GiB | 105.563 s |
| Flats and sills | 6.65 GiB | 12.416 s |
| Presills | 2.86 GiB | 2.591 s |
| Auxiliary surface preparation (`auxtopo=False`) | 4.42 GiB | 1.093 s |
| Cost-distance weights (`MCP_Geometric`) | **17.60 GiB** | **267.054 s** |
| Sorting | 9.66 GiB | 44.014 s |
| Receivers | 11.56 GiB | 32.460 s |
| Connection filtering | 4.32 GiB | 2.725 s |
| Elevation extraction | 5.43 GiB | 2.740 s |
| Unconnected-cell positions | 4.59 GiB | 2.631 s |

These are total worker RSS values within each phase, not additional memory
or allocations that can be added together. Sampling can miss short spikes:
the OS-reported process high-water mark was 18125.5 MiB (17.70 GiB), compared
with the largest sampled phase peak of 18018.8 MiB (17.60 GiB). Phase times
include intervening bookkeeping until the next marker; boundary RSS values
have approximately one sampling interval of uncertainty.

The complete constructor took 473.400 s. The weights phase accounted for
about 56% of that time and set the largest observed memory peak. Receiver
optimization alone would not remove that peak: reducing the cost-distance
phase is the next priority, preserving flat routing and output ordering.

All three output checksums matched the previous uninstrumented 800 MiB run.
The instrumentation also matched ordinary Flow construction on a smoke DEM.
Raw phase measurements are in `flow_memory_800_phases.json`.

```bash
python benchmarks/profile_flow_memory.py --size 14480 --output phases.json
```

This sampler requires Linux `/proc`. It does not change production Flow code.

## Avoiding weight-calculation copies

Three changes preserve the existing MCP routing: release `sills` immediately
after presill extraction, prepare the friction surface as Fortran-contiguous
float64 (MCP's native layout), and add one to cumulative costs in place after
releasing the MCP object. For `auxtopo=True`, auxiliary surface computation
still uses float32 before conversion, preserving the previous numerical costs.
The outside-flat cost remains 99999; no cells are newly blocked.

The same default Flow constructor and 799.83 MiB synthetic float32 DEM were
profiled again in a separate process with 50 ms RSS sampling:

| Measurement | Before (`b83835d`) | After the three changes |
| --- | --- | --- |
| OS-reported peak process RSS | 18125.5 MiB (17.70 GiB) | 16326.6 MiB (15.94 GiB) |
| Sampled peak during weights | 18018.8 MiB (17.60 GiB) | 16295.4 MiB (15.91 GiB) |
| Weights time | 267.054 s | 270.696 s |
| Complete Flow time | 473.400 s | 479.689 s |

The global peak fell by 1798.9 MiB (1.76 GiB), approximately 9.9%. The peak
remains in MCP's weight calculation. Timing did not improve in this run;
one run per version does not establish a repeatable timing difference.
The baseline is the preceding recorded phase measurement, rather than a
newly repeated baseline run. RSS includes the retained original DEM and
the rest of the worker process, not just weights or an individual array.

All three 800 MiB output checksums (`_ix`, `_ixc`, `_zx`) matched exactly.
Arrays and dtypes, including `_nodata_pos`, also matched the baseline on
small25, tunez and jebja30 for every combination of `filled`, `raw_z` and
`auxtopo` (24 configurations). All 108 tests passed with the activated
Conda environment and discovery restricted to `test_*.py`.

Raw after-change measurements are in `flow_memory_800_weights_phases.json`.
To reproduce on this checkout:

```bash
python benchmarks/profile_flow_memory.py --size 14480 --output weights-phases.json
```
