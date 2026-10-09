# Network gradient regression comparison

The second optimization uses a centered linear regression compiled with Numba
for well-conditioned windows. It preserves the 0.001 gradient floor, computes
R2 from the unclamped slope, and keeps NumPy's float32 variance denominator for
Flow elevations loaded from disk. Short, non-finite and poorly conditioned
windows retain the original `np.polyfit` path. Window selection and channel
traversal order are unchanged. Floating-point results are numerically equivalent
within the comparisons below, rather than bit-for-bit identical.

```bash
python benchmarks/validate_network_regression.py --output validation.json
python benchmarks/benchmark_network.py --baseline 6958244 --flow /path/jebja30_flow.tif --gradients --compare-gradients
```

Against `6958244`, all 360 gradient/R2 array comparisons passed at rtol=1e-10
and atol=1e-12. They cover small25, tunez and jebja30; newly created and saved/
reloaded Flow objects; theta 0.25, 0.45 and 0.65; and npoints 1, 2, 5, 10 and 20.
The largest absolute difference was 5.07e-12. Cell IDs, areas, elevations,
distances and chi were compared exactly. Maxima and versions are recorded in
`network_regression_validation.json`.

With `gradients=True` on saved jebja30 Flow, three fresh-process runs had median
construction times of 2.515 s before and 0.826 s after (about 3.0 times faster).
Median whole-process peak RSS was 242.30 versus 241.37 MiB. This compares against
the preceding compact-memory implementation, not the original Network. Imports,
Flow loading and compilation warmup are outside the timer. Raw observations
are in `network_regression_results.json`. These three-run results are indicative
and do not establish statistical significance or performance on a large DEM.

# Network construction comparison

Activate the GDAL/PROJ environment before running:

```bash
python benchmarks/benchmark_network.py --gradients
python benchmarks/benchmark_network.py --flow /path/jebja30_flow.tif --gradients
```

The baseline is `6c351a43f46ac8c73ec98b0b1eab500c6c53318d`. Each pair
runs sequentially in fresh processes, with JIT warmup outside the construction
timer. The script checks SHA-256 equality of all eleven Network arrays,
including their data types. RSS includes imports, warmup and Flow loading.
The default synthetic case is a chain of 4,095 edges in a 4,096 x 4,096 layout;
it exposes memory costs for sparse networks and is not a realistic DEM.

Medians of three runs with `gradients=True`:

| Case | Baseline time | Compact time | Baseline peak RSS | Compact peak RSS |
| --- | --- | --- | --- | --- |
| Synthetic sparse chain | 4.918 s | 2.262 s | 557.13 MiB | 280.75 MiB |
| jebja30, 5,947 network cells | 3.708 s | 2.863 s | 262.89 MiB | 241.84 MiB |

Raw observations and tool versions are in `network_sparse_results.json` and
`network_jebja_results.json`. These small samples are indicative, not a
statistical significance claim or a measurement on the user's large DEM.
The new implementation sorts compact node IDs, so dense networks may have
different time/memory tradeoffs. Flow accumulation still requires its raster.

Separate validation compared 180 arrays and data types exactly against the
baseline across small25, tunez and jebja30, using default, 1, 20 and above-maximum
thresholds. Default-threshold networks included slope/ksn and both R2 arrays;
chi was also recalculated with theta=0.3 and a0=2. The gradient regression and
window rules are preserved; existing Strahler and export defects are outside
this optimization.

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

## Outside-flat barriers and weight precision

Outside-flat friction is now `np.inf` rather than 99999. MCP no longer
traverses non-flat cells to connect flat regions. Returned outside-flat
weights remain -99999, and the no-presill fallback is unchanged. With at
least one presill, unreachable flat cells retain infinite distance.

Stored distances remain float64. A separate float32 experiment rounded
some distinct distances into ties: tunez changed ordering in 8 of the 24
fixture configurations, and three giver cells changed their receiver in
the default configuration. This is a routing change, not merely a tiny
reported-distance difference. Barriers alone, keeping float64, matched
all four output arrays and their dtypes in all 24 configurations. Raw
comparisons are in `weights_precision_comparison.json`; receiver counts
align arrays by giver cell before comparison so ordering differences do
not inflate them. All 111 tests passed, including barrier, diagonal-cost
and no-presill tests.

Production float64 construction was also profiled on the same 799.83 MiB
synthetic DEM. Comparison with the preceding recorded 99999-cost run:

| Measurement | Cost 99999 (`5c9a68a`) | Barrier `np.inf`, float64 |
| --- | --- | --- |
| OS-reported peak process RSS | 16326.6 MiB (15.94 GiB) | 16326.9 MiB (15.94 GiB) |
| Weights time | 270.696 s | 78.952 s |
| Complete Flow time | 479.689 s | 293.112 s |

In this run, weight calculation took about 71% less time and complete Flow
construction took about 39% less time. Global peak memory did not improve.
All three output checksums matched exactly. One run per version is not a
statistical speed benchmark, and results depend on flat extent and geometry.
Raw production measurements are in `flow_memory_800_inf_phases.json`.

```bash
python benchmarks/profile_flow_memory.py --size 14480 --output inf-phases.json
```

On the 799.83 MiB DEM, the experimental barriers-plus-float32 run took
291.839 s overall (80.269 s in weights), with peak RSS 16327.2 MiB
(15.94 GiB). Its `_ix` and `_ixc` checksums changed. Float32 reduced the
sampled sorting peak to 7.28 GiB, but did not reduce the overall peak:
MCP's internal calculation still uses float64. These experimental numbers
do not describe the production float64 constructor. Raw experimental
measurements are in `flow_memory_800_inf_float32_phases.json`.

The profiler exposes float32 only as an experimental worker option:

```bash
python benchmarks/profile_flow_memory.py --size 14480 --weights-dtype float32 --output experimental-phases.json
```

The profiler now defaults to `--weights-dtype native`, using production code
unchanged. Explicit dtype options only cast the stored result; they do not
change the solver's accumulation precision or restore lost precision.

## Compiled Dijkstra replacement

Flow now uses a Numba-compiled eight-neighbour Dijkstra solver instead of
`MCP_Geometric`. Edge costs retain MCP's operation order:
`length * 0.5 * (old_cost + new_cost)`, with cardinal length 1 and diagonal
length sqrt(2). Seeds start at zero; infinite outside-flat friction blocks
propagation. Returned weights still add one and assign -99999 outside flats.
The no-presill fallback and float64 distance precision are unchanged.

The indexed binary heap has one entry per frontier cell, supports decrease-key
and grows only when needed. Heap entries and per-cell heap positions use
int32 when the raster fits, otherwise int64. The solver omits traceback,
full-raster edge maps, separate heap priority arrays and Fortran-layout copies.
Costs and output can use C order directly. Heap growth temporarily holds both
old and new buffers; this remains an in-memory solver with full-raster float64
distances and a full-raster position array.

Validation against the previous MCP revision (`09e481d`):

- 100 random friction surfaces, including barriers, multiple/duplicate seeds
  and uniform surfaces: every distance matched exactly (maximum error zero).
- small25, tunez and jebja30 across every `filled`, `raw_z` and `auxtopo`
  combination: `_ix`, `_ixc`, `_zx` and `_nodata_pos` matched exactly in all
  24 configurations. No giver changed its receiver.
- All 115 library tests passed; existing tests were retained. Four new
  solver tests cover the MCP reference, narrow rasters, unreachable cells,
  duplicate/invalid seeds, heap growth and the int64 index specialization.

The random comparison found exact equality; reference tests additionally
allow 1e-14 tolerance for floating-point differences on other environments.
These checks establish agreement on the tested inputs, not a guarantee for
every terrain. Raw comparison data are in `dijkstra_validation.json`.

The benchmark scripts now warm both compiled kernels before timing. A 256
square smoke DEM matched the MCP reference and the phase profiler, including
array dtypes. The new solver has **not yet been measured on the 800 MiB DEM**;
the large measurements above describe the earlier MCP implementation.

## Real DEM comparison

A user-provided `DEM_30m.tif` was reconstructed from four 7z volumes and
benchmarked without resampling, dtype conversion or cropping. It has
16,627 columns by 14,448 rows (240,226,896 cells), int16 storage and a
458.20 MiB elevation array. NoData is -9999 and occupies 22.32% of cells;
valid elevations range from -32 to 3733. The TIFF and preview remain local.

Default Flow construction (`filled=False`, `raw_z=False`, `auxtopo=False`)
ran sequentially in separate processes: MCP from revision `09e481d`, then
compiled Dijkstra from `4b29508`. Both use the current Priority-Flood fill
and float64 weights, with infinite friction outside flats. NoData behavior
is the existing library behavior in both versions. Input loading and kernel
warmup precede timing; output hashing and saving follow timing.

| Measurement | MCP | Compiled Dijkstra |
| --- | --- | --- |
| OS peak process RSS | 19.51 GiB | 13.71 GiB |
| Complete Flow construction | 218.832 s | 146.491 s |
| Sampled weights peak RSS | 19.45 GiB | 9.07 GiB |
| Weights time | 100.098 s | 28.976 s |
| Sampled sorting peak RSS | 11.05 GiB | 9.24 GiB |
| Sampled receivers peak RSS | 13.61 GiB | 13.62 GiB |

Global peak RSS fell by 5.80 GiB (29.7%), and complete construction time by
33.1%. Weight calculation took 71.1% less time. The largest remaining peak
is receiver construction, so further solver memory reductions alone would
not remove the current global peak. Phase peaks are total process RSS and
are not additive; 50 ms sampling can miss brief peaks. These are one-run
measurements per solver on this cloud machine.

Output arrays were saved outside the timed region, then compared directly
in chunks as well as by checksum. `_ix`, `_ixc` and `_zx` each have
186,573,695 elements and zero differences; `_nodata_pos` has 53,630,198
elements and zero differences. Shapes and dtypes matched exactly. The input
SHA-256, phase measurements and comparison counts are recorded in
`real_dem_mcp_dijkstra_results.json`; no user raster data are committed.

The profiler now accepts a local input file, a reference Git revision and
optional output-array saving:

```bash
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --revision 09e481d --save-arrays /tmp/mcp-arrays --output mcp-phases.json
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --save-arrays /tmp/dijkstra-arrays --output dijkstra-phases.json
```

The current checkout must contain the compiled solver; the reference
revision selects the historical Flow implementation. A small fixture smoke
test also confirmed matching arrays/dtypes in both file-based profiler modes.

## Compiled receiver selection

Receiver calculation now builds the inverse sort ranks with Numba, examines
cardinal and diagonal neighbours directly, computes gradients in scalar
buffers and writes a single uint32 output. It no longer creates full-raster
dilations, copies of those dilations, two candidate arrays, gradient arrays,
selection masks or a final cast copy. For the default contiguous layout, its
main new buffers are the rank map and output: about 1.79 GiB on this DEM.

This preserves the original rank-first candidate rules, centre participation,
SciPy's reflected boundary behaviour, integer subtraction overflow and NumPy's
division promotion rules. It does not replace those rules with a different
steepest-neighbour algorithm. Unsupported floating-point precisions (float16
and extended floats) use bounded NumPy arithmetic blocks with compiled
neighbour selection rather than converting the entire DEM.

Validation includes 288 receiver reference cases across 12 numeric dtypes,
four raster shapes (including single rows/columns), C/F indexing and three
cell-size scalar variants. All 24 Flow fixture configurations matched exactly,
including array dtypes. The full suite passes 117 tests; existing tests were
retained and two receiver regression tests were added.

The same complete `DEM_30m.tif` was profiled again, with all defaults and
kernel warmup outside timing. The baseline is the preceding recorded Dijkstra
run with dilation-based receivers, not a newly repeated baseline run.

| Measurement | Previous receivers | Compiled receivers |
| --- | --- | --- |
| OS peak process RSS | 13.71 GiB | **9.24 GiB** |
| Complete Flow time | 146.491 s | **134.710 s** |
| Sampled receivers peak RSS | 13.62 GiB | **4.31 GiB** |
| Receivers time | 22.099 s | **16.012 s** |

Global peak RSS fell by 4.47 GiB (32.6%). The receivers phase peak fell by
68.4%; its time fell by 27.5%. Total time fell by 8.0% in these single runs,
which also include normal timing variation in other phases. The largest
observed memory peak is now sorting, with weights close behind it.

| Current phase | Sampled peak total process RSS | Time |
| --- | --- | --- |
| Fill | 2.13 GiB | 56.318 s |
| Flats and sills | 7.45 GiB | 11.486 s |
| Presills | 4.96 GiB | 6.175 s |
| Auxiliary surface preparation | 6.01 GiB | 0.414 s |
| Weights | 9.08 GiB | 27.158 s |
| Sorting | **9.24 GiB** | 12.277 s |
| Receivers | 4.31 GiB | 16.012 s |
| Connection filtering | 4.23 GiB | 1.414 s |
| Elevation extraction | 4.56 GiB | 1.233 s |
| Unconnected-cell positions | 4.81 GiB | 1.758 s |

Every element of `_ix`, `_ixc`, `_zx` and `_nodata_pos` was compared against
the previous saved arrays, with zero differences. Checksums and dtypes also
matched. Raw measurements and validation are in `real_dem_receivers_results.json`.
User raster data and full arrays remain outside Git. Phase RSS is total worker
memory, not additive; 50 ms sampling can miss brief spikes.

The profiler and full-construction comparison now also warm the receiver
kernel before timing. To reproduce the before/after comparison on this checkout:

```bash
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --revision a9d0dba --output receivers-before.json
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --output receivers-after.json
```

## Unit-friction mask and native float32 accumulation

At the user's request, the default `auxtopo=False` path now omits the friction
raster and uses the flat mask directly. Cardinal steps cost float32(1),
diagonal steps float32(sqrt(2)); distances are accumulated and rounded to
float32 before heap comparisons. This is not the earlier experiment that
cast a float64 solver result only at the end. `auxtopo=True, filled=False`
retains the previous float64 variable-cost calculation. As before, `filled=True`
disables auxiliary topography, so it uses the new unit-friction path.

Outside-flat cells block entry; returned weights there remain -99999. The
final +1 is float32. Without presills, flat weights remain 1. `_zx` remains
float64. A weights raster is still needed for sorting; the removed raster is
the separate **friction/cost surface**.

All 120 tests pass, retaining the previous tests and adding mask/float32 and
fallback checks. Across 24 fixture configurations, all giver sets, heights
aligned by giver and NoData positions are unchanged. Eighteen configurations
match all arrays exactly. In the default Tunez configuration, four giver cells
choose a different receiver (0.00131%); small25 and jebja30 remain identical.
The auxiliary-topography float64 route remains identical on all three fixtures.

For review, `tunez_unit_float32_changed_receivers.csv` lists the four cells,
coordinates and old/new receiver indices. The accompanying GeoTIFF mask has
1 at those cells and 0 elsewhere, with Tunez's original georeferencing. These
artifacts derive from the existing repository fixture, not the user's DEM.
The full configuration comparison is in `unit_float32_fixture_comparison.json`.

On the full user DEM, compared with the preceding recorded `51b9fe9` run:

| Measurement | Cost raster + float64 | Flat mask + float32 |
| --- | --- | --- |
| OS peak process RSS | 9.24 GiB | **8.34 GiB** |
| Sampled weights peak RSS | 9.08 GiB | **6.17 GiB** |
| Weights time | 27.158 s | **22.882 s** |
| Complete Flow time | 134.710 s | **128.643 s** |

The overall peak fell by about 9.7%; total time fell about 4.5% in these single
runs. Sorting still sets the overall peak. Removing a 1.79 GiB friction raster
and halving distance storage substantially reduces the weights phase, but
sorting still retains the float32 weights and its index arrays.

This mode intentionally changes routing. Of 186,573,695 common giver cells
in `DEM_30m.tif`, 11,466 change receiver (0.0061456%). No giver was added or
removed; heights aligned by giver and NoData positions remain identical.
The `_ix` and `_ixc` hashes change; `_zx` matches. Small direction changes can
also affect downstream accumulation and basin assignments; their consequences
must not be inferred solely from the fraction of changed receivers.

Raw profiles and aligned comparison counts are in
`real_dem_unit_float32_results.json`. Detailed changed directions for the
user's DEM remain local, outside Git. The profiler records the observed native
weights dtype, warms both solver precisions, and uses `native` by default.

```bash
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --revision 51b9fe9 --output unit-before.json
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --output unit-after.json
```

`benchmark_flow_memory.py --allow-differences` reports changed hashes instead
of requiring equality for such deliberate precision experiments; dtype checks
on the final Flow arrays remain enabled.

## Single stable lexicographic ordering

`sort_dem()` now uses `np.lexsort((-weights, -elevations))` instead of two
stable argsorts followed by index gathering. The last key is primary, so
elevations descend first, then weights within each elevation. Lexsort's
stability retains original raster index order for exact ties; no full-raster
`arange` key is needed. Key dtypes and their negation behaviour are preserved.
Output indices remain uint32, with the same C/F indexing semantics.

All 121 tests pass. A new regression test compares the previous implementation
on integer extremes, NaN, infinities and signed zero in both indexing orders.
All four Flow arrays and their dtypes match exactly in all 24 fixture
configurations. A synthetic 256-square smoke comparison also matched.

The full user DEM was profiled again and compared with the preceding recorded
`998176b` run (not a newly repeated baseline):

| Measurement | Two stable argsorts | Lexsort |
| --- | --- | --- |
| OS peak process RSS | 8.34 GiB | **7.67 GiB** |
| Sampled sorting peak RSS | 8.34 GiB | **7.45 GiB** |
| Sorting time | 12.235 s | **9.814 s** |
| Complete Flow time | 128.643 s | **127.692 s** |

Lexsort still requires internal sort workspace and full-raster negated keys;
the memory gain must not be inferred just by counting removed Python arrays.
The sorting phase peak fell by about 0.89 GiB. The largest sampled phase peak
is now flats/sills detection (7.60 GiB), so the overall memory reduction is
smaller than the sorting reduction. Total times are effectively similar in
these single runs; the phase times include normal variation elsewhere.

Direct chunked comparison of every element in `_ix`, `_ixc`, `_zx` and
`_nodata_pos` found zero differences, including shapes and dtypes. This change
introduces no additional direction changes beyond the existing authorized
float32 routing variant. Checksums also match. Raw measurements and fixture
validation are in `real_dem_lexsort_results.json`; raster data and full arrays
remain local. Phase RSS is sampled every 50 ms, so brief peaks can be missed;
the OS high-water mark is the global value reported above.

```bash
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --revision 998176b --output sorting-before.json
python benchmarks/profile_flow_memory.py --dem /path/to/DEM_30m.tif --output sorting-after.json
```
