# Flow.drainageBasins optimization

Last modified: 2026-10-16 11:40:19 CEST (Europe/Madrid).

The benchmark uses the previously prepared Flow for `DEM_30m.tif`: 240,226,896
raster cells and 186,573,485 topologically ordered edges. The DEM SHA-256 is
`53c8d23885e8fff575a6eed5511a535f02b3221be80da7aa617bf84aa04241b9`.
The baseline implementation matches commit `9f1f422`, which already includes
the optimized Flow accumulation and Network construction.

| Mode | Baseline time | Optimized median | Speedup | Baseline peak RSS | Optimized peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| All basins (`min_area=0`) | 91.954 s | 1.040 s | 88.4x | 5.139 GiB | 5.175 GiB |
| Default (`min_area=0.005`) | 111.995 s | 3.662 s | 30.6x | 10.571 GiB | 5.088 GiB |
| Ten specified outlet coordinates | 49.689 s | 0.268 s | 185.6x | 3.543 GiB | 3.584 GiB |

Optimized times are medians of three separate warm processes per mode. The
ranges were 0.870–1.150 s (all), 3.005–3.933 s (default), and 0.258–0.274 s
(specified outlets). Each baseline was measured once. A separate default run
with an empty Numba cache and no warmup took 3.847 s and also matched exactly.
These are measurements from this DEM and machine, not performance guarantees.

The default-mode peak RSS falls by 51.9%. Other modes retain the required dense
int64 output and have similar RSS, with a small increase from the compiled
runtime. RSS includes the loaded Flow and runtime; it is sampled externally
every 10 ms during drainageBasins. Cached-Flow loading, result comparison,
hashing and file writing are excluded. DEM/Flow creation is also excluded.

Machine: Linux, 4-vCPU quota, 32 GiB memory limit, Python 3.13.5, NumPy 2.5.4,
Numba 0.68.0 and GDAL 3.10.3. The Flow cache stores arrays directly as `.npy`
without the elevation quantization of Flow GeoTIFF serialization.

## Implementation and compatibility

- A compiled reverse topological traversal labels basins directly in the
  output array, replacing the Python loop.
- Automatic outlet selection uses a boolean raster and one forward topological
  pass. A receiver becomes a candidate outlet; its later selected outgoing edge
  clears it. This avoids the sparse adjacency matrix, degree arrays and edge
  copies previously constructed through streamPoi.
- The accumulation scratch array is released before allocating basin labels.
- Output remains int64, with zero for unassigned cells and Grid NoData=0.
  Basin numbering, strict accumulation > threshold selection, outlet ordering,
  and the handling of nested and duplicate seeds are unchanged.
- Existing edge-case semantics are retained: one outlet gets ID 1 even if a
  custom ID is supplied; an empty seed set (including no automatically selected
  outlets) falls back to labeling all basins. This optimization does not change
  those policies.

All 240,226,896 labels match the baseline exactly in all nine final runs and
the cold default run. The all-basin result has 23,045 basin IDs; the default
result has 13. Manual coordinates are taken from ten evenly spaced positions
in the cached receiver array and supplied IDs 11–20.

The regression tests cover signed/unsigned edge indices, exact basin numbering,
minimum-area thresholds, nested/duplicate seeds, custom and negative IDs,
single or empty outlet sets, coordinate forms and errors, Grid metadata and
empty Flow objects. The full suite, including these checks, passes:
134 tests and 465 subtests, with 31 warnings.

## Reproduction

Use the development Python environment plus `psutil` for RSS sampling. The
runner is [`drainage_basins.py`](drainage_basins.py). Prepare a Flow cache once
using `network_construction.py --prepare`, as described in the Network benchmark
report. Each dense labels file for this DEM occupies about 1.79 GiB; they and
the DEM/cache are kept outside the repository.

From the repository root, use an external cache and output directories:

```bash
python benchmarks/drainage_basins.py --cache /tmp/landspy-flow \
  --mode default --reference-ref 9f1f422 --warmup --save-labels \
  --output /tmp/basins-before
python benchmarks/drainage_basins.py --cache /tmp/landspy-flow \
  --mode default --warmup --compare /tmp/basins-before/labels.npy \
  --output /tmp/basins-after
```

Repeat with distinct output directories and `--mode all` or `--mode seeded`.
For a first-compilation measurement, omit `--warmup` and set `NUMBA_CACHE_DIR`
to a new empty directory. `--reference-ref` reads an archived source tree from
the local Git revision without changing the working checkout. The runner was
also validated against that archived revision in all three modes on small25.

Raw timing/RSS results and label digests are saved in
[`drainage_basins_results.json`](drainage_basins_results.json). Its modification
metadata is recorded in its `last_modified` field; the interpretation of those
measurements is documented here because JSON does not support comments.
