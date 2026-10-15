# Network construction benchmark

The supplied `DEM_30m.tif` has 16,627 columns and 14,448 rows
(240,226,896 cells, Int16, NoData -9999). Its SHA-256 is
`53c8d23885e8fff575a6eed5511a535f02b3221be80da7aa617bf84aa04241b9`.
The archives, DEM, cached Flow and large output arrays stay outside the repository.

The Flow contains 186,573,485 giver/receiver pairs. With `threshold=1500`,
`thetaref=0.45` and `gradients=False`, the Network has 3,405,387 channel cells.
The baseline is commit `82ffe00` (the DEM/numeric correctness fixes).

| Measurement | Baseline | Optimized |
| --- | ---: | ---: |
| Network construction | 129.423 s | 4.845 s median |
| Range of 3 warm optimized runs | — | 4.760–4.896 s |
| Peak process RSS during construction | 8.320 GiB | 4.556 GiB maximum |
| Accumulation called by Network | 85.762 s | 1.047 s median |
| First execution with an empty Numba cache | — | 7.664 s |

This is a 26.7x speedup for Network construction and a 45.2% reduction in its
peak process RSS. The baseline was measured once; these are observations from
this DEM and machine, not portable performance guarantees. Raw results are in
[`network_results.json`](network_results.json).

## Scope and measurement

- Linux container with a 4-vCPU quota and 32 GiB memory limit, Python 3.13.5,
  NumPy 2.5.4, Numba 0.68.0, GDAL 3.10.3.
- Prepare Flow once and cache its arrays directly as `.npy` plus layout JSON.
  This avoids the elevation rounding of `Flow.save()`/`Flow.load()` and uses
  exactly the same graph and elevations in every run.
- Measure construction in separate processes after loading all Flow arrays
  into RAM. Sample RSS externally every 10 ms; the reusable runner also records
  RSS at phase boundaries. RSS includes Flow and the Python/native runtimes;
  it is not just the additional memory owned by Network.
- Warm runs initialize the compiled kernels on a small DEM before timing.
  The cold run uses a separate empty `NUMBA_CACHE_DIR` and no warmup.
- Flow preparation itself took 147.511 s, with 7.618 GiB peak RSS. These costs
  are separate and are not included in the Network speedup. The whole pipeline
  therefore still has the Flow preparation cost and its memory peak.
- The optional `gradients=True` calculation was exercised by the existing
  small-fixture tests, but not benchmarked on the large DEM. Its regression
  fitting and full-raster scratch arrays have not been optimized here.

## Changes and correctness

Flow accumulation now uses a compiled traversal, retaining the existing uint32
unweighted and float64 weighted behavior. NoData masking uses the already-flat
indices directly instead of allocating separate row/column arrays.

Network selection no longer creates a full-DEM threshold mask or gathers a
full copy of all giver accumulation values. Accumulation and selection scratch
arrays are released before calculating distances. Distances and chi use compact
channel receiver positions and compiled reverse traversals, with scratch storage
proportional to channel count instead of total DEM cells. Coordinate arithmetic
and downstream summation order are preserved; no fast-math flags are used.

All 11 persisted Network arrays match the baseline exactly in each of the three
final large-DEM runs, including giver/receiver IDs, areas, elevations, edge
lengths, distances and chi. The slope/ksn arrays are zero for this benchmark.
The full suite passes: 129 tests and 436 subtests, with 31 warnings.
New tests cover confluences, multiple outlets, fractional and empty thresholds,
weighted accumulation, NoData, repeated chi with custom reference area, and
compact lookups for cell IDs near the uint32 limit.

## Reproduce

Use the repository's Python environment with its normal dependencies plus
`psutil` for external RSS sampling. Run from the repository root. The benchmark
script does not alter the checkout when reading the baseline Git revision.
Allow several GiB of additional disk space for the cache and result arrays.

```bash
python benchmarks/network_construction.py --cache /tmp/landspy-flow \
  --prepare /path/to/DEM_30m.tif --output /tmp/network-prepare
python benchmarks/network_construction.py --cache /tmp/landspy-flow \
  --reference-ref 82ffe00 --output /tmp/network-before
python benchmarks/network_construction.py --cache /tmp/landspy-flow \
  --warmup --output /tmp/network-after
NUMBA_CACHE_DIR=/tmp/landspy-new-numba-cache \
  python benchmarks/network_construction.py --cache /tmp/landspy-flow \
  --output /tmp/network-cold
```

Use an empty directory for `--cache` when preparing and a new, empty Numba
cache directory when measuring first compilation. Repeat the optimized command
with distinct output paths for multiple runs. The default threshold is 1500;
use `--threshold 0` for Network's automatic threshold or `--gradients` to measure
optional slope/ksn calculation separately. Construction timing excludes cache
loading and saving `network.npz`.

Compare outputs independently:

```python
import numpy as np

with np.load('/tmp/network-before/network.npz') as before, \
     np.load('/tmp/network-after/network.npz') as after:
    for name in before.files:
        np.testing.assert_array_equal(before[name], after[name])
```
