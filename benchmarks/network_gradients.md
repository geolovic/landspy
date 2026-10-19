# Network gradient optimization

The baseline is commit `5b80da0`. It already includes the earlier Flow, Network
construction and basin optimizations, but still uses the dense Python gradient
traversal and one `numpy.polyfit` call per window.

The large benchmark uses the existing Network from `DEM_30m.tif`, with
`threshold=1500`, `thetaref=0.45`, `npoints=5`, 3,405,387 channel givers and
240,226,896 DEM cells. DEM SHA-256:
`53c8d23885e8fff575a6eed5511a535f02b3221be80da7aa617bf84aa04241b9`.

| Calculation | Baseline | Optimized |
| --- | ---: | ---: |
| Slope alone | 182.213 s | 2.316 s median |
| Ksn alone | 165.978 s | 2.371 s median |
| Slope then ksn, separate public calls | 348.191 s | 4.687 s median |
| Both, one shared traversal as used by the constructor | 348.191 s | 3.033 s median |
| Peak RSS, both calculations | 7.860 GiB | 0.614 GiB maximum, shared traversal |

Three warm optimized runs gave 2.888–3.185 s for the shared traversal and
4.423–4.852 s for separate calls. The shared traversal is 114.8x faster than the
single baseline run, with a 92.2% reduction in peak process RSS. These are
measurements on this data and machine, not portable performance guarantees.
A separate run with an empty Numba cache, no warmup and the same data took
6.624 s and matched the same numerical reference.

The process holds the loaded Network, not the full Flow, in these gradient-only
measurements. Peak RSS includes the Network and Python/native runtimes. Loading,
comparison and saving outputs are excluded from the timed phase. The runner
samples RSS from a separate process every 10 ms, plus phase-boundary readings.

A separate complete `Network(flow, threshold=1500, gradients=True)` run took
8.636 s and peaked at 4.549 GiB RSS with the original Flow also resident. This
excludes prior Flow construction and cache loading. All nongradient Network
fields matched the earlier construction benchmark exactly; the four gradient
fields matched the optimized standalone shared traversal exactly.

The machine has a 4-vCPU quota and 32 GiB memory limit, using Python 3.13.5,
NumPy 2.5.4, Numba 0.68.0 and GDAL 3.10.3. The new traversal is serial because
channel priorities and downstream windows have dependencies; no fast-math
flags are used.

## Implementation

- The correspondence maps channel givers and distinct terminal receivers to
  compact positions. No DEM-sized lookup, gradient or R2 arrays are allocated.
  Head extraction also avoids the former full-raster sparse adjacency matrix.
- The three original eight-byte DEM-sized arrays alone required 5.369 GiB on
  this raster. Compact equivalents scale with the 3.4 million channel cells,
  plus terminal outlets. Sorting and output storage are also channel-sized.
- A fixed reusable window buffer follows the original window growth,
  shrinking and duplicated samples near outlets. Raster-order heads are sorted
  by elevation using the original tie-breaking input order.
- The default regression uses centered float64 sums for slope and residuals.
  This avoids repeated Vandermonde construction and general least-squares solves.
  The gradient floor of 0.001, constant-elevation R2 policy and unclamped
  regression residuals are preserved.
- `Network.polynomial_fit()` remains public and delegates to the same compiled
  fitting kernel used by the traversal. Nonfinite, undersized and numerically
  degenerate cases fall back to the original NumPy regression, including its
  existing exception behavior. No large list of fallback windows is stored.
- Custom/overridden polynomial fits are honored through the same compact
  traversal in Python, with independent window arrays and the original callback
  sequence. Public gradient overrides also retain constructor dispatch.
- Default `Network(..., gradients=True)` shares one topology/window traversal
  for slope and ksn. Independent `calculateGradients(npoints, 'slp')` and
  `calculateGradients(npoints, 'ksn')` calls remain supported.

## Numerical validation

All six final large-data runs and the cold-cache run passed comparisons with
`rtol=1e-10, atol=1e-11`. Every output has the same zero/nonzero mask as the
baseline. The shared and separate optimized results are bitwise identical.
The optimized and baseline regressions are mathematically equivalent, but are
not generally bitwise identical because their numerical algorithms differ.

| Field | Maximum absolute difference | Maximum relative difference, nonzero reference |
| --- | ---: | ---: |
| Slope | 2.445e-12 | 1.019e-11 |
| Ksn | 6.139e-11 | 2.853e-12 |
| Slope R2 | 5.142e-12 | 1.572e-11 |
| Ksn R2 | 4.059e-12 | 2.514e-11 |

Older serialized Flow files can hold float32 elevations. Their previous R2
calculation used float32 variance, whereas the optimized regression uses
float64. On small25, these legacy comparisons differed by at most about 5.1e-8
in R2; gradient differences were below 1e-12. Tests account for that precision
change separately from the native-float64 reference.

The complete suite passes: 144 tests and 548 subtests, with 31 existing warnings.
New coverage includes confluences, equal-height heads, multiple outlets, window
sizes 0 through very large values, empty and short channels, degenerate fits,
constant elevations, read-only inputs, custom fitting callbacks, constructor
fusion and public overrides. A test sets a fictitious DEM size of one trillion
cells while retaining a tiny channel graph to guard against dense lookup
allocations. Numerical reference arrays and all large caches remain outside
the repository.

## Reproduce

Use the development interpreter plus `psutil` for external RSS sampling. The
runner is [`network_gradients.py`](network_gradients.py). It takes the
`network.npz` output and Flow `layout.json` from the Network construction
benchmark. Put input caches and output arrays outside the repository.

```bash
python benchmarks/network_gradients.py \
  --network /tmp/network-cache/network.npz --layout /tmp/flow-cache/layout.json \
  --reference-ref 5b80da0 --mode both --warmup --output /tmp/gradients-before
python benchmarks/network_gradients.py \
  --network /tmp/network-cache/network.npz --layout /tmp/flow-cache/layout.json \
  --mode both --warmup --compare /tmp/gradients-before/gradients.npz \
  --output /tmp/gradients-after
```

Use `--mode separate` for the two public calls, or another `--npoints` value
for different windows. For first compilation, omit `--warmup` and set
`NUMBA_CACHE_DIR` to a new empty directory. `--reference-ref` reads a temporary
archived source tree without modifying the checkout. The runner was checked
against the archived reference on small25 as well as the large cached Network.

Raw timings, RSS measurements, per-field errors and constructor verification
are in [`network_gradients_results.json`](network_gradients_results.json).
