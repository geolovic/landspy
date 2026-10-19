# BNetwork construction and channel extraction

The benchmark uses the supplied 240,226,896-cell DEM, its cached Network with
threshold 1500 (3,405,387 nodes), and basin 4 of the four user outlets snapped
to channels. This basin has 83,286 network nodes, a 7,437,000-cell bounding
rectangle and 779 channels. The baseline is commit `7dd9a20`.

Three separate processes per version were measured after warming the new Numba
traversal. Times below are medians; RSS values are the largest observed peaks.

| Operation | Before | After | Speedup | Peak process RSS before / after |
| --- | ---: | ---: | ---: | ---: |
| `BNetwork(net, grid, bid=4)` | 3.793 s | 0.209 s | 18.1× | 8.187 / 1.043 GiB |
| `getChannels("ALL")` | 1.743 s | 0.045 s | 38.9× | 1.175 / 0.929 GiB |
| `getChannel(0)` | 1.845 s | 0.026 s | 69.7× | 1.192 / 0.932 GiB |

RSS includes the parent Network, basin label raster, GDAL, Python and Numba.
It is sampled externally every 10 ms and at phase boundaries. Loading inputs,
Flow/Network creation, Numba compilation/cache initialization and correctness
checks are outside the timed phases. Results describe this basin and machine;
other basin sizes and channel counts will have different gains.
Raw observations are in [bnetwork_results.json](bnetwork_results.json).

## Changes

- The constructor scans basin labels in bounded blocks, samples membership
  only at Network cells, and copies selected Network attributes once. A `Basin`
  stays in its cropped grid instead of being expanded to the original DEM.
- Remapping to the cropped grid uses integer row/column offsets.
- Channel extraction uses compact giver/receiver and owner arrays containing
  only network cells, including terminal receivers. It traces paths with Numba
  and gathers attributes directly, removing the per-channel raster mask and
  complete Network scan.
- `getChannel(id)` traverses preceding channels to preserve connection ownership,
  but constructs only the requested Channel. Length filtering also precedes
  materialization while retaining the original traversal state.

Public `Network.streamPoi`, snapping, chi sensitivity, plotting and I/O are
unchanged. The existing final-point exclusion, traversal termination and
confluence ownership conventions are preserved.

## Correctness

The six benchmark runs have the same SHA-256 over every BNetwork and Channel
attribute (normalizing scalar Python/NumPy representation). All 779 output
channels and all BNetwork arrays were also compared directly with the baseline.
Each single-channel result is checked against the first complete-list result.

New regression tests cover Grid and Basin inputs, three fixture basins,
optional heads, channel counts, length filtering, NaN length, negative channel
indices, connection IDs, unchanged inputs, independent attribute storage,
block boundaries and empty basins. A frozen original traversal supplies the
channel reference. A constructor spy verifies that single-channel extraction
constructs exactly one Channel.

## Reproduce

Run from the repository root with the normal dependencies and `psutil`:

```bash
python benchmarks/bnetwork.py \
  --network /path/to/network.npz --layout /path/to/flow/layout.json \
  --basins /path/to/cuencas.tif --bid 4 \
  --reference-ref 7dd9a20 --output /tmp/bnetwork-before
python benchmarks/bnetwork.py \
  --network /path/to/network.npz --layout /path/to/flow/layout.json \
  --basins /path/to/cuencas.tif --bid 4 --output /tmp/bnetwork-after
```

The Network archive and layout are produced by `network_construction.py`.
Use different output directories for each repetition and compare the `digest`
fields in `metrics.json`. Baseline sources are extracted to a temporary
directory without changing the checkout. The benchmark uses `gradients=False`
Network data; existing gradients, if supplied, are copied through unchanged.
