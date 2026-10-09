# API documentation review — 2026-10-09

Reviewed the 164 function/method/class definitions in `src/landspy` against
signatures and implementation, including private compiled helpers. All now
have docstrings. Inherited methods were checked at their defining class.
Tutorial notebooks, examples, README examples, benchmarks and test code were
not edited. This is an API documentation review, not a claim that every
algorithm or every function in test/benchmark scripts has been validated.

## Documentation corrections

- Corrected copy/reference semantics, strictly positive `Grid.find()` values,
  NoData handling, raster layout assumptions and output dtype conversions.
- Corrected parameter names, accepted inputs, ignored legacy arguments,
  callback types, defaults, threshold inequalities and output array shapes.
- Described Flow persistence, including milliscale uint32 elevation storage
  and float32 loading, rather than implying a lossless GIS raster format.
- Documented chi normalization exactly as implemented, its mutation of state,
  and the fact that recalculating chi does not automatically update ksn.
- Corrected Network/BNetwork/Channel headers and `.swt`/`.npy` persistence.
- Documented slope-window clipping, channel getter units, sensitivity-analysis
  side effects, quartile selection and compiled-helper fallback contracts.
- Added missing class, constructor, method and private-helper docstrings.

## Behaviour decisions requiring confirmation

The documentation records existing behaviour; no executable code was changed.
These discrepancies need separate implementation work and appropriate tests.

| API | Finding | Proposed decision |
| --- | --- | --- |
| `HCurve.getKurtosis/getSkewness` and density variants | Getter indices are swapped relative to `_get_moments` | Swap indices; changes existing returned values |
| `HCurve._calculate_hi2` | The last trapezoidal interval is omitted | Include every interval; changes HI2 |
| `Grid(path, band)` / `DEM(path, band)` | `band` is ignored; band 1 is always read | Honour the requested band |
| `Channel.getXY(head)` | `head=False` is ignored | Reverse coordinates consistently with the other getters |
| `Flow.snapPoints` / `Network.snapPoints` | XY-only deduplication compares Y alone and drops Y; extra-column inputs have a different output contract | Deduplicate by target ID, preserve XY and extra columns, decide whether to retain the appended ID when duplicates are allowed |
| `Basin.copy` | Calls `Basin()` without its mandatory `dem` argument | Create a copy without invoking the path-loading constructor |
| `Flow.flowAccumulation(weights)` | Resampling branch calls nonexistent `get_extent` and `get_ncells` and allocates integer weights | Decide whether to support aligned-only weights or repair resampling with float weights |
| `Channel.save/load` | Saves regression ID/start pairs but reloads them as start/end | Store and reload actual bounds, with a decision on legacy-file compatibility |
| `Channel.addRegression/getRegression` | Fit excludes `p2`; coverage lookup includes it, and `p2 == channel size` is rejected | Define inclusive or exclusive end consistently |
| `Channel.getRegression` | `min_dist` is never updated; overlapping fits can return the last qualifying fit rather than the nearest midpoint | Update distance when selecting a candidate |
| `BNetwork.chiPlot(relative=True)` | Recalculates default-theta chi after channels are extracted; plotted chi is not explicitly offset | Define a side-effect-free relative plot and its chi origin |
| `BNetwork.chiSensitivityAnalysis` | Leaves basin chi/theta at last candidate, not best candidate | Restore original state or keep the best fit explicitly |
| `SwathProfile` sampling | Coordinates use fractions `i/n_points`, while distances use `arange`; at exact step multiples their lengths differ | Define one shared sampling-distance array |
| `SwathProfile._get_zi` | Valid zero elevations become NaN | Preserve zero unless it is the actual NoData sentinel |
| `SwathProfile._get_parameters` | All-NaN rows raise; zero relief divides by zero | Define outputs for missing/flat rows |
| `Flow.drainageBasins` | No qualifying outlets falls back to all basins; a single supplied custom ID becomes 1 | Define empty-result and single-ID behaviour |
| `Flow.save/load` | Negative elevations / uint32 limits and float32 decoding can lose information | Decide whether a new backward-compatible persistence format is required |
| `PRaster.xyToCell/isInside` | Truncation toward zero can map just-outside coordinates into edge cells; rotation is ignored | Define boundary rounding and rotated-raster support |
| `Grid.isInside` | Vector NoData path indexes input lists with NumPy masks | Normalize vector inputs, or formally restrict them to arrays |

## Validation

- Parsed and compiled all library source files.
- Compared executable ASTs before/after after removing docstrings: unchanged.
- Verified existing Usage sections and repository examples/tutorials unchanged.
- Checked `git diff --check`.
- Runtime tests were not run: this session lacks GDAL, Numba and Shapely.
  These edits change documentation only; no numerical-equivalence claim is made.
