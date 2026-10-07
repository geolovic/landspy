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
