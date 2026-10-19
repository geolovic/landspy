# Changelog

## 1.5

- Preserve DEM NoData masks during filling and handle unsupported extended precision inputs.
- Accelerate Flow accumulation with compiled topological traversal and reduce temporary memory.
- Optimize Network construction, distances, chi and gradient calculations using compact network topology and compiled kernels.
- Accelerate drainage basin delineation, including point outlets and channel snapping workflows.
- Reduce BNetwork construction memory and optimize complete and individual channel extraction while preserving channel ordering and connections.
- Add numerical regression tests and reproducible time and memory benchmarks.

## 1.4.1

Previous published release.
