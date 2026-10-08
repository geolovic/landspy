# Changelog

## 1.4.0

- Replace DEM filling with Numba-compiled Priority-Flood.
- Fill only once during Flow construction and reuse filled elevations.
- Replace MCP weight calculation with compiled Dijkstra and compile receiver selection.
- Calculate unit-cost distances directly from the flat mask in float32; retain
  float64 for variable costs and stored elevations. This deliberate precision
  change can alter a small number of receivers compared with 1.3.x.
- Reduce temporary arrays and use stable lexicographic DEM ordering.
- On the supplied real DEM, reduce peak Flow construction memory from 19.51
  to 7.67 GiB and construction time from 218.8 to 127.7 seconds. These are
  single-run measurements, excluding input reading and initial JIT compilation.
- Require Python 3.10 or newer, Shapely 2 or newer, and Numba.
- Exclude benchmarks and test data from PyPI distributions.
