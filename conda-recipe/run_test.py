"""Exercise JIT kernels and Flow against a small DEM in the installed package."""
import numpy as np
from landspy import DEM, Flow

dem = DEM()
dem.setArray(np.array([[5, 5, 5], [5, 1, 5], [5, 5, 4]], dtype=np.int16))
filled = dem.fill()
assert filled.readArray()[1, 1] == 4
assert dem.readArray()[1, 1] == 1
for auxtopo in (False, True):
    flow = Flow(dem, auxtopo=auxtopo)
    assert flow._ix.size > 0
    assert flow._ix.shape == flow._ixc.shape == flow._zx.shape
    assert flow._ix.dtype == flow._ixc.dtype == np.uint32
    assert flow._zx.dtype == np.float64
    assert np.all(flow._ix != flow._ixc)
    assert np.isfinite(flow.flowAccumulation().readArray()).all()
