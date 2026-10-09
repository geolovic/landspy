"""Compare gradient arrays against an earlier Network implementation."""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from landspy import DEM, Flow, Network


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', default='6958244af1b90ee78b328ba6245b8d305c6bf3e8')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    comparisons = []
    source = subprocess.check_output(['git', 'show', f'{args.baseline}:src/landspy/network.py'], cwd=repo)
    with tempfile.TemporaryDirectory(prefix='network-regression-') as tmp:
        path = Path(tmp) / 'network.py'
        path.write_bytes(source)
        spec = importlib.util.spec_from_file_location('landspy._reference_network', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for name in ['small25', 'tunez', 'jebja30']:
            fresh = Flow(DEM(str(repo / 'tests/data/in' / f'{name}.tif')))
            path = Path(tmp) / f'{name}_flow.tif'
            fresh.save(str(path))
            for mode, flow in [('fresh', fresh), ('loaded', Flow(str(path)))]:
                for theta in [.25, .45, .65]:
                    before = module.Network(flow, thetaref=theta)
                    after = Network(flow, thetaref=theta)
                    for npoints in [1, 2, 5, 10, 20]:
                        for kind in ['slp', 'ksn']:
                            before.calculateGradients(npoints, kind)
                            after.calculateGradients(npoints, kind)
                            for field in ['_' + kind, '_r2' + kind]:
                                a = getattr(before, field)
                                b = getattr(after, field)
                                assert a.dtype == b.dtype, field
                                np.testing.assert_allclose(b, a, rtol=1e-10, atol=1e-12)
                                comparisons.append({'dem': name, 'flow': mode,
                                    'theta': theta, 'npoints': npoints, 'field': field,
                                    'max_absolute_difference': float(np.max(np.abs(a-b), initial=0.))})
                        for field in ['_ix', '_ixc', '_dd', '_dx', '_chi', '_zx', '_ax']:
                            np.testing.assert_array_equal(getattr(before, field), getattr(after, field))
                print(name, mode, 'passed', flush=True)
    result = {'baseline': args.baseline, 'rtol': 1e-10, 'atol': 1e-12,
              'comparisons': comparisons}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(len(comparisons), 'gradient array comparisons passed')


if __name__ == '__main__':
    main()
