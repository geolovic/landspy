"""Compare accumulation or Network construction using lossless saved Flow arrays.

Pass --dem for raster metadata and --flow-arrays for a directory containing
_ix.npy, _ixc.npy, _zx.npy and _nodata_pos.npy. No raster resampling or Flow
reconstruction occurs. Input loading and JIT warmup are outside the timer;
peak RSS includes them. Each measurement runs in a fresh process.
"""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time
import types

import numpy as np
from osgeo import gdal
from landspy import Flow, Network
from landspy._accumulation import accumulate_flow
from landspy._network import accumulate_downstream, linear_fit


def worker(args):
    cls = Flow
    if args.revision:
        source = subprocess.check_output(
            ['git', 'show', args.revision + ':src/landspy/flow.py'],
            cwd=Path(__file__).resolve().parents[1], text=True)
        module = types.ModuleType('landspy._accumulation_reference')
        module.__package__ = 'landspy'
        exec(compile(source, '<reference-flow>', 'exec'), module.__dict__)
        cls = module.Flow
    dataset = gdal.Open(str(args.dem))
    flow = cls()
    flow._size = (dataset.RasterXSize, dataset.RasterYSize)
    flow._geot = dataset.GetGeoTransform()
    flow._proj = dataset.GetProjection()
    dataset = None
    for name in ['_ix', '_ixc', '_zx', '_nodata_pos']:
        setattr(flow, name, np.load(args.flow_arrays / (name + '.npy')))
    accumulate_flow(np.array([0], np.uint32), np.array([1], np.uint32), np.ones(2, np.uint32))
    accumulate_downstream(np.array([0]), np.array([1]), np.array([1.]), 2)
    linear_fit(np.array([0., 1., 2.]), np.array([0., 1., 2.]))
    start = time.perf_counter()
    if args.network:
        network = Network(flow, threshold=args.threshold, gradients=args.gradients)
        names = ['_ix', '_ixc', '_ax', '_zx', '_dd', '_dx', '_chi',
                 '_slp', '_ksn', '_r2slp', '_r2ksn']
        arrays = {name: getattr(network, name) for name in names}
    else:
        arrays = {'accumulation': flow.flowAccumulation(nodata=args.nodata, asgrid=False)}
    seconds = time.perf_counter() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    # Capture peak before output serialization and correctness checks.
    np.savez(args.arrays, **arrays)
    return dict(revision=args.revision or 'working-tree', seconds=seconds,
                peak_rss_mib=peak, dem_cells=flow.getNCells(),
                flow_edges=flow._ix.size, nodata_cells=flow._nodata_pos.size,
                mode='network' if args.network else 'accumulation',
                nodata=args.nodata if not args.network else False,
                threshold=args.threshold if args.network else None,
                gradients=args.gradients if args.network else None,
                hashes={name: hashlib.sha256(memoryview(array)).hexdigest()
                        for name, array in arrays.items()},
                dtypes={name: str(array.dtype) for name, array in arrays.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dem', type=Path, required=True)
    parser.add_argument('--flow-arrays', type=Path, required=True)
    parser.add_argument('--baseline', default='HEAD')
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--network', action='store_true')
    parser.add_argument('--gradients', action='store_true')
    parser.add_argument('--threshold', type=int, default=0)
    parser.add_argument('--nodata', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--revision', help=argparse.SUPPRESS)
    parser.add_argument('--arrays', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)), flush=True)
        return
    if args.runs < 1:
        parser.error('--runs must be positive')
    repo = Path(__file__).resolve().parents[1]
    baseline = subprocess.check_output(
        ['git', 'rev-parse', args.baseline + '^{commit}'], cwd=repo, text=True).strip()
    common = [sys.executable, str(Path(__file__).resolve()), '--worker',
              '--dem', str(args.dem.resolve()), '--flow-arrays', str(args.flow_arrays.resolve()),
              '--threshold', str(args.threshold)]
    for flag in ['network', 'gradients', 'nodata']:
        if getattr(args, flag):
            common.append('--' + flag)
    results = []
    with tempfile.TemporaryDirectory(prefix='landspy-accumulation-') as temporary:
        for run in range(args.runs):
            paths = []
            for index, revision in enumerate([baseline, None]):
                path = Path(temporary) / f'{index}.npz'
                paths.append(path)
                command = common + ['--arrays', str(path)]
                if revision:
                    command += ['--revision', revision]
                row = json.loads(subprocess.check_output(command, text=True))
                row['run'] = run
                results.append(row)
                print(json.dumps(row), flush=True)
            with np.load(paths[0]) as before, np.load(paths[1]) as after:
                for name in before.files:
                    a, b = before[name].ravel(), after[name].ravel()
                    assert a.dtype == b.dtype and a.shape == b.shape, name
                    # Bound comparison temporaries on large rasters.
                    for start in range(0, a.size, 1_000_000):
                        np.testing.assert_array_equal(a[start:start+1_000_000], b[start:start+1_000_000])
                    assert results[-2]['hashes'][name] == results[-1]['hashes'][name], name
            for path in paths:
                path.unlink()
            if args.output:
                args.output.write_text(json.dumps(dict(baseline=baseline, exact_comparison=True, runs=results), indent=2) + '\n')
    print('Exact arrays, dtypes and hashes matched in every run.', flush=True)


if __name__ == '__main__':
    main()
