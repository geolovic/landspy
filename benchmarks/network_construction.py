"""Measure Network construction in a separate process (requires psutil).

Prepare the exact Flow once, without GeoTIFF elevation quantization:
  python benchmarks/network_construction.py --cache /tmp/net-flow --prepare DEM.tif
Then compare isolated processes against a local Git revision:
  python benchmarks/network_construction.py --cache /tmp/net-flow --output /tmp/before --reference-ref 82ffe00
  python benchmarks/network_construction.py --cache /tmp/net-flow --output /tmp/after --warmup

Defaults: threshold=1500, gradients=False. RSS includes the loaded Flow and
Python/native runtimes. --warmup excludes initial JIT compilation. No DEM or
Flow cache is committed; allow several GiB of disk and RAM for large DEMs.
"""

import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time


def event(phase, state, **extra):
    import psutil
    print(json.dumps(dict(phase=phase, state=state, time=time.perf_counter(),
                          rss_bytes=psutil.Process().memory_info().rss, **extra)), flush=True)


def worker(args):
    import numpy as np
    from landspy import DEM, Flow, Network

    cache = Path(args.cache)
    if args.prepare:
        cache.mkdir(parents=True, exist_ok=True)
        if any(cache.iterdir()):
            raise ValueError('Use an empty cache directory when preparing a Flow')
        event('prepare', 'start')
        flow = Flow(DEM(args.prepare))
        for name in ('_ix', '_ixc', '_zx', '_nodata_pos'):
            np.save(cache / (name + '.npy'), getattr(flow, name))
        (cache / 'layout.json').write_text(json.dumps(dict(size=flow._size, geot=flow._geot, proj=flow._proj)))
        event('prepare', 'end', edges=int(flow._ix.size))
        return

    if args.warmup:
        dem = DEM()
        dem.setArray(np.array([[4., 3., 2.], [3., 2., 1.], [2., 1., 0.]]))
        Network(Flow(dem), threshold=type(args.threshold)(1))
    event('load', 'start')
    flow = Flow()
    layout = json.loads((cache / 'layout.json').read_text())
    flow._size, flow._geot, flow._proj = tuple(layout['size']), tuple(layout['geot']), layout['proj']
    for name in ('_ix', '_ixc', '_zx', '_nodata_pos'):
        setattr(flow, name, np.load(cache / (name + '.npy')))
    event('load', 'end', edges=int(flow._ix.size))

    original_accumulation = flow.flowAccumulation
    def measured_accumulation(*a, **kw):
        event('accumulation', 'start')
        result = original_accumulation(*a, **kw)
        event('accumulation', 'end')
        return result
    flow.flowAccumulation = measured_accumulation
    original_chi = Network.calculateChi
    def measured_chi(self, *a, **kw):
        event('chi', 'start')
        result = original_chi(self, *a, **kw)
        event('chi', 'end')
        return result
    Network.calculateChi = measured_chi
    event('Network', 'start')
    net = Network(flow, threshold=args.threshold, gradients=args.gradients)
    event('Network', 'end', nodes=int(net._ix.size), threshold=net._threshold)
    arrays = ('_ix', '_ixc', '_ax', '_zx', '_dd', '_dx', '_chi', '_slp', '_ksn', '_r2slp', '_r2ksn')
    np.savez(Path(args.output) / 'network.npz', **{name: getattr(net, name) for name in arrays})


def monitor(args):
    import psutil

    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg')
    with tempfile.TemporaryDirectory(prefix='landspy-network-') as temporary:
        if args.reference_ref:
            repo = Path(__file__).resolve().parents[1]
            archive = subprocess.check_output(['git', '-C', str(repo), 'archive', args.reference_ref, 'src'])
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(temporary, filter='data')
            env['PYTHONPATH'] = str(Path(temporary) / 'src')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--cache', str(Path(args.cache).resolve()), '--output', str(output),
                   '--threshold', str(args.threshold)]
        if args.prepare:
            command += ['--prepare', str(Path(args.prepare).resolve())]
        if args.warmup:
            command.append('--warmup')
        if args.gradients:
            command.append('--gradients')
        samples = []
        with (output / 'events.log').open('w') as log:
            child = subprocess.Popen(command, stdout=log, stderr=log, env=env)
            process = psutil.Process(child.pid)
            while child.poll() is None:
                try:
                    samples.append((time.perf_counter(), process.memory_info().rss))
                except psutil.NoSuchProcess:
                    break
                time.sleep(.01)
            code = child.wait()
    events = []
    for line in (output / 'events.log').read_text().splitlines():
        try:
            events.append(json.loads(line))
        except ValueError:
            pass
    summary = dict(exit_code=code, reference_ref=args.reference_ref,
                   warmup=args.warmup, gradients=args.gradients, phases={})
    for end in events:
        if end.get('state') != 'end':
            continue
        start = next(e for e in events if e['phase'] == end['phase'] and e['state'] == 'start')
        rss = [r for t, r in samples if start['time'] <= t <= end['time']]
        summary['phases'][end['phase']] = dict(seconds=end['time'] - start['time'],
            peak_rss_bytes=max([start['rss_bytes'], end['rss_bytes'], *rss]),
            start_rss_bytes=start['rss_bytes'])
        summary['phases'][end['phase']].update({k: v for k, v in end.items() if k not in ('phase', 'state', 'time', 'rss_bytes')})
    (output / 'metrics.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if code:
        print((output / 'events.log').read_text(), file=sys.stderr)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--output', default='network-benchmark-output')
    parser.add_argument('--prepare', metavar='DEM')
    parser.add_argument('--reference-ref', help='Read a local Git revision without changing the checkout')
    parser.add_argument('--threshold', type=float, default=1500)
    parser.add_argument('--warmup', action='store_true')
    parser.add_argument('--gradients', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    options = parser.parse_args()
    if options.worker:
        worker(options)
    else:
        sys.exit(monitor(options))
