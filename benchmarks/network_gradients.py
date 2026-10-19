"""Benchmark Network gradients with external RSS sampling and reference checks.

Last modified: 2026-10-19 10:46:34 CEST (Europe/Madrid).
Requires psutil. Inputs are the Network .npz and Flow layout.json produced by
network_construction.py. Reading/writing/comparison are outside the timed phase.
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


def calculate(net, npoints, mode, measured=True):
    if mode == 'both' and hasattr(net, '_calculate_gradients'):
        net._calculate_gradients(npoints, 'slp', both=True)
    else:
        for kind in ('slp', 'ksn'):
            if measured:
                event(kind, 'start')
            net.calculateGradients(npoints, kind)
            if measured:
                event(kind, 'end')


def worker(args):
    import numpy as np
    from landspy import Network

    if args.warmup:
        small = Network()
        small._size = (4, 1)
        small._ix = np.array([0, 1, 2], dtype=np.uint32)
        small._ixc = np.array([1, 2, 3], dtype=np.uint32)
        small._zx = np.array([3., 2., 1.])
        small._dx, small._chi = small._zx.copy(), small._zx.copy()
        calculate(small, args.npoints, args.mode, measured=False)
    net = Network()
    layout = json.loads(Path(args.layout).read_text())
    net._size, net._geot, net._proj = tuple(layout['size']), tuple(layout['geot']), layout['proj']
    with np.load(args.network) as archive:
        for name in archive.files:
            setattr(net, name, archive[name])
    event('gradients', 'start')
    calculate(net, args.npoints, args.mode)
    event('gradients', 'end', nodes=int(net._ix.size), raster_cells=int(net.getNCells()))
    names = ('_slp', '_ksn', '_r2slp', '_r2ksn')
    if args.compare:
        errors = {}
        with np.load(args.compare) as reference:
            for name in names:
                actual, expected = getattr(net, name), reference[name]
                np.testing.assert_array_equal(actual == 0, expected == 0)
                np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-11)
                finite = np.isfinite(actual) & np.isfinite(expected)
                nonzero = finite & (expected != 0)
                errors[name] = dict(max_abs_error=float(np.max(np.abs(actual[finite] - expected[finite]), initial=0)),
                    max_relative_error=float(np.max(np.abs((actual[nonzero] - expected[nonzero]) / expected[nonzero]), initial=0)),
                    zero_mask_mismatches=0)
        event('comparison', 'verified', errors=errors)
    np.savez(Path(args.output) / 'gradients.npz', **{name: getattr(net, name) for name in names})


def monitor(args):
    import psutil

    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg')
    with tempfile.TemporaryDirectory(prefix='landspy-gradients-') as temporary:
        if args.reference_ref:
            repo = Path(__file__).resolve().parents[1]
            archive = subprocess.check_output(['git', '-C', str(repo), 'archive', args.reference_ref, 'src'])
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(temporary, filter='data')
            env['PYTHONPATH'] = str(Path(temporary) / 'src')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--network', str(Path(args.network).resolve()), '--layout', str(Path(args.layout).resolve()),
                   '--output', str(output), '--mode', args.mode, '--npoints', str(args.npoints)]
        if args.warmup:
            command.append('--warmup')
        if args.compare:
            command += ['--compare', str(Path(args.compare).resolve())]
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
    summary = dict(exit_code=code, reference_ref=args.reference_ref, mode=args.mode,
                   npoints=args.npoints, warmup=args.warmup, phases={})
    for end in events:
        if end.get('state') == 'verified':
            summary['comparison'] = end['errors']
        if end.get('state') != 'end':
            continue
        start = next(e for e in events if e['phase'] == end['phase'] and e['state'] == 'start')
        rss = [r for t, r in samples if start['time'] <= t <= end['time']]
        summary['phases'][end['phase']] = dict(seconds=end['time'] - start['time'],
            peak_rss_bytes=max([start['rss_bytes'], end['rss_bytes'], *rss]), start_rss_bytes=start['rss_bytes'])
        summary['phases'][end['phase']].update({k: v for k, v in end.items()
                                              if k not in ('phase', 'state', 'time', 'rss_bytes')})
    (output / 'metrics.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if code:
        print((output / 'events.log').read_text(), file=sys.stderr)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--network', required=True)
    parser.add_argument('--layout', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference-ref')
    parser.add_argument('--compare')
    parser.add_argument('--npoints', type=int, default=5)
    parser.add_argument('--mode', choices=('both', 'separate'), default='both')
    parser.add_argument('--warmup', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    options = parser.parse_args()
    if options.worker:
        worker(options)
    else:
        sys.exit(monitor(options))
