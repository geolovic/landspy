"""Measure drainageBasins on a cached Flow, with external RSS sampling.

Last modified: 2026-10-16 11:38:32 CEST (Europe/Madrid).
Requires psutil. Prepare the Flow cache with network_construction.py first.
The JSON report and optional labels are written outside the timed interval.
"""

import argparse
import hashlib
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
    from landspy import Flow

    if args.warmup:
        small = Flow()
        small._size = (2, 2)
        small._ix = np.array([0, 1, 2], dtype=np.uint32)
        small._ixc = np.array([1, 3, 3], dtype=np.uint32)
        small.drainageBasins(min_area=0, asgrid=False)
        small.drainageBasins(asgrid=False)
    cache = Path(args.cache)
    event('load', 'start')
    flow = Flow()
    layout = json.loads((cache / 'layout.json').read_text())
    flow._size, flow._geot, flow._proj = tuple(layout['size']), tuple(layout['geot']), layout['proj']
    for name in ('_ix', '_ixc', '_zx', '_nodata_pos'):
        setattr(flow, name, np.load(cache / (name + '.npy')))
    event('load', 'end', edges=int(flow._ix.size))
    kwargs = dict(asgrid=False, min_area=args.min_area)
    if args.mode == 'all':
        kwargs['min_area'] = 0
    elif args.mode == 'seeded':
        indices = flow._ixc[np.linspace(0, flow._ixc.size - 1, 10, dtype=int)]
        x, y = flow.cellToXY(*flow.indToCell(indices))
        kwargs['outlets'] = np.column_stack((x, y, np.arange(11, 21)))
    original_accumulation = flow.flowAccumulation
    def measured_accumulation(*a, **kw):
        event('accumulation', 'start')
        result = original_accumulation(*a, **kw)
        event('accumulation', 'end')
        return result
    flow.flowAccumulation = measured_accumulation
    event('drainageBasins', 'start', mode=args.mode)
    labels = flow.drainageBasins(**kwargs)
    event('drainageBasins', 'end', mode=args.mode)
    flat = labels.ravel()
    if args.compare:
        reference = np.load(args.compare, mmap_mode='r')
        assert reference.shape == labels.shape
        assert reference.dtype == labels.dtype
        reference = reference.ravel()
        for start in range(0, flat.size, 1000000):
            np.testing.assert_array_equal(flat[start:start + 1000000], reference[start:start + 1000000])
    if args.save_labels:
        np.save(Path(args.output) / 'labels.npy', labels)
    event('result', 'verified', sha256=hashlib.sha256(memoryview(flat).cast('B')).hexdigest(),
          dtype=str(labels.dtype), shape=labels.shape, max_id=int(labels.max()),
          compared=bool(args.compare))


def monitor(args):
    import psutil

    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg')
    with tempfile.TemporaryDirectory(prefix='landspy-basins-') as temporary:
        if args.reference_ref:
            repo = Path(__file__).resolve().parents[1]
            archive = subprocess.check_output(['git', '-C', str(repo), 'archive', args.reference_ref, 'src'])
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(temporary, filter='data')
            env['PYTHONPATH'] = str(Path(temporary) / 'src')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--cache', str(Path(args.cache).resolve()), '--output', str(output),
                   '--mode', args.mode, '--min-area', str(args.min_area)]
        if args.compare:
            command += ['--compare', str(Path(args.compare).resolve())]
        for option in ('warmup', 'save_labels'):
            if getattr(args, option):
                command.append('--' + option.replace('_', '-'))
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
                   min_area=args.min_area, warmup=args.warmup, phases={})
    for end in events:
        if end.get('state') == 'verified':
            summary['result'] = {k: v for k, v in end.items()
                                 if k not in ('phase', 'state', 'time', 'rss_bytes')}
        if end.get('state') != 'end':
            continue
        start = next(e for e in events if e['phase'] == end['phase'] and e['state'] == 'start')
        rss = [r for t, r in samples if start['time'] <= t <= end['time']]
        summary['phases'][end['phase']] = dict(seconds=end['time'] - start['time'],
            peak_rss_bytes=max([start['rss_bytes'], end['rss_bytes'], *rss]), start_rss_bytes=start['rss_bytes'])
    (output / 'metrics.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if code:
        print((output / 'events.log').read_text(), file=sys.stderr)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--mode', choices=('all', 'default', 'seeded'), default='default')
    parser.add_argument('--min-area', type=float, default=.005)
    parser.add_argument('--reference-ref')
    parser.add_argument('--warmup', action='store_true')
    parser.add_argument('--save-labels', action='store_true')
    parser.add_argument('--compare', help='Compare every label with a saved .npy baseline')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    options = parser.parse_args()
    if options.worker:
        worker(options)
    else:
        sys.exit(monitor(options))
