"""Compare BNetwork construction and channel extraction with sampled process RSS.

Last modified: 2026-10-19 11:56:48 CEST (Europe/Madrid).
Use --reference-ref HEAD for the committed baseline and omit it for working code.
Inputs: cached Network .npz, Flow layout.json, and a labelled basin raster.
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


def worker(args):
    import numpy as np
    import psutil
    from landspy import Network, BNetwork, Grid

    # Exclude Numba compilation/cache loading from the measured phases.
    try:
        from landspy._network_kernels import basin_channel_path
        basin_channel_path(0, np.array([1, 2]), np.array([0, 1, -1]),
                           np.full(3, -1, dtype=np.int64), 0)
    except ImportError:
        pass

    def event(phase, state, **extra):
        print(json.dumps(dict(phase=phase, state=state, time=time.perf_counter(),
                              rss_bytes=psutil.Process().memory_info().rss,
                              **extra)), flush=True)

    net = Network()
    layout = json.loads(Path(args.layout).read_text())
    net._size, net._geot, net._proj = tuple(layout['size']), tuple(layout['geot']), layout['proj']
    with np.load(args.network) as archive:
        for name in archive.files:
            setattr(net, name, archive[name])
    grid = Grid(args.basins)
    event('constructor', 'start')
    basin = BNetwork(net, grid, bid=args.bid)
    event('constructor', 'end', nodes=int(basin._ix.size), raster_cells=int(basin.getNCells()))
    event('all_channels', 'start')
    channels = basin.getChannels('ALL')
    event('all_channels', 'end', channels=len(channels))
    digest = hashlib.sha256()

    def update(value):
        if isinstance(value, np.ndarray):
            digest.update(str(value.shape).encode())
            digest.update(np.ascontiguousarray(value).tobytes())
        else:
            digest.update(json.dumps(value, default=lambda v: v.item(), sort_keys=True).encode())

    for obj in [basin, *channels]:
        for name in sorted(obj.__dict__):
            update(getattr(obj, name))
    first = channels[0]
    del channels
    event('single_channel', 'start')
    channel = basin.getChannel(0)
    event('single_channel', 'end')
    for name in first.__dict__:
        np.testing.assert_equal(getattr(first, name), getattr(channel, name))
    print(json.dumps(dict(digest=digest.hexdigest())), flush=True)


def monitor(args):
    import psutil
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg')
    with tempfile.TemporaryDirectory(prefix='landspy-bnetwork-') as temporary:
        if args.reference_ref:
            repo = Path(__file__).resolve().parents[1]
            archive = subprocess.check_output(['git', '-C', str(repo), 'archive', args.reference_ref, 'src'])
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(temporary, filter='data')
            env['PYTHONPATH'] = str(Path(temporary) / 'src')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--network', str(Path(args.network).resolve()), '--layout', str(Path(args.layout).resolve()),
                   '--basins', str(Path(args.basins).resolve()), '--bid', str(args.bid), '--output', str(output)]
        samples = []
        with (output / 'events.log').open('w') as log:
            child = subprocess.Popen(command, stdout=log, stderr=log, env=env)
            process = psutil.Process(child.pid)
            while child.poll() is None:
                try:
                    samples.append((time.perf_counter(), process.memory_info().rss))
                except psutil.NoSuchProcess:
                    break
                time.sleep(0.01)
            code = child.wait()
        events = []
        for line in (output / 'events.log').read_text().splitlines():
            try:
                events.append(json.loads(line))
            except json.JSONDecodeError:
                pass
        result = dict(exit_code=code, bid=args.bid, reference_ref=args.reference_ref,
                      peak_rss_bytes=max((rss for _, rss in samples), default=0), phases={})
        for phase in ('constructor', 'all_channels', 'single_channel'):
            records = [e for e in events if e.get('phase') == phase]
            if len(records) != 2:
                continue
            start, end = records
            result['phases'][phase] = dict(seconds=end['time']-start['time'],
                start_rss_bytes=start['rss_bytes'],
                peak_rss_bytes=max([start['rss_bytes'], end['rss_bytes']] +
                    [rss for moment, rss in samples if start['time'] <= moment <= end['time']]),
                **{k:v for k,v in end.items() if k not in ('phase', 'state', 'time', 'rss_bytes')})
        for e in events:
            if 'digest' in e:
                result['digest'] = e['digest']
        (output / 'metrics.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2))
        if code:
            raise SystemExit(code)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('network', 'layout', 'basins', 'output'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--bid', type=int, default=1)
    parser.add_argument('--reference-ref')
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    worker(args) if args.worker else monitor(args)
