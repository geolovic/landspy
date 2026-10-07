"""Sample Linux process RSS by Flow construction phase, outside the worker.

python benchmarks/profile_flow_memory.py --size 14480 --output phases.json
Sampling uses a separate process, so C extensions holding the GIL cannot
block the sampler. Phase times come from the worker's own clock.
"""

import argparse
import hashlib
import json
import os
import queue
import resource
import subprocess
import sys
import threading
import time


def worker(args):
    import numpy as np
    from landspy import DEM, Flow

    def emit(event, **fields):
        print(json.dumps(dict(event=event, timestamp=time.perf_counter(), **fields)), flush=True)

    warm = DEM()
    warm.setArray(np.ones((3, 3), dtype='float32'))
    warm.fill()
    dem = DEM()
    dem._array = np.random.default_rng(38).integers(
        0, 2000, (args.size, args.size), dtype='int16').astype('float32')
    dem._size = (args.size, args.size)
    dem._tipo = 'float32'

    starts = {'Filling DEM ...': 'fill',
              'Identifiying flats and sills ...': 'flats_and_sills',
              'Identifiying presills ...': 'presills',
              'Generating auxiliar topography ...': 'aux_topography',
              'Calculating weights ...': 'weights',
              'Sorting pixels ...': 'sorting',
              'Calculating receivers ...': 'receivers'}

    def progress(message):
        if message in starts:
            emit('phase', name=starts[message])
        elif message.startswith('3/7'):
            emit('phase', name='aux_topography')
        elif message.startswith('5/7'):
            emit('phase', name='cleanup_after_weights')
        elif message == 'Finishing ...':
            emit('phase', name='filtering')
        elif message == 'Flow algorithm successfully completed':
            emit('phase', name='elevations')

    original = Flow._get_nodata_pos

    def nodata_positions(self):
        emit('phase', name='nodata_positions')
        return original(self)

    Flow._get_nodata_pos = nodata_positions
    emit('phase', name='constructor_setup')
    start = time.perf_counter()
    flow = Flow(dem, verbose=True, verb_func=progress)
    elapsed = time.perf_counter() - start
    emit('finished')
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    checksums = {name: hashlib.sha256(memoryview(getattr(flow, name))).hexdigest()
                 for name in ('_ix', '_ixc', '_zx')}
    emit('result', input_MiB=dem.readArray().nbytes / 1024**2,
         seconds=elapsed, process_peak_MiB=peak, checksums=checksums)


def profile(args):
    if not sys.platform.startswith('linux'):
        raise RuntimeError('This sampler requires Linux /proc')
    events = queue.Queue()
    command = [sys.executable, os.path.abspath(__file__), '--worker',
               '--size', str(args.size)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)

    def read_events():
        for line in process.stdout:
            events.put(json.loads(line))

    reader = threading.Thread(target=read_events, daemon=True)
    reader.start()
    rows = []
    current = None
    result = None
    page_mib = os.sysconf('SC_PAGE_SIZE') / 1024**2

    def rss():
        try:
            with open('/proc/{}/statm'.format(process.pid)) as stream:
                return int(stream.read().split()[1]) * page_mib
        except (FileNotFoundError, ProcessLookupError):
            return 0

    def close_phase(timestamp, memory):
        if current is not None:
            current['seconds'] = round(timestamp - current.pop('start'), 3)
            current['rss_exit_MiB'] = round(memory, 1)
            current['sampled_peak_MiB'] = round(current['sampled_peak_MiB'], 1)
            rows.append(current)
            print(json.dumps(current), flush=True)

    while process.poll() is None or reader.is_alive() or not events.empty():
        memory = rss()
        if current is not None:
            current['sampled_peak_MiB'] = max(current['sampled_peak_MiB'], memory)
        while not events.empty():
            event = events.get()
            if event['event'] in ('phase', 'finished'):
                close_phase(event['timestamp'], memory)
                current = None
                if event['event'] == 'phase':
                    current = dict(phase=event['name'], start=event['timestamp'],
                                   rss_entry_MiB=round(memory, 1), sampled_peak_MiB=memory)
            elif event['event'] == 'result':
                result = {key: value for key, value in event.items()
                          if key not in ('event', 'timestamp')}
        # This is a short sampling interval, not a wait for task completion.
        time.sleep(args.interval)
    reader.join()
    if process.wait() != 0 or result is None:
        raise RuntimeError('Flow worker failed or did not return a result')
    result.update(phases=rows, sample_interval_seconds=args.interval,
                  size=args.size, dtype='float32', seed=38)
    if args.output:
        with open(args.output, 'w') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    print(json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=2048)
    parser.add_argument('--interval', type=float, default=0.05)
    parser.add_argument('--output')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.size < 1 or args.interval <= 0:
        parser.error('size and interval must be positive')
    worker(args) if args.worker else profile(args)


if __name__ == '__main__':
    main()
