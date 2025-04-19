##############################################
# Automated ds::compare testing script for all implemented data structures.
# Stores the results as CSV file per data structure in a folder.
# Lots of folders, in fact.
#
# How to add your own structure to this script:
# 1. Make sure your data structure works as intended within the ds::compare framework.
#    Info: See the README for further details, but the tl;dr: It should be registered
#          in the cluster, the include.h file and the SCons file to compile correctly.
#          It should have its own DSType (see ds.h). All its exposed parameters should
#          be in both the DSArguments struct (see ds.h) and the clargs handler function
#          (see dscompare.cpp).
#
# 2. Add the exposed parameters of your data structure to the settings dict below.
#    Info: Verify the settings flags match the ones you specified in the clargs handler
#          function. If you used aliases, pick only one. The order of the data structures
#          (however not the settings) matters -- please make sure you place it in the
#          same order as in the ds.h DSType enum. If you don't have any exposed parameters,
#          please create an empty dict nevertheless. ('your_ds' = { })
#
# 3. Your data structure should now be registered.
#
# Feel free to do whatever you want to this script.
# Michael Eickmeyer, 2025 @ TU Wien.
##############################################

import subprocess, os, time, math, signal, sys, uuid, warnings, csv, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait
from threading import Event, Thread
from collections.abc import Iterable

##############################################
# A bunch of classes handling settings storage.
# Settings may either be fixed or fluid.
#
# Fixed settings (Value) consist of a single value, as well as its command line flag.
#
# Fluid settings (Range, Toggle) on the other hand allow for a parameter to increase,
# decrease or toggle within this test script. For instance, the Range class takes an
# optional 'func' parameter which allows the user to specify how exactly the value
# should increase within the range. The Toggle class allows the user to easily toggle
# between a truthy/falsy value.
#
# Feel free to extend the FluidSetting class for custom behavior if needed.
#
# All fluid settings consist of the functions...
# • get()   -> obtain the current value
# • flag()  -> obtain the command flag
# • reset() -> set the value back to its initial state
# • steps() -> obtain how many steps the setting needs from start to finish for completion
# • next()  -> change the value to its next state
# • final() -> check whether the value has reached its final state (upper bound, !current, etc.) 
##############################################

class Value:
    def __init__(self, flag, value):
        self._flag = f'-{flag}' if len(flag) == 1 else f'--{flag}'
        self._value = value

    def get(self, as_str=False):
        if not as_str:
            return self._value

        val_type = type(self._value)
        if (val_type is bool):
            return str(int(self._value))
        if (val_type is list):
            return [str(v) for v in self._value]

        return str(self._value)
    
    def flag(self):
        return self._flag

class FluidSetting(Value):
    def __init__(self, flag, value):
        super().__init__(flag, value)
        self._default = value

    def reset(self):
        self._value = self._default

    def steps(self):
        raise NotImplementedError

    def next(self):
        raise NotImplementedError

    def final(self):
        raise NotImplementedError

class Range(FluidSetting):
    def __init__(self, flag, start, end, func = None):
        super().__init__(flag, start)
        self._end = end
        
        if not func:
            func = lambda x: x + start
        self._func = func

    def steps(self):
        curr = self._default
        end = self._end
        fn = self._func

        s = 1
        while not fn(curr) > end:
            curr = fn(curr)
            s += 1
        return s
    
    def next(self):
        self._value = self._func(self._value)

    def final(self):
        next_value = self._func(self._value)
        return (next_value > self._end)
    
class Toggle(FluidSetting):
    def __init__(self, flag, state):
        super().__init__(flag, state)

    def steps(self):
        return 2
    
    def next(self):
        self._value = not self._value

    def final(self):
        return (self._value != self._default)

##############################################
# Base config.
##############################################

settings = {
    'testing': {
        # Enables multithreading for concurrent evaluation of environment maps.
        'multithreading': True,
        # (Max.) Number of batches the envmaps get divided into. Set to -1 to disable & use the provided folder structure.
        'batches': 8,
        # Metrics to store in the benchmark.csv files. Names must match the metrics specified in ds::compare.
        'metrics': ['MD', 'Memory', 'store() (s)']
    },
    'general': {
        # Path to folder containing the envmaps.
        'envmap_path': Value('p', './data/tests/envmaps/'),
        # Path to folder containing the results.
        'result_path': Value('rp', './data/results/'),
        # Number of learning samples per data structure. These samples are used to create a guiding distribution.
        'samples_learning': Range('sl', start=64, end=65536, func=lambda x: x * 2),
        # Number of guiding samples per data structure. These samples are used to recreate the sampled distribution from the guiding distribution.
        'samples_guiding': Value('sg', 65536),
        # A blacklist specifying which data structures should be skipped in the overall test.
        'blacklist': Value('b', [0]),
        # Whether to normalize the envmap each data structure is 'learning' with. (Keep this true unless you know what you're doing.)
        'normalize': Value('n', True)
    },
    'noise': {
        # Whether the initial envmap should be noisified.
        'envmap': Value('ne', False),
        # Whether a variable amount of noise should be introduced to each learning sample.
        'samples': Value('ns', False)
    },
    'structures': {
        'Unidirectional': { },
        'Spherical Harmonics': {
            'bands': Range('shb', start=1, end=10),
            'depth': Range('shd', start=1, end=20),
            'use_offset': Toggle('sho', True)
        },
        'D-Tree': {
            'frac_loss': Value('dtl', 'none'),
            'dir_filter': Value('dtf', 'nearest'),
            'threshold': Range('dtt', start=0.01, end=0.5, func=lambda x: x + 0.01),
            'iterations': Value('dti', -1),
            'max_depth': Range('dtd', start=2, end=20)
        },
        'Tile Coding': {
            'tilings': Range('t', start=1, end=8),
            'tiles_x': Range('tx', start=2, end=32),
            'tiles_y': Range('ty', start=2, end=32)
        },
        'Binary Tile Coding': {
            'tilings': Range('bt', start=1, end=8),
            'tiles_x': Range('btx', start=1, end=8),
            'tiles_y': Range('bty', start=1, end=8),
            'max_depth': Range('btd', start=2, end=20),
            'subdiv_thresh': Range('btt', start=0.0001, end=0.01, func=lambda x: x * 2)
        },
        'von Mises-Fisher Mixtures': {
            'components': Range('vc', start=1, end=32),
            'use_ruppert': Toggle('vr', True)
        }
    }
}

##############################################
# Stat Tracker.
##############################################

class StatTrak:
    run = {
        'curr': 0,
        'total': 0,
        'ds': 0
    }
    time = {
        'avg': 0,
        'count': 0
    }
    progress = { }

    def reset(self):
        self.time['avg'] = 0
        self.time['count'] = 0
        self.run['curr'] = 0

    def increment(self, what: str):
        if what == 'run':
            self.run['curr'] += 1
        else:
            self.time['count'] += 1

    def store_time(self, elapsed):
        avg_delta = elapsed - self.time['avg']
        self.time['avg'] += avg_delta / self.time['count']

##############################################
# Start of code.
##############################################

stats = StatTrak()
stop_event = Event()
pause_event = Event()

timestamp = time.time()
base_pairing = {}
batch_paths = []
commands = []
ds_keys = list(settings['structures'].keys())

def build_argvals(s, a):
    """
    Builds a list of setting values required for cl call.

    Parameters
    ----------
    a : list
      List for accumulating all relevant (flag, value) pairs
    s : obj
      The settings object
    """

    for k, v in s.items():
        if isinstance(v, dict):
            build_argvals(v, a)
        elif k.lower() in ['multithreading', 'batches', 'metrics', 'envmap_path', 'result_path']:
            continue
        else:
            a.append((v.flag(), v.get(as_str=True)))

    return a

def print_status():
    """
    Prints the current progress of the benchmark.
    """

    if os.name == 'nt':
        _ = os.system('cls')
    else:
        _ = os.system('clear')
    
    curr_ds = ds_keys[stats.run['ds']]
    settings_kv = settings['structures'][curr_ds].items()

    print(f"DS: {curr_ds} | Run: {stats.run['curr'] + 1} / {stats.run['total']}")

    output = ''
    for i, kv in enumerate(settings_kv):
        name, value = kv
        output += f'{name} » {value.get(as_str=True)}'
        if not len(settings_kv) == i + 1:
            output += ' | '
    print(output)
    print('=' * len(output))

    for k, v in stats.progress.items():
        completion_rate = v[0] / v[1]
        full_bars = math.floor(completion_rate * 20)
        empty_bars = 20 - full_bars
        
        percentage = math.floor(completion_rate * 100)
        percentage = ' ' * (3 - len(str(percentage))) + str(percentage)
        progress_bar = '|' + '█' * full_bars + ' ' * empty_bars + '|'
        
        print(f'{percentage}%{progress_bar} {v[0]}/{v[1]} ({k})')

    remaining = stats.time['avg'] * (stats.run['total'] - stats.run['curr'])
    m, s = divmod(remaining, 60)
    h, m = divmod(m, 60)
    print('Est. Time Remaining: ' + (f'{int(h):d}:{int(m):02d}:{int(s):02d}' if remaining != 0 else 'Calculating...'))
   
def watch_folder():
    """
    Tracks the current progress of the benchmark in a given folder.
    """

    res_name = settings['general']['result_path'].get()

    while not stop_event.is_set():
        time.sleep(1)

        if pause_event.is_set():
            continue

        for k, v in stats.progress.items():
            _, f_names, _ = next(os.walk(res_name + k))
            f_count = len(f_names)
            if (f_count != v[0]):
                stats.progress[k] = [f_count, v[1]]
                print_status()

def collect_args():
    """
    Accumulates a list of args for every sub-folder / worker based on the global settings.
    """

    combined = build_argvals(settings, [])
    for path in batch_paths:
        args = ['mtsutil', 'dscompare', '-p', path, '--sm', 'sphere']
        
        for flag, value in combined:
            args.append(flag)

            if isinstance(value, Iterable) and not isinstance(value, str):
                for v in value:
                    args.append(v)
            else:
                args.append(value)
        
        commands.append(args)

def create_batches():
    """
    Splits all test envmaps into batches of roughly even size (based on the file size) for more efficient computation.
    """

    print('Creating batches...')

    b_id = 0
    curr_bytes = 0
    curr_files = 0

    base_name = settings['general']['envmap_path'].get()
    res_name = settings['general']['result_path'].get()
    batch_count = settings['testing']['batches']
    folders = os.listdir(os.fsencode(base_name))
    total_bytes = sum(f.stat().st_size for f in Path(base_name).rglob('*') if f.is_file())

    def batch_key() -> str:
        return f'Batch {str(b_id)}'

    def create_folder(p) -> str:
        if not os.path.exists(p):
            os.makedirs(p)
        return p

    for i, path in enumerate(folders):
        folder_name = os.fsdecode(path)
        full_path = os.path.join(base_name, folder_name)
        _, _, files = next(os.walk(full_path))

        if batch_count < 0:
            stats.progress[folder_name] = [0, len(files)]
            batch_paths.append(full_path)
            continue

        for file in files:
            # generate file uid
            file_id = uuid.uuid4()
            _, file_ext = os.path.splitext(file)
            new_name = str(file_id) + file_ext

            # store old folder path
            base_pairing[new_name] = os.path.join(full_path, file)

            # create new results and testing folders
            create_folder(os.path.join(res_name, batch_key()))
            new_path = create_folder(os.path.join(base_name, 'testing', batch_key()))

            # move to new testing folder
            Path(os.path.join(full_path, file)).rename(os.path.join(new_path, new_name))

            curr_bytes += os.path.getsize(os.path.join(new_path, new_name))
            curr_files += 1

            if (curr_bytes >= (total_bytes / batch_count)) and (b_id < batch_count):
                stats.progress[batch_key()] = [0, curr_files]
                batch_paths.append(new_path)
                
                b_id += 1
                curr_bytes = 0
                curr_files = 0

        if (i == len(folders) - 1) and (batch_key() not in stats.progress.keys()):
            stats.progress[batch_key()] = [0, curr_files]
            batch_paths.append(f'{base_name}testing/{batch_key()}/')

def collect_data(curr_settings, wipe=False):
    """
    Collects the data from each csv and stores it into the respective data structure csv.
    If wipe is enabled, the folder is wiped after storing for the next iteration.
    """

    res_path = settings['general']['result_path'].get()
    metrics = settings['testing']['metrics']
    collector = { m: [] for m in metrics }
    ds_index = -1

    # Iterate over all output batches
    for folder_name in stats.progress.keys():
        path = os.path.join(res_path, folder_name)
        for _, dirs, _ in os.walk(path):
            # Iterate over all envmap folders inside a batch
            for out_folder in dirs:
                # Store results from metrics.csv into the collector
                csv_path = os.path.join(path, out_folder, 'metrics.csv')
                with open(csv_path) as file:
                    reader = csv.reader(file, delimiter=',', quotechar='"')

                    data = [row for row in reader if row[1]]
                    indices = [data[0].index(metric) for metric in metrics if metric in data[0]]
                    data = data[-1]
                    ds_index = int(data[0])

                    for metric, index in zip(metrics, indices):
                        base_path = base_pairing.get(f'{out_folder}.exr', base_pairing.get(f'{out_folder}.hdr'))
                        name = Path(base_path).stem
                        collector[metric].append((name, data[index]))
    
    # Wipe folders if requested
    if wipe:
        pause_event.set()

        for folder_name in stats.progress.keys():
            path = os.path.join(res_path, folder_name)
            shutil.rmtree(path)
            os.makedirs(path)

        pause_event.clear()

    # Write to benchmark csv
    ds_name = list(settings['structures'].keys())[ds_index]
    benchmark_path = os.path.join(res_path, 'benchmark', str(int(timestamp)))
    
    for metric in metrics:
        data = collector[metric]

        ds_path = os.path.join(benchmark_path, f'{ds_name}_{metric}.csv')
        payload = {}
        with open(ds_path) as file:
            reader = csv.reader(file, delimiter=',', quotechar='"')
            
            first_row = True
            for row in reader:
                if first_row:
                    payload[None] = row[1:]
                    first_row = False
                else:
                    payload[row[0]] = row[1:]

        if None not in payload:
            payload[None] = ['+'.join(str(s.get()) for s in curr_settings)]
        else:
            payload[None].append('+'.join(str(s.get()) for s in curr_settings))
        
        for envmap, value in data:
            if envmap not in payload:
                payload[envmap] = [value]
            else:
                payload[envmap].append(value)

        with open(ds_path, 'w+') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"')

            for envmap, values in payload.items():
                writer.writerow([envmap, *values])

def start_comparer(command):
    """
    Calls the command to start Mitsuba.
    """

    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def rcall(data, pos=0):
    """
    Recursively goes through all possible setting permutations and creates multiple workers that each start
    Mitsuba with their respective envmap batch for the current settings state.
    """
    
    if pos == len(data):
        now = time.perf_counter()
        collect_args()

        if settings['testing']['multithreading']:
            with ThreadPoolExecutor(max_workers=len(commands)) as executor:
                futures = [executor.submit(start_comparer, c) for c in commands]
                wait(futures)
        else:
            for command in commands:
                start_comparer(command)
        
        commands.clear()
        collect_data(data, wipe=True)

        stats.increment('run')
        stats.increment('time')
        stats.store_time(time.perf_counter() - now)
        return

    while not data[pos].final():
        rcall(data, pos + 1)
        data[pos].next()

    rcall(data, pos + 1)
    data[pos].reset()

def run_test():
    """
    Performs a full benchmark run.
    """

    print('Running test...')
    
    blacklist = settings['general']['blacklist'].get()
    if not all(isinstance(v, int) for v in blacklist):
        warnings.warn('Invalid blacklist content. Falling back to using all data structures.')
        blacklist = []

    sl = settings['general']['samples_learning']

    # Create directories for benchmark results
    res_path = settings['general']['result_path'].get()
    benchmark_path = os.path.join(res_path, 'benchmark', str(int(timestamp)))
    if not os.path.exists(benchmark_path):
        os.makedirs(benchmark_path)

    for ds_name in ds_keys:
        for metric in settings['testing']['metrics']:
            ds_path = os.path.join(benchmark_path, f'{ds_name}_{metric}.csv')
            if not os.path.exists(ds_path):
                with open(ds_path, 'w') as _:
                    pass

    # Start watch thread
    watcher = Thread(target=watch_folder)
    watcher.start()

    # Run benchmark per data structure
    for ds_i in range(len(ds_keys)):
        if ds_i in blacklist:
            continue

        settings['general']['blacklist']._value = list(range(len(ds_keys)))
        settings['general']['blacklist']._value.remove(ds_i)

        fluid_settings = [v for v in settings['structures'][ds_keys[ds_i]].values() if isinstance(v, FluidSetting)]
        stats.run['ds'] = ds_i
        stats.run['total'] = math.prod(fs.steps() for fs in fluid_settings)

        while not sl.final():
            rcall(fluid_settings)
            stats.reset()
            sl.next()

        rcall(fluid_settings)
        stats.reset()
        sl.reset()

def restore_old_folders():
    """
    Restores the batchified folder structure back to its initial state.
    """
    
    if settings['testing']['batches'] < 0:
        return

    base_name = settings['general']['envmap_path'].get()
    test_path = os.path.join(base_name, 'testing')
    for path in os.listdir(os.fsencode(test_path)):
        folder_name = os.fsdecode(path)
        full_path = os.path.join(test_path, folder_name)
        _, _, files = next(os.walk(full_path))

        for file in files:
            old_path = base_pairing.pop(file)
            Path(os.path.join(full_path, file)).rename(old_path)

        if os.path.isfile(full_path):
            os.remove(full_path)

    res_path = settings['general']['result_path'].get()
    for folder_name in stats.progress.keys():
        path = os.path.join(res_path, folder_name)
        shutil.rmtree(path)

def shutdown(signum, frame):
    """
    Gracefully stops script interruptions and restores the initial folder structure.
    """
    
    stop_event.set()
    print('Shutting down...')
    signal.signal(signum, signal.SIG_IGN)
    restore_old_folders()
    sys.exit(0)

signal.signal(signal.SIGINT, shutdown)

if __name__ == '__main__':
    create_batches()
    run_test()
    shutdown(signal.SIGINT, None)