##############################################
# Automated DS::compare testing script for all implemented data structures.
# Stores the results as CSV file per data structure, per configuration, in a folder.
# Lots of folders, in fact.
#
# How to add your own structure to this script:
# TBA!
#
# Feel free to do whatever you want to this script.
# Michael Eickmeyer, 2025 @ TU Wien.
##############################################

import subprocess
import os
import time
import math
import signal
import sys
import uuid
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

##############################################
# A bunch of classes handling fluid setting parameters.
#
# Fluid settings allow for a parameter to increase, decrease or toggle
# within this test script. For instance, the Range class takes an optional
# 'func' parameter which allows the user to specify how exactly the value
# should increase within the range.
#
# The Toggle class allows the user to easily toggle between two truthy/falsy values.
#
# All fluid settings consist of the functions...
# • get() -> obtain the current value
# • reset() -> set the value back to its initial state
# • next() -> increase the value to its next state
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
            return ' '.join([str(v) for v in self._value])

        return str(self._value)
    
    def flag(self):
        return self._flag

class FluidSetting(Value):
    def __init__(self, flag, value):
        super().__init__(flag, value)
        self._default = value

    def reset(self):
        self._value = self._default

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
    
    def next(self):
        self._value = self._func(self._value)

    def final(self):
        next_value = self._func(self._value)
        return (next_value > self._end)
    
class Toggle(FluidSetting):
    def __init__(self, flag, state):
        super().__init__(flag, state)
    
    def next(self):
        self._value = not self._value

    def final(self):
        return (self._value != self._default)

##############################################
# Base config.
##############################################

settings = {
    "testing": {
        # Enables multithreading for concurrent evaluation of environment maps.
        "multithreading": True,
        # Set a time limit in seconds after which the application will stop running. If set to '-1', the time limit will be ignored.
        "time_limit": -1,
        # (Max.) Number of batches the envmaps get divided into. Set to -1 to disable & use the provided folder structure.
        "batches": 12
    },
    "general": {
        # Path to folder containing the envmaps.
        "envmap_path": Value('p', './data/tests/envmaps/'),
        # Path to folder containing the results.
        "result_path": Value('rp', './data/results/'),
        # Number of learning samples per data structure. These samples are used to create a guiding distribution.
        "samples_learning": Range('sl', start=64, end=65536, func=lambda x: x * 2),
        # Number of guiding samples per data structure. These samples are used to recreate the sampled distribution from the guiding distribution.
        "samples_guiding": Value('sg', 65536),
        # A blacklist specifying which data structures should be skipped in the overall test.
        "blacklist": Value('b', [0]),
        # Whether to normalize the envmap each data structure is 'learning' with. (Keep this true unless you know what you're doing.)
        "normalize": Value('n', True)
    },
    "noise": {
        # Whether the initial envmap should be noisified.
        "envmap": Value('ne', False),
        # Whether a variable amount of noise should be introduced to each learning sample.
        "samples": Value('ns', False)
    },
    "structures": {
        "unidir": { },
        "sh": {
            "bands": Range('shb', start=1, end=10),
            "depth": Range('shd', start=1, end=20),
            "use_offset": Toggle('sho', True)
        },
        "dt": {
            "frac_loss": Value('dtl', 'none'),
            "dir_filter": Value('dtf', 'nearest'),
            "threshold": Range('dtt', start=0.01, end=0.5, func=lambda x: x + 0.01),
            "iterations": Value('dti', -1),
            "max_depth": Range('dtd', start=2, end=20)
        },
        "tc": {
            "tilings": Range('t', start=1, end=8),
            "tiles_x": Range('tx', start=2, end=32),
            "tiles_y": Range('ty', start=2, end=32)
        },
        "btc": {
            "tilings": Range('bt', start=1, end=8),
            "tiles_x": Range('btx', start=1, end=8),
            "tiles_y": Range('bty', start=1, end=8),
            "max_depth": Range('btd', start=2, end=20),
            "subdiv_thresh": Range('btt', start=0.0001, end=0.01, func=lambda x: x * 2)
        },
        "vmf": {
            "components": Range('vc', start=1, end=32),
            "use_ruppert": Toggle('vr', True)
        }
    }
}

##############################################
# Start of code.
##############################################

base_pairing = {}
batch_paths = []
commands = []
progress = {}

class FHolder:
    futures = None
fholder = FHolder()

def build_argvals(s, a = []):
    """
    Builds a list of setting values required for cl call

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
        elif k.lower() in ["multithreading", "time_limit", "batches", "envmap_path", "result_path"]:
            continue
        else:
            a.append((v.flag(), v.get(as_str=True)))

    return a


def print_status():
    """
    Prints the current progress of the benchmark
    """

    if os.name == 'nt':
        _ = os.system('cls')
    else:
        _ = os.system('clear')
        
    for k, v in progress.items():
        completion_rate = v[0] / v[1]
        full_bars = math.floor(completion_rate * 20)
        empty_bars = 20 - full_bars
        
        percentage = math.floor(completion_rate * 100)
        percentage = ' ' * (3 - len(str(percentage))) + str(percentage)
        progress_bar = '|' + '█' * full_bars + ' ' * empty_bars + '|'
        
        print(f'{percentage}%{progress_bar} {v[0]}/{v[1]} ({k})')
   
def watch_folder():
    """
    Tracks the current progress of the benchmark in a given folder
    """

    res_name = settings['general']['result_path'].get()
    
    while True:
        if fholder.futures != None:
            tasks_finished = True
            for future in fholder.futures:
                if not future.done():
                    tasks_finished = False
                    break
                    
            if tasks_finished:
                break
        
        # TODO: Fix
        # TODO: Make sure old res folders are used, even if batched up!
        # TODO: Fix description of f-settings
        for k, v in progress.items():
            _, f_count, _ = next(os.walk(res_name + k))
            f_count = len(f_count)
            if (f_count != v[0]):
                progress[k] = [f_count, v[1]]
                print_status()

def collect_args():
    """
    Accumulates a list of args for every sub-folder / worker based on the global settings
    """

    combined = build_argvals(settings)
    for path in batch_paths:
        args = ["mtsutil", "dscompare", "-p", path, "--sm", "sphere"]
        
        for flag, value in combined:
            args.append(flag)
            args.append(value)
        
        commands.append(args)

def create_batches():
    b_id = 0
    curr_bytes = 0
    curr_files = 0

    base_name = settings['general']['envmap_path'].get()
    batch_count = settings['testing']['batches']
    folders = os.listdir(os.fsencode(base_name))
    total_bytes = sum(os.path.getsize(base_name + os.fsdecode(f)) for f in folders)

    def batch_key() -> str:
        return f'Batch {str(b_id)}'

    for i, path in enumerate(folders):
        folder_name = os.fsdecode(path)
        full_path = os.path.join(base_name, folder_name)
        _, _, files = next(os.walk(full_path))
        
        if batch_count < 0:
            progress[folder_name] = [0, len(files)]
            batch_paths.append(full_path)
            continue
        
        for file in files:
            # generate file uid
            file_id = uuid.uuid4()
            _, file_ext = os.path.splitext(file)
            new_name = str(file_id) + file_ext

            # store old folder path
            base_pairing[new_name] = os.path.join(full_path, file)

            # create new folder
            new_path = os.path.join(base_name, "testing", batch_key())
            if not os.path.exists(new_path):
                os.makedirs(new_path)

            # move to new folder
            Path(os.path.join(full_path, file)).rename(os.path.join(new_path, new_name))

            curr_bytes += os.path.getsize(os.path.join(new_path, new_name))
            print(curr_bytes)
            curr_files += 1

            if (curr_bytes >= (total_bytes / batch_count)) and (b_id < batch_count):
                progress[batch_key()] = [0, curr_files]
                batch_paths.append(new_path)
                b_id += 1
                curr_files = 0
        
        if (i == len(folders) - 1) and (batch_key() not in progress.keys()):
            progress[batch_key()] = [0, curr_files]
            batch_paths.append(base_name + "testing/" + batch_key() + "/")
            
def start_comparer(command):
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def rcall(data, pos=0):
    if pos == len(data):
        commands.clear()
        collect_args()

        with ProcessPoolExecutor() as executor:
            executor.submit(watch_folder)
            fholder.futures = executor.map(start_comparer, commands)

        return

    while not data[pos].final():
        rcall(data, pos + 1)
        data[pos].next()

    rcall(data, pos + 1)
    data[pos].reset()

def run_test():
    blacklist = settings['general']['blacklist'].get()
    if not all(isinstance(v, int) for v in blacklist):
        warnings.warn('Invalid blacklist content. Falling back to using all data structures.')
        blacklist = []

    sl = settings['general']['samples_learning']
    ds_keys = list(settings["structures"].keys())

    for ds_i in range(len(ds_keys)):
        if ds_i in blacklist:
            continue

        settings['general']['blacklist'].value = list(range(len(ds_keys)))
        settings['general']['blacklist'].value.remove(ds_i)

        fluid_settings = [v for v in settings["structures"][ds_keys[ds_i]].values() if isinstance(v, FluidSetting)]
        while not sl.final():
            rcall(fluid_settings)
            sl.next()

        rcall(fluid_settings)
        sl.reset()

def restore_old_folders():
    if settings['testing']['batches'] < 0:
        return

    base_name = settings['general']['envmap_path']
    test_path = os.path.join(base_name.get(), 'testing')
    for path in os.listdir(os.fsencode(test_path)):
        folder_name = os.fsdecode(path)
        full_path = os.path.join(test_path, folder_name)
        _, _, files = next(os.walk(full_path))

        for file in files:
            old_path = base_pairing.pop(file)
            Path(os.path.join(full_path, file)).rename(old_path)

        if os.path.isfile(full_path):
            os.remove(full_path)

def sighandler(signum, frame):
    signal.signal(signum, signal.SIG_IGN)
    restore_old_folders()
    sys.exit(0)

signal.signal(signal.SIGINT, sighandler)

if __name__ == '__main__':
    create_batches()
    run_test()
    restore_old_folders()

# TODO:
# - Store results in CSV; one folder per ds, one csv for every permutation
# - Add time limit
