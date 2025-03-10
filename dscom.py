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
# • finished() -> check whether the value has reached its final state (upper bound, !current, etc.) 
##############################################

class FluidSetting:
    def __init__(self, value):
        self.value = value
        self.default = value

    def get(self):
        return self.value

    def reset(self):
        self.value = self.default

class Range(FluidSetting):
    def __init__(self, start, end, func = None):
        super().__init__(start)
        self.end = end
        
        if not func:
            func = lambda x: x + 1
        self.func = func
    
    def next(self):
        self.value = self.func(self.value)

    def finished(self):
        next_value = self.func(self.value)
        return (next_value > self.end)
    
class Toggle(FluidSetting):
    def __init__(self, state):
        super().__init__(state)
    
    def next(self):
        self.value = not self.value

    def finished(self):
        return (self.value != self.default)

##############################################
# Base config.
##############################################

settings = {
    "testing": {
        # Enables multithreading for concurrent evaluation of environment maps.
        "multithreading": True,
        # Set a time limit in seconds after which the application will stop running. If set to '-1', the time limit will be ignored.
        "time_limit": -1,
        # Max. number of envmaps per folder. Set to -1 to disable.
        "batch_size": 100
    },
    "general": {
        # Path to folder containing the envmaps. Set to 'None' to use the default folder. If multithreading is enabled, this setting is ignored.
        "envmap_path": None,
        # Path to folder containing the results. Set to 'None' to use the default folder. If multithreading is enabled, this setting is ignored.
        "result_path": None,
        # Number of learning samples per data structure. These samples are used to create a guiding distribution.
        "samples_learning": Range(start=64, end=65536, func=lambda x: x * 2),
        # Number of guiding samples per data structure. These samples are used to recreate the sampled distribution from the guiding distribution.
        "samples_guiding": 65536,
        # A blacklist specifying which data structures should be skipped in the overall test.
        "blacklist": [],
        # Whether to normalize the envmap each data structure is 'learning' with. (Keep this true unless you know what you're doing.)
        "normalize": True
    },
    "noise": {
        # Whether the initial envmap should be noisified.
        "envmap": False,
        # Whether a variable amount of noise should be introduced to each learning sample.
        "samples": False
    },
    "structures": {
        "sh": {
            "bands": Range(start=1, end=10),
            "depth": Range(start=1, end=20),
            "use_offset": Toggle(True)
        },
        "dt": {
            "frac_loss": "none",
            "dir_filter": "nearest",
            "threshold": Range(start=0.01, end=0.5, func=lambda x: x + 0.01),
            "iterations": -1,
            "max_depth": Range(start=2, end=20)
        },
        "tc": {
            "tilings": Range(start=1, end=8),
            "tiles_x": Range(start=2, end=32),
            "tiles_y": Range(start=2, end=32)
        },
        "btc": {
            "tilings": Range(start=1, end=8),
            "tiles_x": Range(start=1, end=8),
            "tiles_y": Range(start=1, end=8),
            "max_depth": Range(start=2, end=20),
            "subdiv_thresh": Range(start=0.0001, end=0.01, func=lambda x: x * 2)
        },
        "vmf": {
            "components": Range(start=1, end=32),
            "use_ruppert": Toggle(True)
        }
    }
}

base_name = "./data/tests/envmaps/"
res_name = "./data/results/"

##############################################
# Start of code.
##############################################

base_pairing = {}
batch_paths = []
commands = []
progress = {}
futures = None
i = 0

def build_argvals(s, a = []):
    """
    Builds a list of setting values required for cl call

    Parameters
    ----------
    a : list
      List for accumulating all relevant values
    s : obj
      The settings object
    """

    for k, v in s.items():
        if isinstance(v, dict):
            build_argvals(v, a)
        elif k.lower() in ["multithreading", "time_limit", "batch_size", "envmap_path", "result_path"]:
            continue
        else:
            a.append(v.get() if isinstance(v, FluidSetting) else v)

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
    
    while True:
        if futures != None:
            tasks_finished = True
            for future in futures:
                if not future.done():
                    tasks_finished = False
                    break
                    
            if tasks_finished:
                break
        
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

    commands = ["--sl", "--sg", "-b", "-n", "--ne", "--ns", "--shb", "--shd", "--sho", "--dtl", "--dtf", "--dtt", "--dti", "--dtd", "-t", "--tx", "--ty", "--bt", "--btx", "--bty", "--btd", "--btt", "--vc", "--vr"]
    commands = zip(commands, build_argvals(settings))
    
    for path in os.listdir(os.fsencode(base_name)):
        folder_name = os.fsdecode(path)
        full_path = base_name + folder_name
        _, _, files = next(os.walk(full_path))
        progress[folder_name] = [0, len(files)]
        
        args = ["mtsutil", "dscompare", "-p", full_path, "--sm", "sphere"]
        for prefix, value in commands:
            args.append(prefix)
            args.append(value)
        
        commands.append(args)

def create_batches():
    batch_size = settings["testing"]["batch_size"]
    batch_id = 0
    file_count = 0

    for i, path in enumerate(os.listdir(os.fsencode(base_name))):
        folder_name = os.fsdecode(path)
        full_path = base_name + folder_name
        _, _, files = next(os.walk(full_path))
        
        if batch_size < 0:
            progress[folder_name] = [0, len(files)]
            batch_paths.append(full_path)
            continue
        
        for file in files:
            # generate file uid
            file_id = uuid.uuid4()
            _, file_ext = os.path.splitext(file)

            # store old folder path
            new_name = str(file_id) + file_ext
            base_pairing[new_name] = full_path + "/" + file

            # create new folder
            new_path = base_name + "testing/batch_" + str(batch_id) + "/"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
            
            # move to new folder
            Path(full_path + "/" + file).rename(new_path + new_name)
            file_count += 1

            if file_count % batch_size == 0 and file_count != 0:
                progress["Batch " + str(batch_id)] = [0, batch_size]
                batch_paths.append(new_path)
                batch_id += 1
        
        # Unless there was exactly a multiple of 'batch_size' files, the last batch hasn't been tracked yet
        if i == (len(os.listdir(os.fsencode(base_name))) - 1) and ("Batch " + str(batch_id)) not in progress.keys():
            progress["Batch " + str(batch_id)] = [0, file_count % batch_size]
            batch_paths.append(base_name + "testing/batch_" + str(batch_id) + "/")

def restore_old_folders():
    if settings["testing"]["batch_size"] < 0:
        return

    test_path = base_name + "testing/"
    for path in os.listdir(os.fsencode(test_path)):
        folder_name = os.fsdecode(path)
        full_path = test_path + folder_name
        _, _, files = next(os.walk(full_path))

        for file in files:
            old_path = base_pairing.pop(file)
            Path(full_path + "/" + file).rename(old_path)

        if os.path.isfile(full_path):
            os.remove(full_path)
    
def start_comparer(command):
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def all_finished(s, c = True):
    for _, v in s.items():
        if isinstance(v, dict):
            c = all_finished(v, c)
        elif isinstance(v, FluidSetting) and not v.finished():
            return False

    return c

def run_test():
    sl = settings["general"]["samples_learning"]
    sg = settings["general"]["samples_guiding"]

    while not all_finished(settings):
        
        collect_args()

        with ProcessPoolExecutor() as executor:
            executor.submit(watch_folder)
            futures = executor.map(start_comparer, commands)

def sighandler(signum, frame):
    signal.signal(signum, signal.SIG_IGN)
    restore_old_folders()
    sys.exit(0)

signal.signal(signal.SIGINT, sighandler)

if __name__ == '__main__':
    create_batches()
    while True:
        pass
    
    exit(0)
    run_test()

# TODO:
# - Change file number approach to size approach for even distribution
# - Iterate over settings and increase each time
# - Store results in CSV; one folder per ds, one csv for every permutation