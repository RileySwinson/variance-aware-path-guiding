import subprocess
import os
import time
import math
from concurrent.futures import ProcessPoolExecutor

# TODO:
# - Split up folders so they're 100 images at most for faster processing
# - Iterate over settings and increase each time
# - Store results in CSV; one folder per ds, one csv for every permutation

class FluidSetting:
    pass

class Range(FluidSetting):
    def __init__(self, start, end, func = None):
        self.bounds = (start, end)
        self.value = start
        
        if not func:
            func = lambda x: x + 1
        self.func = func

    def get(self):
        return self.value
    
    def next(self):
        self.value = self.func(self.value)

    def finished(self):
        next_value = self.func(self.value)
        return (next_value > self.bounds[-1])
    
    def reset(self):
        self.value = self.bounds[0]
    
class Toggleable(FluidSetting):
    def __init__(self, state):
        self.state = state
        self.initial = state

    def get(self):
        return self.state
    
    def next(self):
        self.state = not self.state

    def finished(self):
        return (self.state != self.initial)
    
    def reset(self):
        self.state = self.initial

###################################################

settings = {
    "testing": {
        # Enables multithreading for concurrent evaluation of environment maps.
        "multithreading": True,
        # Set a time limit in seconds after which the application will stop running. If set to '-1', the time limit will be ignored.
        "time_limit": -1
    },
    "general": {
        # Path to folder containing the envmaps. Set to 'None' to use the default folder. If multithreading is enabled, this setting is ignored.
        "envmap_path": None,
        # Path to folder containing the results. Set to 'None' to use the default folder. If multithreading is enabled, this setting is ignored.
        "result_path": None,
        #
        "samples_learning": Range(start=64, end=65536, func=lambda x: x * 2),
        #
        "samples_guiding": 65536,
        #
        "blacklist": [],
        #
        "normalize": True
    },
    "noise": {
        #
        "envmap": False,
        #
        "samples": False
    },
    "structures": {
        "sh": {
            "bands": Range(start=1, end=10),
            "depth": Range(start=1, end=20),
            "use_offset": Toggleable(True)
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
            "use_ruppert": Toggleable(True)
        }
    }
}

base_name = "./data/tests/envmaps/"
res_name = "./data/results/"

###################################################
################# END CONFIG ######################
###################################################

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
        elif k.lower() in ["multithreading", "time_limit", "envmap_path", "result_path"]:
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

    # dont forget to add sm at end!
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
    
def start_comparer(command):
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def all_finished(s, c = True):
    for _, v in s.items():
        if isinstance(v, dict):
            c = all_finished(v, c)
        elif isinstance(v, FluidSetting) and not v.finished():
            return False

    return c

def run():
    sl = settings["general"]["samples_learning"]
    sg = settings["general"]["samples_guiding"]

    while not all_finished(settings):
        
        collect_args()

        with ProcessPoolExecutor() as executor:
            executor.submit(watch_folder)
            futures = executor.map(start_comparer, commands)

if __name__ == '__main__':
    run()