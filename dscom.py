import subprocess
import os
import time
import math
from concurrent.futures import ProcessPoolExecutor

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
        #
        "result_path": None,
        #
        "samples_learning": [64, 65536],
        #
        "samples_guiding": [64, 65536],
        #
        "blacklist": [],
        #
        "normalize": True,
    },
    "noise": {
        #
        "envmap": False,
        #
        "samples": False
    },
    "structures": {
        "sh": {
            "bands": [1, 10, 1],
            "depth": [1, 20, 1],
            "use_offset": True
        },
        "dt": {
            "frac_loss": "none",
            "dir_filter": "nearest",
            "threshold": [0.01, 0.5, 0.01],
            "iterations": -1,
            "max_depth": [2, 20, 1]
        },
        "tc": {
            "tilings": [1, 8, 1],
            "tiles_x": [2, 32, 1],
            "tiles_y": [2, 32, 1]
        },
        "btc": {
            "tilings": [1, 8, 1],
            "tiles_x": [1, 8, 1],
            "tiles_y": [1, 8, 1],
            "max_depth": [2, 20, 1],
            "subdiv_thresh": [0.0001, 0.01, 0.0001]
        },
        "vmf": {
            "components": [1, 32, 1],
            "use_ruppert": True
        }
    }
}

base_name = "./data/tests/envmaps/"
res_name = "./data/results/"

commands = []
progress = {}
futures = None
i = 0

def print_status():
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
    for path in os.listdir(os.fsencode(base_name)):
        folder_name = os.fsdecode(path)
        full_path = base_name + folder_name
        _, _, files = next(os.walk(full_path))
        progress[folder_name] = [0, len(files)]
        
        args = ["mtsutil", "dscompare", "-p", full_path, "--sl", "65536", "--sm", "sphere"]
        commands.append(args)
    
def start_comparer(command):
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
if __name__ == '__main__':
    collect_args()
    
    with ProcessPoolExecutor(max_workers=7) as executor:
        executor.submit(watch_folder)
        futures = executor.map(start_comparer, commands)