import subprocess, threading, queue
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd
import json
import h5py

def read_solution(file):
    f = h5py.File(file, "r")
    time_duration = len(f['solution'].keys())
    sol_shape = np.array(list(f['solution']['1'])).T.shape

    solution = np.empty((tuple([time_duration] + list(sol_shape))))
    solution.shape
    for i in range(1,time_duration):
        solution[i,:] = np.array(list(f['solution'][str(i)])).T
    time_ = list(f['solution'].keys())
    time_ = [int(i) for i in time_]
    time_.sort()
    time_ = np.array(time_)
    return solution, time_

def save_to_json(file, inputs):
    (mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z,
    bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z, 
    mold_material, die_material, epoxy_material, bumps_materials, substrate_materials, solder_materials,
    mold_surface_material, die_surface_material, epoxy_surface_material, bumps_surface_materials, substrate_surface_materials, solder_surface_materials, 
    ambient_temp, start_time, end_time, top_bc_type, top_bc_value, top_bc_ref_temp, bottom_bc_type, bot_bc_value, bot_bc_ref_temp) = inputs
    print(file)
    with open(file, 'r') as f:
        settings_data = json.load(f)

    print(settings_data)

class StdoutCapture:
    def __init__(self):
        self.buffer = []

    def write(self, text):
        self.buffer.append(text)

def run_julia(working_dir, scenario_name, run_name, output_queue):
    try:
        command = ["julia", "--project=/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment", 
                             "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment/precompile_.jl", 
                             "-t", "4", "-working_dir", working_dir, "-power", scenario_name, "-run_name", 
                             run_name]
        stdout_capture = StdoutCapture()

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )

        def capture_output(stream, output_type):
            for line in stream:
                #print(f"{output_type}: {line.strip()}")
                stdout_capture.write(f"{output_type}: {line.strip()}\n")

        stdout_thread = threading.Thread(target=capture_output, args=(process.stdout, f"__{scenario_name}__ "))
        stderr_thread = threading.Thread(target=capture_output, args=(process.stderr, "STDERR"))

        stdout_thread.start()
        stderr_thread.start()

        process.wait()

        stdout_thread.join()
        stderr_thread.join()

        if output_queue is not None:
            output_queue.put(''.join(stdout_capture.buffer))

    except Exception as e:
        print(f"Error: {str(e)}")