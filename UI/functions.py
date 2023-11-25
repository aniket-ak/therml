import subprocess, threading, queue
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd
import json
import h5py
import os

def read_solution(file):
    f = h5py.File(file, "r")
    time_duration = len(f['solution'].keys())
    time_vals = [float(i) for i in list(f['solution'].keys())]
    time_vals.sort()
    keys = [str(i) for i in time_vals]

    first_key = keys[0]
    sol_shape = np.array(list(f['solution'][first_key])).shape

    solution = np.empty((tuple([time_duration] + list(sol_shape))))
    for i,key in enumerate(keys):
        solution[i,:] = np.array(list(f['solution'][key]))
    time_ = list(f['solution'].keys())
    time_ = [float(i) for i in time_]
    time_.sort()
    time_ = np.array(time_)
    f.close()
    return solution, time_

def get_properties(material_data, key_name):
    for material in material_data['materials']:
        if key_name in material:
            return material[key_name]
    return None

def save_to_json(file, inputs):
    (mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z,
    bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z, 
    mold_material, die_material, underfill_material, bumps_materials, substrate_materials, solder_materials,
    mold_surface_material, die_surface_material, underfill_surface_material, bumps_surface_materials, substrate_surface_materials, solder_surface_materials, 
    ambient_temp, start_time, end_time, saveat, top_bc_type, top_bc_value, top_bc_ref_temp, bottom_bc_type, bot_bc_value, bot_bc_ref_temp) = inputs

    settings_file = os.path.join("../", "3d/settings.json")
    with open(settings_file, 'r') as settings_f:
        settings_data = json.load(settings_f)

    material_file = os.path.join("../", "3d/common/materials.json")
    with open(material_file, 'r') as file_:
        materials = json.load(file_)
    
    settings_data["BC"]["Z+"]["type"] = top_bc_type
    settings_data["BC"]["Z+"]["value"]["t_amb"] = top_bc_ref_temp
    settings_data["BC"]["Z+"]["value"]["value"] = top_bc_value
    
    settings_data["BC"]["Z-"]["type"] = bottom_bc_type
    settings_data["BC"]["Z-"]["value"]["t_amb"] = bot_bc_ref_temp
    settings_data["BC"]["Z-"]["value"]["value"] = bot_bc_value

    settings_data["IC"] = ambient_temp

    settings_data["model"]["bodies"]["die"]["size"]["X"] = die_x/1e3
    settings_data["model"]["bodies"]["die"]["size"]["Y"] = die_y/1e3
    settings_data["model"]["bodies"]["die"]["size"]["Z"] = die_z/1e3

    settings_data["model"]["bodies"]["die"]["material"]["name"] = die_material
    settings_data["model"]["bodies"]["die"]["material"]["k"] = get_properties(materials, die_material)["k"]
    settings_data["model"]["bodies"]["die"]["material"]["rho"] = get_properties(materials, die_material)["rho"]
    settings_data["model"]["bodies"]["die"]["material"]["cp"] = get_properties(materials, die_material)["cp"]

    settings_data["model"]["bodies"]["mold"]["size"]["X"] = mold_x/1e3
    settings_data["model"]["bodies"]["mold"]["size"]["Y"] = mold_y/1e3
    settings_data["model"]["bodies"]["mold"]["size"]["Z"] = mold_z/1e3

    settings_data["model"]["bodies"]["mold"]["material"]["name"] = mold_material
    settings_data["model"]["bodies"]["mold"]["material"]["k"] = get_properties(materials, mold_material)["k"]
    settings_data["model"]["bodies"]["mold"]["material"]["rho"] = get_properties(materials, mold_material)["rho"]
    settings_data["model"]["bodies"]["mold"]["material"]["cp"] = get_properties(materials, mold_material)["cp"]

    settings_data["model"]["bodies"]["underfill"]["size"]["X"] = underfill_x/1e3
    settings_data["model"]["bodies"]["underfill"]["size"]["Y"] = underfill_y/1e3
    settings_data["model"]["bodies"]["underfill"]["size"]["Z"] = underfill_z/1e3

    settings_data["model"]["bodies"]["underfill"]["material"]["name"] = underfill_material
    settings_data["model"]["bodies"]["underfill"]["material"]["k"] = get_properties(materials, underfill_material)["k"]
    settings_data["model"]["bodies"]["underfill"]["material"]["rho"] = get_properties(materials, underfill_material)["rho"]
    settings_data["model"]["bodies"]["underfill"]["material"]["cp"] = get_properties(materials, underfill_material)["cp"]

    settings_data["model"]["bodies"]["bumps"]["size"]["X"] = bumps_x/1e3
    settings_data["model"]["bodies"]["bumps"]["size"]["Y"] = bumps_y/1e3
    settings_data["model"]["bodies"]["bumps"]["size"]["Z"] = bumps_z/1e3

    settings_data["model"]["bodies"]["bumps"]["material"]["name"] = bumps_materials
    settings_data["model"]["bodies"]["bumps"]["material"]["k"] = get_properties(materials, bumps_materials)["k"]
    settings_data["model"]["bodies"]["bumps"]["material"]["rho"] = get_properties(materials, bumps_materials)["rho"]
    settings_data["model"]["bodies"]["bumps"]["material"]["cp"] = get_properties(materials, bumps_materials)["cp"]

    settings_data["model"]["bodies"]["solder"]["size"]["X"] = solder_x/1e3
    settings_data["model"]["bodies"]["solder"]["size"]["Y"] = solder_y/1e3
    settings_data["model"]["bodies"]["solder"]["size"]["Z"] = solder_z/1e3

    settings_data["model"]["bodies"]["solder"]["material"]["name"] = solder_materials
    settings_data["model"]["bodies"]["solder"]["material"]["k"] = get_properties(materials, solder_materials)["k"]
    settings_data["model"]["bodies"]["solder"]["material"]["rho"] = get_properties(materials, solder_materials)["rho"]
    settings_data["model"]["bodies"]["solder"]["material"]["cp"] = get_properties(materials, solder_materials)["cp"]

    settings_data["model"]["bodies"]["substrate"]["size"]["X"] = substrate_x/1e3
    settings_data["model"]["bodies"]["substrate"]["size"]["Y"] = substrate_y/1e3
    settings_data["model"]["bodies"]["substrate"]["size"]["Z"] = substrate_z/1e3
    
    settings_data["model"]["bodies"]["substrate"]["material"]["name"] = substrate_materials
    settings_data["model"]["bodies"]["substrate"]["material"]["k"] = get_properties(materials, substrate_materials)["k"]
    settings_data["model"]["bodies"]["substrate"]["material"]["rho"] = get_properties(materials, substrate_materials)["rho"]
    settings_data["model"]["bodies"]["substrate"]["material"]["cp"] = get_properties(materials, substrate_materials)["cp"]

    settings_data["start_time"] = start_time
    settings_data["end_time"] = end_time
    settings_data["dt"] = saveat

    with open(file, 'w') as f:
        json.dump(settings_data,f, indent=2)

class StdoutCapture:
    def __init__(self):
        self.buffer = []

    def write(self, text):
        self.buffer.append(text)

def run_julia(working_dir, scenario_name, run_name, output_queue):
    try:
        command = ["julia", "--project=/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment", 
                             "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment/precompile_.jl", 
                             working_dir, scenario_name, run_name]
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