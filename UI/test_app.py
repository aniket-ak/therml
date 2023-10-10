import subprocess
working_dir = "/Users/aniket/Documents/MarlinSim/04_testing/scenarios/"
power_file = "file_1.csv"
simulation_name = "name__1"
result = subprocess.run(["julia", "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/3d_fvm.jl", "-t", "4", "-working_dir", working_dir, "-power", power_file, "-run_name", simulation_name], capture_output=True, text=True)
    
# Check for errors
if result.returncode != 0:
    print("Error:", result.stderr)

print(result.stdout.strip())