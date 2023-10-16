import subprocess, threading, queue

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