import argparse
import logging
import subprocess
import time
import numpy as np
from pathlib import Path


def measure_runtime(command):
    start_time = time.time()
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command '{command}' failed with error: {e}")
        return
    end_time = time.time()

    runtime = end_time - start_time
    logging.info(f"The command '{command}' took {runtime} seconds to run.")


def runtime_benchmark(path_to_folder: Path) -> np.ndarray:
    # List of commands to run
    commands = [
        "docker run --rm -v $(pwd):/data nanodisco",  # Nanodisco
        "docker run --rm -v $(pwd):/data meme-chip",  # MEME-CHIP
        "perl /path/to/cisfinder.pl",  # CisFinder
        "python2 /path/to/vconv.py",  # vConv
    ]

    # Initialize an empty list to store the runtimes
    runtimes = []

    # Iterate over each file in the folder
    for file in Path(path_to_folder).iterdir():
        if file.is_file():
            # Iterate over each command
            for command in commands:
                start_time = time.time()
                try:
                    # Run the command with the file as input
                    subprocess.run(f"{command} {file}", shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Command '{command}' failed with error: {e}")
                    continue
                end_time = time.time()

                # Calculate the runtime and add it to the list
                runtime = end_time - start_time
                runtimes.append(runtime)

    # Convert the list of runtimes to a numpy array and return it
    return np.array(runtimes)


def main():
    parser = argparse.ArgumentParser(
        description="This script measures the runtime of the specified tools. "
        "Each command should be enclosed in quotes. "
        'Example: python3 script.py "docker run --rm -v $(pwd):/data nanodisco" '
        '"docker run --rm -v $(pwd):/data meme-chip" "perl /path/to/cisfinder.pl" '
        '"python2 /path/to/vconv.py"'
    )
    parser.add_argument(
        "commands",
        nargs="+",
        help="The commands to measure the runtime of. Each command should be enclosed in quotes.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    for command in args.commands:
        measure_runtime(command)


if __name__ == "__main__":
    main()
