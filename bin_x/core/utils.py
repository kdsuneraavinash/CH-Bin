import os
import pstats
import subprocess
from cProfile import Profile
from functools import wraps
from pathlib import Path
from typing import List, Optional, Union

import pyprof2calltree


def run_command(*command: Union[str, Path], env_paths: Optional[List[str]] = None) -> None:
    """
    Credits: https://github.com/hannosch/profilestats

    :param command: Command to run. Provide each argument one by one.
    :param env_paths: PATHs to append to the system PATH variable.
    """
    environment = os.environ.copy()
    if env_paths is not None:
        env_paths.append(environment["PATH"])
        environment["PATH"] = ":".join(env_paths)
    process = subprocess.Popen(command, env=environment)
    process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command terminated with non-zero exit code. (code {process.returncode})\n" f"\t{command}")


def profile(
    cumulative=True,
    print_stats=0,
    sort_stats="cumulative",
    dump_stats=False,
    profile_filename="profilestats.out",
    callgrind_filename="callgrind.out",
):
    def closure(func):
        @wraps(func)
        def decorator(*args, **kwargs):
            result = None
            if cumulative:
                global profiler
            else:
                profiler = Profile()
            try:
                result = profiler.runcall(func, *args, **kwargs)
            finally:
                if dump_stats:
                    profiler.dump_stats(profile_filename)
                stats = pstats.Stats(profiler)
                conv = pyprof2calltree.CalltreeConverter(stats)
                with open(callgrind_filename, "w") as fd:
                    conv.output(fd)
                if print_stats:
                    stats.strip_dirs().sort_stats(sort_stats).print_stats(print_stats)
            return result

        return decorator

    return closure


profiler = Profile()
