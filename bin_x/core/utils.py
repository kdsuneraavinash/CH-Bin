import os
import subprocess
from pathlib import Path
from typing import List, Optional, Union


def run_command(*command: Union[str, Path], env_paths: Optional[List[str]] = None) -> None:
    """

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
