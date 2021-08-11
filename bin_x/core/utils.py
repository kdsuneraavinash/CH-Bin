import os
from collections import defaultdict

_WARNED: defaultdict = defaultdict(int)


def run_command(command: str) -> None:
    exit_code = os.system(command)
    if exit_code != 0:
        raise Exception(f"Command terminated with non-zero exit code. (code {exit_code})\n" f"\t{command}")


def warn_once(key: str, warning: str):
    if _WARNED[key] == 0:
        print(f"Warning: {warning}")
    _WARNED[key] += 1
