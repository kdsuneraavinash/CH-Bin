import os


def run_command(command: str) -> None:
    exit_code = os.system(command)
    if exit_code != 0:
        raise Exception(f"Command terminated with non-zero exit code. (code {exit_code})\n" f"\t{command}")
