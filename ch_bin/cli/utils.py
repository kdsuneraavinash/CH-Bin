import sys
import traceback

import click


def handle_error(e: Exception):
    click.secho(f"Error: {e}", fg="red", bold=True)
    with open("errors.log", "a") as fa:
        fa.write(f"{e}\n{traceback.format_exc()}\n")
    sys.exit(1)
