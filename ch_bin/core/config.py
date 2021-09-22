"""
Module with the configuration related utilities and helper functions.

Will also contain the USER_CONFIG variable, which is a
global state variable containing all the currently active user configuration.
To load a new configuration, use `USER_CONFIG.read_file()`.
"""

import logging
from configparser import ConfigParser

import click

"""
This represents the configuration of the user.
This will include information such as tool locations as well as
specific parameters set by the user.
"""
USER_CONFIG = ConfigParser()


class CHBinLogHandler(logging.FileHandler):
    """
    Custom log handler which logs to a file but also outputs
    INFO/ERROR messages to the terminal using click.
    """

    def __init__(self, filename, mode="w", encoding="UTF-8", delay=False):
        super().__init__(filename, mode, encoding, delay)

    def emit(self, record: logging.LogRecord):
        if record.levelno == logging.ERROR:
            click.secho(f"Error: {record.getMessage()}", fg="red", bold=True)
        elif record.levelno == logging.INFO:
            message = record.getMessage()
            if message.startswith(">>"):
                click.secho(message, bold=True, fg="green")
            else:
                click.secho(message)
        return super().emit(record)


def initialize_logger():
    """
    Set the logging functionality default.
    Will output to a file name ch-bin.log in the current directory.
    """
    logging.basicConfig(
        handlers=[CHBinLogHandler("ch-bin.log")],
        format="%(asctime)s.%(msecs)03d %(name)-40s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )
    logging.getLogger("numba").setLevel(logging.WARNING)
