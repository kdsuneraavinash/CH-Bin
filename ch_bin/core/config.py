"""
Module with the configuration related utilities and helper functions.

Will also contain the USER_CONFIG variable, which is a
global state variable containing all the currently active user configuration.
To load a new configuration, use `USER_CONFIG.read_file()`.
"""

import logging
from configparser import ConfigParser

"""
This represents the configuration of the user.
This will include information such as tool locations as well as
specific parameters set by the user.
"""
USER_CONFIG = ConfigParser()


def initialize_logger():
    """
    Set the logging functionality default.
    Will output to a file name ch-bin.log in the current directory.
    """
    logging.basicConfig(
        filename="ch-bin.log",
        filemode="w",
        format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        level=logging.DEBUG,
    )
