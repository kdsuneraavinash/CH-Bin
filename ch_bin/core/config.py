"""
Module with the configuration related utilities and helper functions.

Will also contain the USER_CONFIG variable, which is a
global state variable containing all the currently active user configuration.
To load a new configuration, use `USER_CONFIG.read_file()`.
"""

from configparser import ConfigParser

"""
This represents the configuration of the user.
This will include information such as tool locations as well as
specific parameters set by the user.
"""
USER_CONFIG = ConfigParser()
