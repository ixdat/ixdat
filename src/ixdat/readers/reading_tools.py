"""Module with possibly general-use tools for readers"""

import urllib.request
from ..config import CFG


def url_to_file(url, file_name="temp", directory=None):
    """Copy the contents of the url to a file and return its Path."""
    directory = directory or CFG.ixdat_temp_dir
    suffix = "." + str(url).split(".")[-1]
    path_to_file = (directory / file_name).with_suffix(suffix)
    urllib.request.urlretrieve(url, path_to_file)
    return path_to_file
