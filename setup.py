"""This module is used for building the package for distribution"""

import os
import re
import setuptools

META_PATH = os.path.join("src", "ixdat", "__init__.py")
HERE = os.path.abspath(os.path.dirname(__file__))
PACKAGES = setuptools.find_packages(where="./src")


# The `read` and `find_meta` functions are shamelessly stolen from
# https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
# and their doc strings adapted


def read(*parts):
    """Build an absolute path from *parts* and return the contents of
    the resulting file as a string.

    Assume UTF-8 encoding.

    """
    path_to_file = os.path.join(HERE, *parts)
    with open(path_to_file, "r") as f:
        return f.read()


META_FILE = read(META_PATH)


def find_meta(meta):
    """Extract a piece of double underscore defined metadata from the
    global META_FILE (defined above).

    If e.g. in META_FILE there is the code::

     __author__ = "Guido"

    then find_meta("author") will return "Guido"

    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta), META_FILE, re.M
    )
    if meta_match:
        print(f"found {meta}: '{meta_match.group(1)}'")  # debugging
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


setuptools.setup(
    name=find_meta("title"),
    version=find_meta("version"),
    license=find_meta("license"),
    author=find_meta("author"),
    author_email=find_meta("email"),
    description=find_meta("description"),
    long_description=read("README.rst"),
    long_description_content_type="text/x-rst",
    url="https://github.com/ixdat/ixdat",
    packages=PACKAGES,
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=read("requirements.txt").split("\n"),
    python_requires=">=3.6",
)
