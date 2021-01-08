
"""Initial setup.py

TODO: This file is very rudimentary and is setup purely to enable tox to
run. None of the information in it reflects necessarily what it should be and
will be updated to proper values by project owner.

"""

import os
import re
import codecs
import setuptools

META_PATH = os.path.join("src", "ixdat", "__init__.py")
HERE = os.path.abspath(os.path.dirname(__file__))
PACKAGES = setuptools.find_packages(where="./src")


def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE.
    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
        META_FILE, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError(
        "Unable to find __{meta}__ string.".format(meta=meta)
    )


META_FILE = read(META_PATH)

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
    python_requires='>=3.6',
)
