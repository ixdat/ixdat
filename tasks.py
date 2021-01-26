"""Definition of invoke tasks"""

import sys
import configparser
import platform
from shutil import rmtree
from invoke import task
from pathlib import Path
from subprocess import check_call, CalledProcessError, check_output, DEVNULL


THIS_DIR = Path(__file__).parent
SOURCE_DIR = THIS_DIR / "src" / "ixdat"
# Patterns to match for files of directories that should be deleted in the clean task
CLEAN_PATTERNS = ("__pycache__", "*.pyc", "*.pyo", ".mypy_cache")


tox_config = configparser.ConfigParser()
tox_config.read(THIS_DIR / "tox.ini")


# ### QA tasks


@task(aliases=["lint"])
def flake8(context):
    """Run the flake8 task

    The ``context`` argument is automatically passed in by invoke and
    represents the context the commands is to be "invoked" in. See
    http://docs.pyinvoke.org/en/stable/api/context.html for details.

    """
    print("# flake8")
    return context.run("flake8").return_code


@task(aliases=["test", "tests"])
def pytest(context):
    """Run the pytest task

    See docstring of :func:`flake8` for explanation of `context` argument

    """
    print("# pytest")
    return context.run("pytest").return_code


@task(aliases=["QA", "qa", "check"])
def checks(context):
    """Run all QA checks

    See docstring of :func:`flake8` for explanation of `context` argument

    """
    combined_return_code = flake8(context)
    combined_return_code += pytest(context)
    if combined_return_code == 0:
        print()
        print(r"+----------+")
        print(r"| All good |")
        print(r"+----------+")


@task
def tox(context, single=False):
    """Run tox for the python interpreters available on your system

    See docstring of :func:`flake8` for explanation of `context` argument

    Args:

        single (bool): Whether or not to restrict the tox run to the
            first python environment found, instead of all of
            the. This option is available at the command line as the
            ``single`` or ``-s`` option and used like so:
            `invoke tox --single`

    """
    environments = tox_config["tox"]["envlist"].split(", ")

    # Check which pythons are available on the system
    if platform.system() == "Windows":
        environments_to_run = filter_tox_environments_windows(
            environments, single=single
        )
    else:
        environments_to_run = filter_tox_environments_linux(environments, single=single)

    context.run("tox -p auto -e " + ",".join(environments_to_run))


def filter_tox_environments_linux(environments, single=False):
    """Filter tox environments to only those available on this Linux system

    Args:
        single (bool): Whether or not to return just a single **python** environment,
            no-matter how many are available

    """
    environments_to_run = []
    found_at_least_one_python = False
    for environment in environments:
        if environment.startswith("py"):
            if environment.startswith("pypy"):
                command = environment
            else:
                # The environments look like: py36
                command = "python{}.{}".format(*environment[2:4])

            # Check that the executable exists
            try:
                check_call([command, "--version"], stdout=DEVNULL)
            except (CalledProcessError, FileNotFoundError):
                continue

            # Certain version of Ubuntu may have a old "reduced"
            # Python version, with an "m" suffix. It is insufficient
            # for tox so check that it isn't there.
            try:
                check_call([command + "m", "--version"], stdout=DEVNULL)
                continue
            except (CalledProcessError, FileNotFoundError):
                pass

            if found_at_least_one_python and single:
                continue

            environments_to_run.append(environment)
            found_at_least_one_python = True

        else:
            environments_to_run.append(environment)
    return environments_to_run


def filter_tox_environments_windows(environments, single=False):
    """Filter tox environments to only those available on this Windows system


    Args:
        single (bool): Whether or not to return just a single **python** environment,
            no-matter how many are available

    """
    python_versions = []
    for path in sys.path:
        python_path = Path(path) / "python.exe"
        if python_path.is_file():
            version_string = (
                check_output([str(python_path), "--version"]).decode("ascii").strip()
            )
            # The version string comes back as: Python 3.7.2
            # hence the parsing below
            version_numbers = version_string.replace("Python ", "").split(".")
            version = "".join(version_numbers[:2])
            python_versions.append("py" + version)

    environments_to_run = []
    found_at_least_one_python = False
    for environment in environments:
        if environment.startswith("py"):
            if environment in python_versions:
                if found_at_least_one_python and single:
                    continue

                environments_to_run.append(environment)
                found_at_least_one_python = True
        else:
            environments_to_run.append(environment)

    return environments_to_run


# ### Maintenance tasks


@task
def clean(context, dryrun=False):
    """Clean the repository of temporary files and caches

    See docstring of :func:`flake8` for explanation of `context` argument

    Arguments:

        dryrun (bool): Whether to perform a dryrun, meaning only show
            what would be deleted, but don't actually delete it. This
            option can accessed on the command line as the
            ``--dryrun`` or ``-d`` option like so:
            `invoke clean --dryrun`

    """
    if dryrun:
        print("CLEANING DRYRUN")
    for clean_pattern in CLEAN_PATTERNS:
        for cleanpath in THIS_DIR.glob("**/" + clean_pattern):
            if cleanpath.is_dir():
                print("DELETE DIR :", cleanpath)
                if not dryrun:
                    rmtree(cleanpath)
            else:
                print("DELETE FILE:", cleanpath)
                if not dryrun:
                    cleanpath.unlink()
