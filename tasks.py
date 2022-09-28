"""Definition of invoke tasks

All functions in this file that is decorated with the `@task` decorator,
constitutes an invoke task, which can be execute from the command line by
calling invoke with the name of the task, which is the same as the name of the
function (or using one of the aliasses possible listed in the @task decorator).
E.g. to install all normal dependencies and development dependencies run::

    invoke dependencies

or::

    invoke deps

To see a full list of the available tasks and a brief summary of what they do,
call::

    invoke --list

To get the full description of a task, call invoke --help on it::

    invoke --help deps

Read more about invoke here: https://www.pyinvoke.org/

"""

import os
import sys
import configparser
import platform
from shutil import rmtree
from invoke import task
from pathlib import Path
from subprocess import check_call, CalledProcessError, check_output, DEVNULL


THIS_DIR = Path(__file__).parent
SOURCE_DIR = THIS_DIR / "src" / "ixdat"
TESTS_DIR = THIS_DIR / "tests"
# NOTE The development_scripts folder below is only used in the
# actions for black formatting, but not linting etc.
DEV_SCRIPTS_DIR = THIS_DIR / "development_scripts"
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
    with context.cd(THIS_DIR):
        return context.run(f"flake8 {SOURCE_DIR} {TESTS_DIR}").return_code


@task(
    aliases=["test", "tests"],
    help={
        "color": "Whether to display pytest output in color, 'yes' or 'no'",
        "external": "Also run external tests from submodules, disabled by default",
    },
)
def pytest(context, color="yes", external=False):
    """Run the pytest task

    See docstring of :func:`flake8` for explanation of `context` argument

    """
    print("# pytest")
    if platform.system() == "Windows":
        color = "no"
    args = []
    if external:
        args.append("--external")
    with context.cd(THIS_DIR):
        return context.run(f"pytest tests --color {color} {' '.join(args)}").return_code


@task(
    aliases=(
        "check_black",
        "black_check",
        "bc",
    )
)
def check_code_format(context):
    """Check that the code, tests and development_scripts are black formatted

    See docstring of :func:`flake8` for explanation of `context` argument

    """
    print("### Checking code style ...")
    with context.cd(THIS_DIR):
        result = context.run(f"black --check {SOURCE_DIR} {TESTS_DIR} {DEV_SCRIPTS_DIR}")
    return result.return_code


@task(aliases=["QA", "qa", "check"])
def checks(context):
    """Run all QA checks

    See docstring of :func:`flake8` for explanation of `context` argument

    """
    combined_return_code = flake8(context)
    combined_return_code += pytest(context)
    combined_return_code += check_code_format(context)
    if combined_return_code == 0:
        print()
        print(r"+----------+")
        print(r"| All good |")
        print(r"+----------+")


@task(aliases=("black",))
def format_code(context):
    """Format all spitze and tools code with black"""
    context.run(f"black {SOURCE_DIR} {TESTS_DIR} {DEV_SCRIPTS_DIR}")


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
    with context.cd(THIS_DIR):
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


@task(aliases=["deps"])
def dependencies(context):
    """Install all dedencendies required for development

    This corresponds to:
     * Upgrade pip
     * Install/upgrade all normal dependencies
     * Install/upgrade all development dependencies

    See docstring of :func:`flake8` for explanation of `context` argument
    """
    # See https://stackoverflow.com/a/1883251/11640721 for virtual env detection trick
    conda_environment = os.environ.get("CONDA_PREFIX")
    # conda_environment is None if not using Anaconda python
    if conda_environment is None and sys.prefix == sys.base_prefix:
        raise RuntimeError(
            "Current python does not seem to be in a virtual environment, which is the "
            "recommended way to install dependencies for development. Please "
            "consider using an virtual environment for development."
        )
    context.run("python -m pip install --upgrade pip")
    command = "python -m pip install --upgrade -r"
    context.run(command + " requirements.txt")
    context.run(command + " requirements-dev.txt")
