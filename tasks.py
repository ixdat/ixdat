
"""Definition of invoke tasks"""

import sys
import configparser
import platform
from shutil import rmtree
from invoke import task
from pathlib import Path
from subprocess import check_call, CalledProcessError, check_output


THIS_DIR = Path(__file__).parent
SOURCE_DIR = THIS_DIR / "src" / "ixdat"
CLEAN_PATTERNS = ("__pycache__", "*.pyc", "*.pyo", ".mypy_cache")


tox_config = configparser.ConfigParser()
tox_config.read(THIS_DIR / "tox.ini")


def extract_first_command_from_tox_env(environment):
    """Extract the first command from list of commands from tox.ini"""
    commands_string = tox_config[environment]["commands"]
    commands = commands_string.split("\n")
    if commands[0].strip() == "" and len(commands) > 0:
        commands = commands[1:]
    return commands[0]


# ### QA tasks


@task(aliases=["lint"])
def flake8(context):
    """Run the flake8 task"""
    command = extract_first_command_from_tox_env("testenv:flake8")
    print("# flake8")
    return context.run(command).return_code


@task(aliases=["test", "tests"])
def pytest(context):
    """Run the pytest task"""
    command = extract_first_command_from_tox_env("testenv")
    print("# pytest")
    return context.run(command).return_code


@task(aliases=["QA", "qa", "check"])
def checks(context):
    """Run all QA checks"""
    combined_return_code = flake8(context)
    combined_return_code += pytest(context)
    if combined_return_code == 0:
        print()
        print(r"+----------+")
        print(r"| All good |")
        print(r"+----------+")


@task
def tox(context):
    """Run tox for the python interpreters available on your system"""
    environments = tox_config["tox"]["envlist"].split(", ")

    # Check which pythons are available on the system
    if platform.system() == "Windows":
        environments_to_run = filter_tox_environments_windows(environments)
    else:
        environments_to_run = filter_tox_environments_linux(environments)

    context.run("tox -p auto -e " + ",".join(environments_to_run))


def filter_tox_environments_linux(environments):
    """Filter tox environments to only those available on this Linux system"""
    environments_to_run = []
    for environment in environments:
        if environment.startswith("py"):
            if environment.startswith("pypy"):
                command = environment
            else:
                # The environments look like: py36
                command = "python{}.{}".format(*environment[2: 4])

            # Check that the executable exists
            try:
                check_call([command, "--version"])
            except (CalledProcessError, FileNotFoundError):
                continue

            # Certain version of Ubuntu may have a old "reduced"
            # Python version, with an "m" suffix. It is unsifficient
            # for tox so check that it isn't there.
            try:
                check_call([command + "m", "--version"])
                continue
            except (CalledProcessError, FileNotFoundError):
                pass

            environments_to_run.append(environment)
        else:
            environments_to_run.append(environment)
    return environments_to_run


def filter_tox_environments_windows(environments):
    """Filter tox environments to only those available on this Linux system"""
    python_versions = []
    for path in sys.path:
        ppath = Path(path) / "python.exe"
        if ppath.is_file():
            version_string = check_output([str(ppath), "--version"]).decode("ascii").strip()
            version_numbers = version_string.replace("Python ", "").split(".")
            version = "".join(version_numbers[:2])
            python_versions.append("py" + version)

    environments_to_run = []
    for environment in environments:
        if environment.startswith("py"):
            if environment in python_versions:
                environments_to_run.append(environment)
        else:
            environments_to_run.append(environment)

    return environments_to_run


# ### Maintenance tasks


@task
def clean(c, dryrun=False):
    """Clean the repository of temporary files and caches"""
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
