
"""Definition of invoke tasks"""

import configparser
from shutil import rmtree
from invoke import task
from pathlib import Path
from subprocess import check_call, CalledProcessError


THIS_DIR = Path(__file__).parent
SOURCE_DIR = THIS_DIR / "src" / "ixdat"
CLEAN_PATTERNS = ("__pycache__", "*.pyc", "*.pyo", ".mypy_cache")


tox_config = configparser.ConfigParser()
tox_config.read("tox.ini")


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
    environments_to_run = []

    # Check which pythons are available on the system

    # TODO this almost certainly is Linux and possibly MacOS specific
    # and will need a Windows port
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

    context.run("tox -e " + ",".join(environments_to_run))


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
