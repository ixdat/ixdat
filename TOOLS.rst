This file explains the development tool chain around ixdat

Tools used
==========

The following is a list of specific tools used during ixdat
development:

* **pytest** is used for running tests
* **flake8** is used for linting
* **sphinx** is used to build documentation

The following is a list of "tool and commands runners":
   
* **invoke** is used during development, to run tools and other
  pre-configured maintenance tasks inside the existing development
  environment
* **tox** is used as a tool to run all **QA** tools and commands
  across all supported enviroments, typically before push or by
  continuous integration (CI) tools

Tool and command runners in more detail
=======================================

The reason there are two different tool runners is that they have very
different purposes:

**invoke** is used as general purpose, platform independent way of
defining and running "tasks relevant for development" (reminiscent of
the *makefiles* of olden times, commonly used for software development
on Unix systems). These tasks can be something like running a linter,
type checker, test runner or code formatter, building documentation,
cleaning the repository of temporary files, building a package or
updating the list of contributors from the git log.

**tox** is tuned towards a specific purpose, which is to run
specifically QA tools (something like linter, type checker, test
runner or code formatter) in a clean environment, with a fresh install
of code, across as many supported python versions as possible. `tox`
can also be used in CI tools, such a Github Actions, to specify what
to run as part of the CI, which eliminates the need to an extra
configuration.

Both tools share the purpose of hiding the details of the way the
tools should be run (arguments, settings etc.) and defining a single
place to define those.

From the definitions above, however, it is obvious that the QA tools
may be run by both. This is ok, because the purposes are different.
But it raises the important concern of *where* the one true place to
define the arguments and settings are for those. See definitions of
that in the following section.

Where are settings and arguments kept
-------------------------------------

For settings and arguments for tools the following rules apply:

* **setup.cfg** If the tool support defining them in the shared
  configuration file `setup.cfg`, then they belong there no-matter how
  they are run. If not (settings may only be available at the command
  line), see the next points.
* **tox.ini** If the tool is run by both `tox` and `invoke`, then the
  settings belong in the `tox` configuration file and `invoke` will
  have to read them from there
* **tasks.py** If the tool is `invoke` specific, then obviously the
  configuration goes into its "configuration file" `tasks.py`

Non-tool related "settings"
```````````````````````````

**Package metadata** goes into src/ixdat/__init__.py and README.rst
and all remaining information about how to create a package goes into
setup.py

**requirements** goes into the two requirements files;
`requirements.txt` for package requirements and `requirements-dev.txt`
for development requirements.

Git hooks
=========

The version control system `git` has the option of having check done
before certain important actions, like a commit or a push. These are
referred to as `hooks`. For xidat a pre-push is recommended (as
opposed to a pre-commit hook, which take time to run the tests on
every commit). The pre-push hook is located in `tools/hooks`.

Git hooks cannot be distributed automatically, so it will need to be
"installed". On Linux systems this is done by sym-linking it into the
proper folder:

$ cd .git/hooks
$ ln -s ../../tools/hooks/pre-push .

Windows instructions TODO.

Command quick tips
==================

tox
---

To run all `test environments` simply run

$ tox

To pick a specific one to run, use the `-e` flag followed by the name:

$ tox -e flake8

To list the enviroments do:

$ tox -l

To force recreation of all tox' virtual environments do

$ tox -r

invoke
------

To see a list of all tasks:

$ invoke --list

To run a specific task, say linting, run:

$ invoke lint

To run a command with an **invoke specific** argument do:

$ invoke clean --dryrun

To get help on a command, e.g. `clean`, do:

$ invoke --help clean
