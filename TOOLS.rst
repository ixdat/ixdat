This file explains the development tool chain around ixdat

Tools used
==========

The following is a list of specific tools used during ixdat
development:

* **black** is used for formatting code
* **pytest** is used for running tests
* **flake8** is used for linting
* **sphinx** is used to build documentation

The following is a list of "tool and commands runners":

* **invoke** is used during development, to run tools and other
  pre-configured maintenance tasks inside the existing development
  environment
* **tox** is used as a tool to run all quality assurance (**QA**)
  tools and commands across all supported enviroments, typically
  before push or by continuous integration (CI) tools

Install instructions
====================

To setup a development environment (virtual environment or anaconda
environment) the first time, run the following 2 commands within the
active development environment::

  python -m pip install invoke
  invoke depencencies 

After the initial setup above, the development environment can at any
time be kept up to date with changes to tooling by re-running the
command::

  invoke dependencies


Tool and command runners in more detail
=======================================
Each of these tools and tool runners can be installed with
`$ pip install <tool_name.`,
or all at once with `$ pip install -r requirements-dev.txt` as
recommended above.

For each of the tools, we refer to its own documentation but offer a quick
summary here with our suggestions for usage.

Tools
-----
black
.....
**black** is an autoformatter, which fixes your white space usage etc.
https://black.readthedocs.io/en/stable/
It's nice to run it from the same terminal
that you use to `git commit`, so that you can black-format your code right
before committing (and avoid a "fix formatting" commit later). To get black
and the other tools available in git bash, you have to tell ~\.bashrc where
it is (see Windows instructions below under git hooks).

flake8
......
**flake8** is a linter, which checks the code for errors.
https://flake8.pycqa.org/en/latest/
This includes
syntax errors, but also programming errors like using a variable which
has not been defined, and style errors. flake8 can be a bit over-zeleous,
so even perfectly formatted (i.e. black-formatted) code can sometimes
trigger it. An example is unused import statements in __init__.py, which
are fine practice to make it quicker to import the most important parts of
a package. To allow these, we have added to the project's setup.cfg, the line
`--per-file-ignores = src/*/__init__.py:F401`
Additional allowances may need to be added there.
flake8 also enforces a maximum line length of 89, chosen to match the
default setting of black (+/- 1 char).

pytest
......
**pytest** is a suite of stuff used to write and run software tests.
https://docs.pytest.org/en/stable/

sphinx
......
**sphinx** is used to generate the beautiful documentation on
https://ixdat.readthedocs.io from ReStructuredText and ixdat source code.
To set it up just install sphinx, if you haven't already. In your terminal or Anaconda prompt, type::

  $ pip install sphinx

Then, to build the documentation, navigate to ``ixdat/docs``, and run in your terminal or Anaconda prompt::

  $ ./make html

Note, if you get an "access denied" error, you will just need to run the terminal as an administrator.

Then you can see the built documentation in your browser by double-clicking
``ixdat/docs/build/html/index.html``


**sphinx** is a tool for building the documentation into html and other
formats from restructured text (.rst). It also enables automatic documentation
generation from the doc-strings in the code. Read The Docs uses sphinx to compile
the documentation (ixdat/docs/source/index.rst) to make the will-be-beautiful
ixdat documentation at https://ixdat.readthedocs.io/

Tool runners
------------

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
of code, across as many supported python versions as possible. ``tox``
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

* **setup.cfg**) (not to be confused with setup.py, which is used for
  building the package) If the tools support configuration in the
  shared configuration file ``setup.cfg``, then it belong there
  always!
* **pyproject.toml**) This file is another common configuration file
  much like setup.cfg, which is used by black. Unfortunately, black
  and flake8 (which has settings in setup.cfg) do not share a single
  shared settings file
* **tox.ini** and **tasks.py**) If the tools only support
  configuration on the command line, then it goes into the
  configuration files of the tool runners, ``tox.ini`` for tox and
  ``tasks.py`` for invoke.

Non-tool related "settings"
```````````````````````````

**Package metadata** goes into src/ixdat/__init__.py and README.rst
and all remaining information about how to create a package goes into
setup.py

**requirements** goes into the two requirements files;
``requirements.txt`` for package requirements and
``requirements-dev.txt`` for development requirements.

Git hooks
=========

The version control system ``git`` has the option of having check done
before certain important actions, like a commit or a push. These are
referred to as ``hooks``. For ixdat, a pre-push is recommended (as
opposed to a pre-commit hook, which take time to run the tests on
every commit). The pre-push hooks are located in ``tools/hooks``.

While the pre-push hook is optional for developers, we will require
that any pull request to a branch supported for users passes linting
and pytesting. If you are using a pre-push hook, you enforce with
yourself that your work passes these checks the whole time, avoiding
extra work later.

Git hooks cannot be distributed automatically, so it will need to be
"installed".

Linux instructions
------------------

On Linux systems this is done by sym-linking it into the proper
folder. Starting at the base folder of the git archive it looks like
this::

 $ cd .git/hooks
 $ ln -s ../../tools/hooks/pre-push .

Windows instructions
--------------------

The pre-push hooks on Windows is dependent on a full Git installation,
because it depends on the shell that is embedded in Git Bash. So if
there isn't already a full Git for Windows installed, do that
(https://git-scm.com/downloads) and afterwards confirm that the file
``C:/Program\ Files/Git/usr/bin/sh.exe`` exists.

The next step is to copy the Windows specific hook into place, so
locate the folder ``.git/hooks`` in the main folder of the git archive
and copy the file ``tools/hooks/pre-push_windows_git_bash`` into that
folder and shorten its name to ``.git/hooks/pre-push``. (The interested
reader will notice that the only difference between the Linux and the
Windows version of the hook, is the line at the top that indicates the
executable that should run the program).

Then there is a bit of difference depending on how git is used. If the
main interface to git is via powershell (or maybe Command Prompt) and
a virtual environment is being used and is active, then the commit
hook will just work.

If the development is done elsewhere, but still somehow relies on a
virtual environment (either a separate virtual environment or an
anaconda environment), and git is used via the Git Bash program that
is installed along with Git for Windows, then a bit more work is
required. This can be buggy.

The problem is that Git Bash does not know about the path of
the development tools, so that will have to be set manually.

First we should locate the path of the tools. In case a separate
virtual environment is used, located e.g. in ``c:\venv\ixdat``, then
``c:\venv\ixdat\Scripts`` will be the path of tools. If you are
using conda, look below.

Having found the path of the tools it needs to be added to the Bash
configuration file like so::

 $ cd
 $ pwd
 /c/Users/Kenneth Nielsen
 $ nano .bashrc

The last line will open the ``.bashrc`` file in the terminal editor
``nano`` (https://www.nano-editor.org/dist/latest/cheatsheet.html). In
the editor, add this line to end of the file::

  export PATH=/c/venv/ixdat/Scripts:$PATH

and save the file and exit by pressing ``Ctrl-x`` followed by ``Y``
(if your editor is in english) and ``Enter``. The
``/c/venv/ixdat/Scripts`` part of that line is the path located
before, but converted for Git Bash notation, where the C-drive is
called ``/c`` and ``/`` is used for directory separation.

After having done this all the development tools and the git hook
should work. You can test the development tools by starting a new Git
Bash shell, navigation to the git archive and executing the command::

 $ invoke tox
 GLOB sdist-make: C:\venv\ixdat\ixdat\setup.py
 ___________________________________ summary ___________________________________
   py39: commands succeeded
   flake8: commands succeeded
   congratulations :)

Anaconda Instructions
---------------------

If you are using Conda (Anaconda/Minicond) with Windows, you will need to
add two directories either to your system PATH variable (se explanation
`here <https://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them>`_
) or to your ``.bash_profile`` or ``.bashrc`` (see explanation
`here <https://stackoverflow.com/questions/6883760/git-for-windows-bashrc-or-equivalent-configuration-files-for-git-bash-shell>`_
).

For a full Anaconda installation with the default location for User "scott",
the two paths to add are::

  - ``C:\Users\scott\Anaconda3\Scripts``
  - ``C:\Users\scott\Anaconda3``

And the two corresponding lines to add to ``~/.bash_profile``, using any text editor, are::

  export PATH="$PATH:/c/ProgramData/Anaconda3/Scripts/"
  export PATH="$PATH:/c/ProgramData/Anaconda3/"

Note that ``C:\`` becomes ``/c/``.

After adding these to your PATH variable, they should appear in $PATH in git bash.
To check, type::

  echo $PATH

in git bash and see if they appear in the output. Once that is there you should
be able to run tox and the other tools and tool runners from git bash just by e.g.::

  $ tox

, and it should work in a pre-push hook.

NOTE: There can still be bugs. See https://github.com/ixdat/ixdat/issues/10


Command quick tips
==================

tox
---

To run all ``test environments`` simply run::

 $ tox

To pick a specific one to run, use the ``-e`` flag followed by the
name::

 $ tox -e flake8

To list the enviroments do::

 $ tox -l

To force recreation of all tox' virtual environments do::

 $ tox -r

invoke
------

To see a list of all tasks::

 $ invoke --list

To run a specific task, say linting, run::

 $ invoke lint

To run all (quality assurance) checks, run::

 $ invoke checks

To update the development environment with correct version of tools
and settings, run::

 $ invoke dependencies

To run a command with an **invoke specific** argument do::

 $ invoke clean --dryrun

To get help on a command, e.g. ``clean``, do::

 $ invoke --help clean
