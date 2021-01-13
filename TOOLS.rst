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
* **tox** is used as a tool to run all quality assurance (**QA**)
  tools and commands across all supported enviroments, typically
  before push or by continuous integration (CI) tools

Install instructions
--------------------

Install all development tools by executing::

  pip install -r requirements-dev.txt

while being at the git archive root in your virtual environment or
anaconda environment.

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
referred to as ``hooks``. For xidat a pre-push is recommended (as
opposed to a pre-commit hook, which take time to run the tests on
every commit). The pre-push hooks are located in ``tools/hooks``.

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
folder. (The interested reader will notice that the only difference
between the Linux and the Windows version of the hook, is the line at
the top that indicates the executable that should run the program).

Then there is a bit of difference depending on how git is used. If the
main interface to git is via powershell (or mayby Command Prompt) and
a virtual environment is being used and is active, then the commit
hook will just work.

If the development is done elsewhere, but still somehow relies on a
virtual environment (either a separate virtual environment or an
anaconda environment), and git is used via the Git Bash program that
is installed along with Git for Windows, then a bit more work is
required. The problem is that Git Bash does not know about the path of
the development tools, so that will have to be set manually.

First we should locate the path of the tools. In case a separate
virtual environment is used, located e.g. in ``c:\venv\ixdat``, then
``c:\venv\ixdat\Scripts`` will be the path of tools. TODO anaconda.

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

To run a command with an **invoke specific** argument do::

 $ invoke clean --dryrun

To get help on a command, e.g. ``clean``, do::

 $ invoke --help clean
