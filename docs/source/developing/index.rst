.. _developing:

================
Developing ixdat
================

If there's an experimental technique or analysis procedure or database that ixdat
should support and doesn't, it might be because **you** haven't coded it yet.

Here are a few resources to help you get started developing ixdat.

Git and Github
**************

The source code for ixdat (and this documentation) lives at:
https://github.com/ixdat/ixdat

To develop ixdat, you will need to use git and github. This means

- `install git <https://git-scm.com/downloads>`_. Git bash works well for Windows.

- Create an account at https://github.com

- Make a fork of https://github.com/ixdat/ixdat so that you have your own version of the
  ixdat repository on your own github account.

- Clone your repository. Navigate in the terminal which you will use for git (e.g. git bash) to
  the location that you want the repository (e.g. ``/c/Users/your_user_name/git/``), and type::

    git clone https://github.com/ixdat/ixdat

- Install ixdat from the repository to use the ixdat code you're working on. You can do this in a virtual environment,
  but it is simpler to just install it dynamically. In your terminal or Anaconda Prompt, navigate
  to the folder which contains the ixdat project folder (e.g. ``/c/Users/your_user_name/git/``)
  and type::

    pip uninstall ixdat
    pip install -e ixdat

  If you want to go back to using the released version later, just re-install it from PyPi::

    pip install --upgrade ixdat

- Switch to the branch you want to work from. Note, this is also how to *use* a feature that is under development.::

    git switch branch_to_work_from

  Usually it makes most sense to start work by branching from the main branch::

    git switch main

- Make and switch to a branch on which to develop your feature::

    git switch -c my_feature_branch_name


- Develop your feature, committing regularly and pushing regularly to your github account.

- When it's ready (i.e., works like you want, and passes linting and testing), make a pull request!
  A pull request (PR) is an awesome open review process that gives others a chance to comment and suggest
  improvements in your code before it's merged into the main ixdat package. You can see
  existing pull requests at https://github.com/ixdat/ixdat/pulls


NEXT_CHANGES.rst
****************

When you contribute to ixdat, add descriptions of your changes to the file
**NEXT_CHANGES.rst** in the main project folder. This makes it easy for us to keep
track of everything that will be included in the next release of ixdat. When you make a
pull request, we will remind you to update NEXT_CHANGES.rst if you haven't done so.

When we release a version of ixdat, we do the following:

- Increment the version number according to the semantic versioning `conventions <https://semver.org>`_

- Copy-paste your descriptions and others from NEXT_CHANGES.rst to the main changelog, CHANGES.rst

- Build ixdat using setup.py

- Upload a copy to PyPI so that the new version installs with ``pip``

style
*****

We do our best to follow the conventions at

- code style guide: https://www.python.org/dev/peps/pep-0008/
- docstring style guide: https://www.python.org/dev/peps/pep-0257/

Exceptions include

- It's fine to capitalize names for a quantity that is conventionally capitalized in equations (`T` for temperature, for example).

The tools **black** and **flake8** help us keep the style up to standards.

tools
*****

We use tools to make sure that our code is both functional and pretty. This makes it
easier to work together. See instructions for the tools in `tools.rst <https://github.com/ixdat/ixdat/blob/main/TOOLS.rst>`_


Testing
*******

Software tests are welcomed! The basic test suite is included in this
repository in the ``tests`` directory.

However, if your tests have a specific purpose and/or depend on additional
data, then the data should be initialized in a **submodule**.  Additionally,
the tests should be **marked** as ``external`` with the pytest mark decorator.
It is done like this, so other developers do not have to download your data
through submodules, when they do not wish to test their new feature with your
tests and data.

If you want to run the complete test suite that includes external tests from
all contributing developers, you need to initialize necessary submodules::

    git submodule update --init --recursive

In order to download a specific submodule (e.g. ``ixdat-large-test-files``), the command is::

    git submodule update --init submodules/ixdat-large-test-files

And then run command for testing with ``--external`` option::

    # the 'tests' here is an invoke command
    invoke tests --external

    # or just with pytest, if you wish:
    # the 'tests' here is a project directory with the tests
    pytest tests --external

Tutorials
*********

The tutorials for ixdat is developed in a
`separate repository <https://github.com/ixdat/tutorials>`_. But these tutorials are copied into
the ixdat repository in order to be able to generate docs from the Jupyter notebooks for the documentation:
https://ixdat.readthedocs.io/

We are trying to figure out the best way to automate this process (see https://github.com/ixdat/ixdat/pull/133),
but for now just copy the updated files of the tutorials repository into the corresponding place in
docs/source/tutorials/tutorials_repo an compile them manually.

All files in docs/source/tutorials/tutorials_repo are ignored by the ixdat repo's .gitignore except
for the .ipynb files.

The tutorials should be reviewed and merged on the tutorials repo (in un-compiled state) before being
added to the tutorials page of the documentation.


Write to us
***********
We'd love to know what you're working on and help with any issues developing, even
before you make a PR.
One great way to do so is through `github discussions <https://github.com/ixdat/ixdat/discussions>`_
