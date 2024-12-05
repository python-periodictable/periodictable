.. _contributing:

********************
Contributing Changes
********************

The best way to contribute to the periodic table package is to work
from a copy of the source tree in the revision control system.

The source is available via git at `<https://github.com/python-periodictable/periodictable>`_.
To make changes, create a fork of the project on github, then do following,
with your github user name substituted for *GITNAME*::

    git clone https://github.com/GITNAME/periodictable.git
    cd periodictable
    python setup.py develop

By using the *develop* keyword on setup.py, changes to the files in the
package are immediately available without the need to run setup.py
install each time.

Track updates to the original package using::

    git remote add upstream https://github.com/python-periodictable/periodictable.git
    git remote -v   # check that it is set
    git fetch upstream
    git merge upstream/master

Please update the documentation and add tests for your changes.  We use
doctests on all of our examples that we know our documentation is correct.
More thorough tests are found in test directory.  Using the the nose test
package, you can run both sets of tests::

    pip install pytest pytest-cov
    pytest

When all the tests run, create a pull request (PR) on the github page.

Windows user can use `TortoiseGit <http://code.google.com/p/tortoisegit/>`_
package which provides similar operations.

You can build the documentation as follows (linux, mac)::

    pip install sphinx
    (cd doc/sphinx && PYTHONPATH=../.. make clean html pdf)

You can see the result by pointing your browser to::

    periodictable/doc/sphinx/_build/html/index.html
    periodictable/doc/sphinx/_build/latex/PeriodicTable.pdf

ReStructured text format does not have a nice syntax for superscripts and
subscripts.  Units such as |g/cm^3| are entered using macros such as
\|g/cm^3| to hide the details.  The complete list of macros is available in

        doc/sphinx/rst_prolog

In addition to macros for units, we also define cdot, angstrom and degrees
unicode characters here.  The corresponding latex symbols are defined in
doc/sphinx/conf.py.
