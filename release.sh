#!/bin/bash

# Usage: ./release.sh [sync|test|doc|notes|web|pypi]
#
# Build an official elements package release
#
# Releasing this package requires some setup in your local environment
#    pytest
#    sphinx and latex
#    reflectometry.org server key for updating docs
#    ~/.pypirc should be defined
#

# Adapt the following to your own username/password for pypi, and get yourself
# added to the periodictable package.
# == ~/.pyirc ==
# [distutils]
# index-servers =
#     pypi
#
# [pypi]
# username:...
# password:...

# Set the python processor to use for the build tests
# Note: you may need to change this for your own environment
#PYTHON33=python3.3
#PYTHON27=python2.7
#PYTHON26=python2.6
PYTHON33=~/anaconda/envs/bumps3x/bin/python
PYTHON27=~/anaconda/envs/bumps/bin/python
PYTHON26=~/anaconda/envs/bumps26/bin/python
#BROWSER=firefox # linux
#PDFREADER=evince # linux
#RST2HTML=rst2html # linux
BROWSER=open # OSX
PDFREADER=open # OSX
RST2HTML=rst2html.py

case "$1" in
sync)   step=0;;
test)   step=1;;
doc)    step=2;;
notes)  step=3;;
web)    step=4;;
pypi)   step=5;;
*)      step=0;;
esac


function ready() {
  stepname=$1; shift
  echo -n "$* [y/n] "
  read ans && test "$ans" != "y" \
    && echo Restart with ./release.sh $stepname && exit
}

if [ $step -le 0 ]; then
  echo === Version control status ===
  set -x
  git pull
  git status
  set +x
  ready sync Is the repository up to date?
fi

if [ $step -le 1 ]; then
  echo === Tests ===
  set -x
  $PYTHON27 test.py -q --with-coverage
  $PYTHON26 test.py -q
  $PYTHON33 test.py -q
  set +x
  if false; then
    echo
    # Ask hudson build server if package is working on all platforms
    hudson_server="localhost:8080"
    elements_job="elements"
    url="http://$hudson_server/job/$elements_job/lastBuild"
    jsonurl="$url/api/json?depth=0"
    if curl --silent $jsonurl | grep -q SUCCESS ; then
        echo latest hudson build was successful
    else
        echo **** latest hudson build failed ... see $url
        $BROWSER $url &
    fi
  fi
  ready test Are the tests okay?
fi

if [ $step -le 2 ]; then
  echo === Documentation ===
  (cd doc/sphinx && make clean html pdf)
  $BROWSER doc/sphinx/_build/html/index.html >/dev/null 2>&1 & 
  $PDFREADER doc/sphinx/_build/latex/PeriodicTable.pdf >/dev/null 2>&1 &
  ready doc Does the documentation build correctly?
fi

if [ $step -le 3 ]; then
  echo === Release notes ===
  $RST2HTML README.rst > /tmp/README.html
  $BROWSER /tmp/README.html >/dev/null 2>&1 &
  git log --format="%Cred%ad%Creset %s %Cred%an%Creset" --date=short
  version=$(grep __version__ periodictable/__init__.py | sed -e's/^.*= *//;s/"//g')
  echo *** Version is $version
  ready notes Are the release notes shown in the browser up to date?
fi

if [ $step -le 4 ]; then
  ready web: Push docs to the web?
  ssh reflectometry.org rm -r web/danse/docs/elements
  find doc/sphinx/_build/html | xargs chmod ug+rw
  find doc/sphinx/_build/html -type d | xargs chmod g+x
  (cd doc/sphinx/_build && scp -r html reflectometry.org:web/danse/docs/elements)
  ready web Documentation upload successful?
fi

if [ $step -le 5 ]; then
  ready Push package to pypi?
  git tag -a v$version -m "Release $version"
  git push --tags
  python setup.py sdist
  #twine upload dist/periodictable-$version.tar.gz -p <token>
  twine upload dist/periodictable-$version.tar.gz
  ready pypi Package upload successful?
fi

echo == All done! ==
