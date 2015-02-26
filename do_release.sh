# Release script
# Copyright (C) 2013-2015 Mitchell Jon Stanton-Cook
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# m.stantoncook@gmail.com
# School of Chemistry & Molecular Biosciences
# The University of Queensland
# Brisbane, QLD 4072.
# Australia


#VERSION=1.0.3

# Perform an install-uninstall cycle
pip uninstall Contiguity 
python setup.py install
pip uninstall Contiguity
python setup.py clean


# Do all the versioning stuff here..
bumpversion patch


# Clean, test, build the source distribution & pip install it
# Need to get exit statuses here...
python setup.py clean
#python setup.py test
#STATUS=`echo $?`
#if [ $STATUS -eq 0 ]; then
#    echo ""
#else
#    echo "Tests failed. Will not release"
#    exit
#fi 

python setup.py sdist bdist_wheel
pip install dist/Contiguity-$VERSION.tar.gz
STATUS=`echo $?`
if [ $STATUS -eq 0 ]; then
    echo ""
else
    echo "Package is not pip installable. Will not release"
    exit
fi 


# Docs
# Need to get exit statuses here...
cd docs
make clean
sphinx-apidoc -o API ../Contiguity
mv API/* .
rmdir API
make html
cd ..

git push
# tag & push the tag to github
GIT=`git status`
CLEAN='# On branch master nothing to commit, working directory clean'
if [ "$s1" == "$s2" ]; then
    git tag v$VERSION
    git push --tags
else
    echo "Git not clean. Will not release"
    exit
fi 


# Upload to PyPI & clean
twine upload -u mscook -p $PYPIPASS dist/* && python setup.py clean
