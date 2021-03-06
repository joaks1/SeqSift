#! /usr/bin/env python

import sys
import os

__project__ = "SeqSift"
__version__ = "0.1.2"
__author__ = "Jamie Oaks"
__copyright__ = "Copyright 2013 Jamie Oaks."
__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""

from seqsift.utils import dataio
from seqsift.utils import alphabets
from seqsift.utils import stats
from seqsift.utils import functions

def _get_git_data(repo_path):
    try:
        import subprocess
        import datetime

        p = subprocess.Popen(
                ["git", "rev-parse", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        commit = stdout.strip()[0:12]

        p = subprocess.Popen(
                ["git", "name-rev", "--name-only", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        branch = stdout.strip()

        p = subprocess.Popen(
                ["git", "show", "--quiet", "--pretty=format:'%at'", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        t = stdout.strip().replace("'", "").replace('"', '')
        commit_time = datetime.datetime.fromtimestamp(float(t))

        return branch, commit, commit_time
    except:
        return None, None, None

__homedir__ = None
try:
    try:
        __homedir__ = __path__[0]
    except AttributeError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except IndexError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
except:
    pass

__gitinfo__ = ""

__branch__, __commit__, __committime__ = _get_git_data(__homedir__)

def get_description():
    d = "{0} version {1}".format(__project__, __version__)
    if __branch__:
        d += " {0}".format(__branch__)
    if __commit__:
        d += " {0}".format(__commit__)
    if __committime__:
        d += " {0}".format(__committime__)
    return d
