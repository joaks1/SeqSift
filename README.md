Overview
========

SeqSift is a python package for vetting, manipulating, and converting molecular
sequence data. 

The package is written by Jamie Oaks of the University of Kansas Department of
Ecology and Evolutionary Biology.

Requirements
============

SeqSift has only been tested on version 2.7 of Python.

It requires the Biopython: http://biopython.org/wiki/Biopython

Some components of Biopython require NumPy (http://numpy.scipy.org). However,
none of the Biopython code currently used by SeqSift requires NumPy. Thus,
SeqSift will work if you intall Biopython without NumPy.

Installation
============

You should have Python 2.7 and Biopython installed.

Clone Git repository and install using:

    $ git clone git://github.com/joaks1/SeqSift.git
    $ cd SeqSift
    $ python setup.py install

If the install fails due to a permission error, try:

    $ sudo python setup.py install

To run the test suite, use:

    $ python setup.py test

If you plan to develop the code, install via:

    $ python setup.py develop

Citing SeqSift
==============

Please cite the Biopython library upon which SeqSift depends.

Acknowledgements
================

This software greatly benefited from funding provided to Jamie Oaks from the
National Science Foundation (DEB 1011423 and DBI 1308885), University of Kansas
(KU) Office of Graduate Studies, Society of Systematic Biologists, Sigma Xi
Scientific Research Society, KU Ecology and Evolutionary Biology Department,
and the KU Biodiversity Institute.

License
=======

This program is distributed under the GPL license in the hope that it will be
useful, but WITHOUT ANY WARRANTY.  See the GNU General Public License for more
details.

