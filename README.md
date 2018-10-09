pipelign
========

[![PyPI](https://img.shields.io/pypi/v/pipelign.svg)](https://pypi.python.org/pypi/pipelign)
[![GitHub license](https://img.shields.io/github/license/asmmhossain/pipelign.svg)](./LICENSE.txt)
[![Requires.io](https://img.shields.io/requires/github/asmmhossain/pipelign.svg)](https://requires.io/github/asmmhossain/pipelign/requirements/)
[![Travis](https://img.shields.io/travis/asmmhossain/pipelign.svg)](https://travis-ci.org/asmmhossain/pipelign)

\$

:   Built from
    [makenew/python-package](https://github.com/makenew/python-package).

\$

Description
-----------

\$ A pipeline for automated alignment\$

Installation with pip
---------------------

Install it directly using pip with

    $ pip install pipelign

Installation with setuptools
----------------------------

    $ python3 setup.py install

Development and Testing
-----------------------

### Source Code

The [pipelign source](https://github.com/asmmhossain/pipelign) is hosted
on GitHub. Clone the project with

    $ git clone https://github.com/asmmhossain/pipelign.git

### Requirements

You will need [Python 3](https://www.python.org/) with
[pip](https://pip.pypa.io/).

Install the development dependencies with

    $ pip install -r requirements.devel.txt

### Building a conda package

Installation with conda
-----------------------

First create an environment.

    $ conda create -n pipelign python=3

Activate the environment.

    $ source activate pipelign

Install conda-build:

    $ conda install conda-build

Run conda.

    $ conda-build . -c bioconda

### Tests

Lint code with

    $ python setup.py lint

Run tests with

    $ python setup.py test

Contributing
------------

Please submit and comment on bug reports and feature requests.

To submit a patch:

1.  Fork it (<https://github.com/asmmhossain/pipelign/fork>).
2.  Create your feature branch (`git checkout -b my-new-feature`).
3.  Make changes. Write and run tests.
4.  Commit your changes (`git commit -am 'Add some feature'`).
5.  Push to the branch (`git push origin my-new-feature`).
6.  Create a new Pull Request.

License
-------

This Python package is licensed under the MIT license.

Warranty
--------

This software is provided \"as is\" and without any express or implied
warranties, including, without limitation, the implied warranties of
merchantibility and fitness for a particular purpose.
