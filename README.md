# Amber Parser
## Version 0.0.1a
This is the parser for SANDER and PMEMD Molecular Dynamics codes in [AMBER](http://ambermd.org).
The official version lives at:

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-amber.git

You can browse it at:

    https://gitlab.rzg.mpg.de/nomad-lab/parser-amber

It relies on having the nomad-meta-info and the python-common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

This parser will be in the directory parsers/amber of this repository.

# Running and Testing the Parser
## Requirements
The required python packages can be installed with (see [python-common](https://gitlab.rzg.mpg.de/nomad-lab/python-common)):

    pip install -r nomad-lab-base/python-common/requirements.txt

## Usage
Amber (SANDER/PMEMD) output files can be parsed with:

    python AMBERParser.py [path/toFile]

## Test Files
Example output files of Amber (SANDER/PMEMD) can be found in the directory test/examples.
More details about the calculations and files will be explained in this README file.

# Documentation

## Installation of NOMAD Infrastructure
More information on the installation of NOMAD infrastructure can be found [here](https://workflowy.com/s/DKhN.2pfL6VMJzn)
