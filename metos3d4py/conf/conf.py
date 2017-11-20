#
# Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
# Copyright (C) 2017  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sys
import yaml
from metos3d4py.util.util import _print_message, _print_error, _print_usage

class Conf:
    """
        An instance of the Conf class is used to read in YAML configuration file.
        
        Attributes:
            dict        # the dictionary where the parsed YAML file is stored

        """
    def init(self, comm, argv):
        # the first command line argument is always the name of the current executable,
        # expect the file path as second argument,
        if len(argv) > 1:

            _print_message(comm, "parsing configuration file: {}".format(argv[1]))

            try:
                # open file
                f = open(argv[1])
            except Exception as e:
                _print_error(comm, "Cannot open file: {}".format(argv[1]))
                _print_message(comm, e)
                _print_usage(comm)
                sys.exit(1)

            try:
                # parse yaml content
                self.dict = yaml.load(f)
            except Exception as e:
                _print_error(comm, "Cannot parse as YAML file: {}".format(argv[1]))
                _print_message(comm, e)
                _print_usage(comm)
                sys.exit(1)

            f.close()

        else:
            _print_error(comm, "No configuration file given.")
            _print_usage(comm)
            sys.exit(0)

    def __str__(self):
        return "Conf:\n  Keys: {}".format(list(self.dict.keys()))


