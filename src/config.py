"""Simple config file

"""

import os
import sys

class configuration(object):
    pass

config = configuration()
config.base = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
config.A549_test = os.path.join(config.base, "data", "test_data", "input_data", 'A549')
