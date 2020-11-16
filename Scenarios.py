# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
#
# Author:
#
#     Weiming Hu <weiming@psu.edu>
#
# Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
# Department of Geography and Institute for CyberScience
# The Pennsylvania State University
#

import os
import math
import yaml


class Scenarios(dict):
    """
    Scenarios is the class that stores and define multiple scenarios that will be simulated through the renewable
    simulator. It is inherited from the python built-in class, dict.

    The dictionary values only accept tuples or lists.

    A single scenario will be defined by a certain combination of values from all keys. For example, scenario #1 is
    defined using the first value in all keys; scenario #2 is defined using the second value in the first key and
    the first value in all other keys.
    """

    def total_scenarios(self):
        """
        Counts the total number of scenarios that could be defined using the current object of the class Scenarios.
        The total number is calculated from the multiplication of the length of each dictionary entry.

        :return: the number of total scenarios
        """
        total = 1

        for key in self:
            if isinstance(self[key], (list, tuple)):
                total *= len(self[key])
            else:
                raise Exception('Scenarios only allow lists or tuples as values')
        return total

    def get_scenario(self, index):

        # Check out-of-bound issues
        total_scenarios = self.total_scenarios()
        if index >= total_scenarios:
            raise Exception("Index {} is out of bound. The toal number of scenarios is {}".format(index, total_scenarios))

        # Initialize variables
        accumulate_len = 1
        current_scenario = {}

        for key, value in self.items():

            # Calculate the current index for this configuration
            if accumulate_len == 1:
                accumulate_len = len(value)
                current_index = index % accumulate_len
            else:
                current_index = math.floor(index / accumulate_len) % len(value)
                accumulate_len *= len(value)

            current_scenario[key] = value[current_index]

        return current_scenario

    def write_scenarios(self, file, overwrite=True):

        # Check output file
        if os.path.exists(file):
            if overwrite:
                os.remove(file)
            else:
                raise Exception("Output file {} exists and overwriting is not permitted".format(file))

        # Initialize a dictionary
        scenarios_dict = {}

        # Convert scenarios to dictionary with names
        for scenario_index in range(self.total_scenarios()):

            # Get the current scenario and the name
            scenario_name = 'scenario_{:05d}'.format(scenario_index)
            scenarios_dict[scenario_name] = self.get_scenario(scenario_index)

        # Write the current scenario
        with open(file, 'w') as f:
            yaml.dump(scenarios_dict, f)

    def __str__(self):
        msg = '******************** Scenarios ********************'
        msg += '\nTotal scenarios: {}'.format(self.total_scenarios())
        for k, v in self.items():
            msg += '\n  - {}: {} scenario(s)'.format(k, len(v))
        msg += '\n**************** End of Scenarios *****************'
        return msg
