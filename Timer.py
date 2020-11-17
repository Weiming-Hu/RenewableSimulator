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
# A simple timer for profiling
#

import numpy as np

from time import time


class Timer(list):

    # The timer log consists of a list of 3 values, namely the name, the start time, and the end time
    _GROUP_LENGTH = 3

    def start(self, tag):
        self.extend((tag, time()))

    def stop(self):
        self.append(time())

    def check_sanity(self):
        assert len(self) % Timer._GROUP_LENGTH == 0, 'Timer has unexpected logs'
        for index in range(0, len(self), Timer._GROUP_LENGTH):
            assert isinstance(self[index], str), 'Timer has unexpected logs'
            assert isinstance(self[index + 1], float), 'Timer has unexpected logs'
            assert isinstance(self[index + 2], float), 'Timer has unexpected logs'

    def __str__(self):
        self.check_sanity()

        msg = '*************** Simple Clock ***************\n'

        session_names = []
        session_costs = []

        for index in range(0, len(self), Timer._GROUP_LENGTH):
            session_names.append(self[index])
            session_costs.append(self[index + 2] - self[index + 1])

        total_cost = np.sum(session_costs)
        session_portions = np.array(session_costs) / total_cost * 100
        session_portions = session_portions.tolist()

        msg += 'Total wall time: {:.2f} s\n'.format(total_cost)

        for index in range(len(session_names)):
            msg += '{}: {:.2f} s ({:.2f}%)\n'.format(
                session_names[index], session_costs[index], session_portions[index])

        msg += '*********** End of Simple Clock ************'

        return msg
