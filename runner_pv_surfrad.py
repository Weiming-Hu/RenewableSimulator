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

import argparse
import warnings

from Simulator import SimulatorSolarSurfrad

# Ignore warnings, in this case, coming from Numba
warnings.filterwarnings('ignore')


if __name__ == '__main__':

    scenario_message = "A dictionary defines key values for scenarios"
    welcome_msg = "***********************************************************************\n" + \
                  "Our pursuit of a more sustainable future with solar photovoltaic energy\n" + \
                  "***********************************************************************\n\n" + \
                  "This utility works with the SURFRAD dataset.\n"

    # Define arguments
    parser = argparse.ArgumentParser(description=welcome_msg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--root', help="Data folder", required=True)
    parser.add_argument('--output', help='Output NetCDF file name', required=True)
    parser.add_argument('--scenario', help=scenario_message, required=False, default="scenarios.yaml")
    parser.add_argument('--ext', help='Data file extension', required=False, default=".dat")
    parser.add_argument('--solar', help="Method for solar position calculation", required=False, default="nrel_numpy")
    parser.add_argument('--cores', help="How many cores to use for each process", required=False, default=1, type=int)
    parser.add_argument('--verbose', help="Print verbose messages", action='store_true')
    parser.add_argument('--no-bars', help="Disable progress bars", action='store_true', dest='no_bars')
    parser.add_argument('--no-profile', help="Turn off simple clock profiling", action='store_true', desk='no_profile')

    # Parse arguments
    args = parser.parse_args()

    if args.verbose:
        print(welcome_msg)
        print('Argument summary:')
        for arg in vars(args):
            print('{}: {}'.format(arg, getattr(args, arg)))
        print()
    
    # Create an object of the class Scenarios
    simulator = SimulatorSolarSurfrad(args.root, args.output, args.scenario, args.ext,
                                      args.solar, args.cores, args.verbose, args.no_bars)

    print('\n' + simulator.summary() + '\n')

    simulator.simulate()

    if not args.no_profile:
        print()
        print(simulator.timer)
