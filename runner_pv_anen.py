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
import datetime
import warnings

from mpi4py import MPI
from Simulator import SimulatorSolarAnalogs
from Functions import get_start_index, get_end_index, get_nc_dim_size

# Ignore warnings, in this case, coming from Numba
warnings.filterwarnings('ignore')


if __name__ == '__main__':

    ###################
    # Parse arguments #
    ###################

    # Get MPI info
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_size = MPI.COMM_WORLD.Get_size()

    welcome_msg = "***********************************************************************\n" + \
                  "Our pursuit of a more sustainable future with solar photovoltaic energy\n" + \
                  "***********************************************************************\n\n" + \
                  "This utility works with Analog Ensemble output.\n"
    nc_message = "An NetCDF file with weather analogs generated by PAnEn (https://weiming-hu.github.io/AnalogsEnsemble)"
    map_message = "A variable map that translates variable names in the NetCDF file to internally required names"
    scenario_message = "A dictionary defines key values for scenarios"
    profiler_message = "One of [pyinstrument, yappi, line_profiler, memory, simple]"

    # Define arguments
    parser = argparse.ArgumentParser(description=welcome_msg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--nc', help=nc_message, required=True)
    parser.add_argument('--map', help=map_message, required=False, default="variable-map.yaml")
    parser.add_argument('--scenario', help=scenario_message, required=False, default="scenarios.yaml")
    parser.add_argument('--solar', help="Method for solar position calculation", required=False, default="nrel_numpy")
    parser.add_argument('--cores', help="How many cores to use for each process", required=False, default=1, type=int)
    parser.add_argument('--verbose', help="Print verbose messages", action='store_true')
    parser.add_argument('--no-bars', help="Disable progress bars", action='store_true', dest='no_bars')
    parser.add_argument('--profile', help="Turn on profiling", action='store_true')
    parser.add_argument('--profiler', default='simple', help=profiler_message)
    parser.add_argument('--stations-index', help='A list of station indices to simulate',
                        nargs='*', type=int, default=None, dest='stations_index')
    parser.add_argument('--re-simulate-sky-conditions', action='store_true', dest='re_simulate_sky_conditions',
                        help='Whether to force simulating sky conditions even when it exists')

    # Parse arguments
    args = parser.parse_args()

    if args.verbose and mpi_rank == 0:
        print(welcome_msg)
        print('Argument summary:')
        for arg in vars(args):
            print('{}: {}'.format(arg, getattr(args, arg)))
        print()

    # Start a profiler
    if args.profile and mpi_rank == 0:

        if args.profiler == "yappi":
            import yappi
            yappi.start()

        elif args.profiler == "pyinstrument":
            from pyinstrument import Profiler
            profiler = Profiler()
            profiler.start()

        elif args.profiler == "memory":
            from guppy import hpy

        elif args.profiler == "line_profiler":
            try:
                # Use the decorator from kernprof to profile the main simulation function
                SimulatorSolarAnalogs = profile(SimulatorSolarAnalogs)
            except NameError:
                raise Exception("Failed with function decorating. Did you properly use kernprof ?")

        elif args.profiler == "simple":
            pass

        else:
            raise Exception("Unsupported profiler: {}".format(args.profiler))

    ##############
    # Simulation #
    ##############
    if args.stations_index is None:
        num_stations = get_nc_dim_size(args.nc, 'num_stations', mpi_size > 1)
        args.stations_index = list(range(num_stations))
    else:
        num_stations = len(args.stations_index)

    start = get_start_index(num_stations, mpi_size, mpi_rank)
    end = get_end_index(num_stations, mpi_size, mpi_rank)

    if mpi_size > 1 and args.verbose:
        if mpi_rank == 0:
            print('{} workers (processes) have been initialized.'.format(mpi_size))

        print('Worker {}/{} processes stations index [{}, {}] out of {} indices in total.'.format(
              mpi_rank + 1, mpi_size, start, end, num_stations))

    simulator = SimulatorSolarAnalogs(args.nc, args.map, args.scenario, args.solar,
                                      mpi_size > 1, args.stations_index[start:end],
                                      not args.re_simulate_sky_conditions, args.cores,
                                      args.verbose if mpi_rank == 0 else False,
                                      args.no_bars if mpi_rank == 0 else True)

    if mpi_rank == 0:
        print('\n' + simulator.summary() + '\n')

    simulator.simulate()

    ############
    # Finalize #
    ############

    if args.profile and mpi_rank == 0:

        if args.profiler == "yappi":
            current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

            func_stats = yappi.get_func_stats()

            func_stats.save("yappi_{}_rank-{}.log".format(current_time, mpi_rank), type="callgrind")

            print("You have chosen yappi as the profiler tool.")
            print("yappi_* files are created. You can use callgrind tools for visualization")

        elif args.profiler == "pyinstrument":
            profiler.stop()
            print(profiler.output_text(unicode=True, color=True))

        elif args.profiler == 'memory':
            print("Master process heap usage:")
            heap_usage = hpy().heap()
            print(heap_usage)
            print()
            print(heap_usage.more)

        elif args.profiler == "simple":
            print()
            print(simulator.timer)
