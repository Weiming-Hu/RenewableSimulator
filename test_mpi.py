from mpi4py import MPI
from netCDF4 import Dataset
from Functions import simulate_sun_positions, simulate_power, get_start_index, get_end_index


if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()

    nc = Dataset("/glade/scratch/wuh20/data/AnEn/analogs_domain-USA_chunk-01.nc", "r", parallel=True)
    num_stations = nc.dimensions['num_stations'].size

    start = get_start_index(num_stations, num_procs, rank)
    end = get_end_index(num_stations, num_procs, rank)

    print("Rank #{} starts reading from {} to {} ...".format(rank, start, end))
    nc_ghi = nc.variables["Albedo"][:, :, :, start:end]
    print("Rank #{} finished")
    nc.close()

