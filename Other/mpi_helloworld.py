#!/home/mlam/anaconda/bin/
"""
Parallel Hello World
"""

from   mpi4py     import MPI
import sys

comm = MPI.COMM_WORLD
size = comm.size
rank = comm.Get_rank()
name = MPI.Get_processor_name()

sys.stdout.write("Hello, World! I am process %d of %d on %s.\n" % (rank, size, name))

