import argparse
import ssdb
from ssdb import Habit

#
# Parse arguments
#

parser = argparse.ArgumentParser(description='Download single scattering data.')
parser.add_argument('names', metavar='N', type=str, nargs='*',
                    help='Names of the habits to download.')
args = parser.parse_args()

for n in args.names:
    h = Habit(n)
    if h.local:
        print("Habit {0} found locally.".format(n))
    else:
        print("Downloading habit {0} ...".format(n))
        h.download()
