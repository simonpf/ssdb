import argparse
import ssdb
from ssdb import Habit, get_standard_habits

#
# Parse arguments
#

parser = argparse.ArgumentParser(description='Download single scattering data.')
parser.add_argument('names', metavar='N', type=str, nargs='*',
                    help='Names of the habits to download.')

args = parser.parse_args()
names = args.names

if len(names) == 1 and names[0] in ["all", "All", "ALL"]:
    names = get_standard_habits()

for n in names:
    h = Habit(n)
    if h.local:
        print("Habit {0} found locally.".format(n))
    else:
        print("Downloading habit {0} ...".format(n))
        h.download()
