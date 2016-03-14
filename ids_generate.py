import os, itertools, math, collections
from utils import *


indir = "."
outdir = "."


filename = "data_exp_tracked.txt"
outfile = "data_exp_tracked_ids.txt"

L = 512

"""
0: positions x
0: angles x
1: positions
1: angles
_
1 x
1 x
2
2
_
"""

"""
after first compare-loop:
0: [0,1,2,3,4]
1: [8,2,6,3,5,7]
1: [-1,4,2,-1,-1,-1,-1]
2: 
"""

with open(os.path.join(indir, filename), "r") as fr, open(os.path.join(outdir, outfile), "w") as fw:
    fw.write("W: {} H: {}\n".format(L, L))
    for ids, __ in tracking_ids(fr, 0, 4000):
        fw.write(" ".join((str(id) for id in ids)))
        fw.write("\n")

